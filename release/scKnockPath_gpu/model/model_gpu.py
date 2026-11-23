"""GPU (PyTorch) version of scKnockPath.

`scKnockPathGPU` subclasses `scKnockPath` and replaces the three costly steps with
PyTorch code that runs on a CUDA device:

1. pathway unfolding (`duplicate_genes`): one gather with a precomputed column index;
2. Gaussian model-X knockoffs (`knockoff_sampler`): Ledoit-Wolf covariance, the SDP for
   the S-matrix (log-barrier interior-point Newton method instead of DSDP/cvxpy), and
   the knockoff draw;
3. logistic sparse-group-lasso fit with K-fold CV over the lambda path (`fit`): all
   (lambda, fold) problems are solved at the same time as one batched problem.

The solver is the one groupyr's `LogisticSGL` uses (copt's non-accelerated proximal gradient
with backtracking, max_iter=1000, tol=1e-7, intercept initialised at 0.5), but by default
with the correct proximal operator of the sparse group lasso penalty

    alpha * [ l1_ratio * ||beta||_1 + (1 - l1_ratio) * sum_g sqrt(|g|) * ||beta_g||_2 ].

groupyr 0.3.4 (`SparseGroupL1.prox`) deviates from it in two ways:
  * the group shrinkage is applied only when every group survives the group threshold;
    otherwise the surviving groups are only l1-thresholded (hard group gating), so the
    iteration does not minimise the objective above;
  * the feature -> group map assumes that each group is a contiguous block in feature
    order, which is wrong for `merge_groups=True` (groups [g_i, g_i + L]).
`exact_prox=False` reproduces the first behaviour (and therefore the published CPU results
of the default penalty) to machine precision; groups are always indexed correctly.
`exact_prox='hard'` is the same operator written as a deliberate algorithm (group
hard-thresholding: l1 soft-threshold, then a group is either kept unchanged or zeroed);
it differs from groupyr only in iterations where every group survives.

The knockoff noise is drawn with `np.random.randn`, so `np.random.seed(...)` keeps the same
meaning as in the CPU version. W statistics, selection and saving are inherited unchanged.
"""
import sys
import time
import warnings

import numpy as np
import torch
from sklearn.metrics import log_loss
from sklearn.model_selection import KFold

from model_lambda import scKnockPath


def _to_t(a, device, dtype):
    if isinstance(a, torch.Tensor):
        return a.to(device=device, dtype=dtype)
    return torch.tensor(np.asarray(a), device=device, dtype=dtype)


def _sync(device):
    if torch.device(device).type == 'cuda':
        torch.cuda.synchronize(device)


# ----------------------------------------------------------------------------------------
# preprocessing / covariance
# ----------------------------------------------------------------------------------------
def standardize_t(X):
    """sklearn StandardScaler (ddof=0; near-constant columns are only centred)."""
    mean = X.mean(0)
    std = X.std(0, unbiased=False)
    eps = torch.finfo(X.dtype).eps
    std = torch.where(std < 10 * eps, torch.ones_like(std), std)
    return (X - mean) / std


def ledoit_wolf_t(X):
    """sklearn.covariance.LedoitWolf().fit(X): returns (covariance, precision, shrinkage)."""
    n, p = X.shape
    Xc = X - X.mean(0)
    emp = Xc.T @ Xc / n
    X2 = Xc * Xc
    tr = X2.sum(0) / n
    mu = tr.sum() / p
    beta_ = (X2.T @ X2).sum()
    delta_ = (emp * emp).sum()
    beta = (beta_ / n - delta_) / (p * n)
    delta = (delta_ - 2.0 * mu * tr.sum() + p * mu ** 2) / p
    beta = torch.minimum(beta, delta)
    shrink = torch.zeros((), dtype=X.dtype, device=X.device) if beta == 0 else beta / delta
    cov = (1.0 - shrink) * emp
    cov.diagonal().add_(shrink * mu)
    prec = torch.cholesky_inverse(torch.linalg.cholesky(cov))
    return cov, prec, float(shrink)


def _shift_until_psd(M, tol):
    mineig = torch.linalg.eigvalsh(M)[0]
    if mineig < tol:
        M = M + (tol - mineig) * torch.eye(M.shape[0], dtype=M.dtype, device=M.device)
    return M


def _chol_reg(H):
    """Cholesky of a PD matrix that may be numerically indefinite: add a growing multiple of
    its diagonal until the factorisation succeeds (Levenberg-type damping of the Newton step)."""
    L, info = torch.linalg.cholesky_ex(H)
    tau = 0.0
    while int(info) != 0:
        tau = 1e-12 if tau == 0.0 else tau * 100
        if tau > 1e3:
            raise RuntimeError('Hessian could not be regularised')
        L, info = torch.linalg.cholesky_ex(H + torch.diag(tau * H.diagonal()))
    return L


def _is_pd(M):
    return int(torch.linalg.cholesky_ex(M)[1]) == 0


# ----------------------------------------------------------------------------------------
# SDP knockoffs:  max sum(s)  s.t.  2*Sigma - diag(s) >= 0,  0 <= s <= 1   (Sigma = corr.)
# ----------------------------------------------------------------------------------------
def solve_sdp_t(Corr, gap_tol=1e-6, mu=8.0, max_newton=60, margin=0.05, verbose=False):
    """Log-barrier interior-point Newton method for the MAC/SDP S-matrix.

    phi_t(s) = -t*sum(s) - logdet(2C - diag s) - sum log s - sum log(1-s).
    Each Newton step is damped so that the new iterate keeps a fraction `margin` of the
    current slack in all three constraints (M_new >= margin*M in the PSD order), which keeps
    M = 2C - diag(s) well inside the cone and its inverse numerically meaningful.
    Stops when the duality-gap bound per coordinate, 3/t, is below gap_tol. Returns s (p,)
    and diagnostics, including the duality gap `sdp_gap` of sum(s) certified by the
    dual-feasible point Y = M^{-1}/t.
    """
    p = Corr.shape[0]
    dt, dev = Corr.dtype, Corr.device
    twoC = 2.0 * Corr
    lam_min = float(torch.linalg.eigvalsh(Corr)[0])
    s = torch.full((p,), 0.5 * min(1.0, 2.0 * lam_min), dtype=dt, device=dev)
    t = 1.0
    n_newton = 0

    def phi(s_, t_):
        L_, info = torch.linalg.cholesky_ex(twoC - torch.diag(s_))
        if int(info) != 0:
            return None, None
        val = -t_ * s_.sum() - 2.0 * torch.log(L_.diagonal()).sum() - torch.log(s_).sum() - torch.log1p(-s_).sum()
        return (val if bool(torch.isfinite(val)) else None), L_

    M = twoC - torch.diag(s)
    f, L = phi(s, t)
    while True:
        for _ in range(max_newton):
            Minv = torch.cholesky_inverse(L)
            g = -t + Minv.diagonal() - 1.0 / s + 1.0 / (1.0 - s)
            H = Minv * Minv
            H.diagonal().add_(1.0 / s ** 2 + 1.0 / (1.0 - s) ** 2)
            ds = -torch.cholesky_solve(g[:, None], _chol_reg(H))[:, 0]
            gd = float((g * ds).sum())
            n_newton += 1
            if -gd / 2.0 < 1e-7:              # Newton decrement^2 / 2
                break
            # box constraints: keep `margin` of the slack
            a = 1.0
            neg, pos = ds < 0, ds > 0
            if bool(neg.any()):
                a = min(a, (1.0 - margin) * float((-s[neg] / ds[neg]).min()))
            if bool(pos.any()):
                a = min(a, (1.0 - margin) * float(((1.0 - s[pos]) / ds[pos]).min()))
            fn = None
            while a > 1e-12:
                # PSD constraint: M - a*diag(ds) >= margin*M
                if _is_pd((1.0 - margin) * M - torch.diag(a * ds)):
                    fn, Ln = phi(s + a * ds, t)
                    if fn is not None and float(fn) <= float(f) + 0.25 * a * gd:
                        break
                    fn = None
                a *= 0.5
            if fn is None:
                break
            s, f, L = s + a * ds, fn, Ln
            M = twoC - torch.diag(s)
        if verbose:
            print(f'  t={t:.3g} newton={n_newton} sum(1-s)={float((1 - s).sum()):.6f} dec={-gd / 2:.2e}')
        if 3.0 / t < gap_tol:
            break
        t *= mu
        f, L = phi(s, t)
    Y = torch.cholesky_inverse(L) / t
    dual = 2.0 * (Corr * Y).sum() + torch.clamp(1.0 - Y.diagonal(), min=0).sum()
    return s, dict(sdp_obj=float((1.0 - s).sum()), sdp_gap=float(dual - s.sum()), sdp_newton=n_newton)


def scale_until_psd_t(Sigma, s, tol, num_iter=10):
    """knockpy.utilities.scale_until_PSD for diagonal S (binary search on gamma)."""
    lo, hi = 0.0, 1.0
    I = torch.eye(Sigma.shape[0], dtype=Sigma.dtype, device=Sigma.device)
    for _ in range(num_iter):
        g = (lo + hi) / 2
        if _is_pd(2 * Sigma - g * torch.diag(s) - tol * I):
            lo = g
        else:
            hi = g
    return lo * s, lo


def compute_s_sdp_t(Sigma, tol=1e-5, **kw):
    """knockpy.smatrix.compute_smatrix(Sigma, method='sdp') for p <= max_block, on GPU."""
    scale = torch.sqrt(Sigma.diagonal())
    Corr = Sigma / torch.outer(scale, scale)
    s, info = solve_sdp_t(Corr, **kw)
    s, gamma = scale_until_psd_t(Corr, s, tol=tol)
    info['sdp_gamma'] = gamma
    return s * scale * scale, info


def gaussian_knockoffs_t(X, mu, invSigma, s, Z, sample_tol=1e-5):
    """knockpy.knockoffs.produce_MX_gaussian_knockoffs for diagonal S."""
    invS_S = invSigma * s[None, :]
    mu_k = X - (X - mu[None, :]) @ invS_S
    Vk = 2.0 * torch.diag(s) - s[:, None] * invS_S
    Lk, info = torch.linalg.cholesky_ex(Vk)
    if int(info) != 0:
        warnings.warn('Cholesky of Vk failed, shifting until PSD (FDR violations possible)')
        Lk = torch.linalg.cholesky(_shift_until_psd(Vk, sample_tol))
    return mu_k + Z @ Lk.T


# ----------------------------------------------------------------------------------------
# batched logistic sparse group lasso (groupyr LogisticSGL replica)
# ----------------------------------------------------------------------------------------
class BatchedLogisticSGL:
    """Solve K logistic SGL problems that share X but differ in alpha and sample weights.

    Problem k:  min_{beta,c} sum_i w_ik [(1-b_i) z_ik - log sigmoid(z_ik)] + P_k(beta),
    z_ik = x_i.beta + c, where w_k = 1/n_train on the training rows of problem k (0 elsewhere).
    """

    def __init__(self, X, groups, l1_ratio, exact_prox=True, feat_w=None):
        """feat_w: optional per-feature multiplier of the l1 penalty (weighted l1: sum_j feat_w_j |beta_j|)."""
        self.X = X
        self.feat_w = None if feat_w is None else torch.as_tensor(np.asarray(feat_w), dtype=X.dtype, device=X.device)[:, None]
        n, d = X.shape
        dev = X.device
        g_of = torch.full((d,), -1, dtype=torch.long)
        for j, g in enumerate(groups):
            g_of[torch.as_tensor(np.asarray(g), dtype=torch.long)] = j
        if bool((g_of < 0).any()):
            raise ValueError('every feature must belong to a group')
        self.g_of = g_of.to(dev)
        self.G = len(groups)
        sizes = torch.tensor([len(g) for g in groups], dtype=X.dtype)
        self.gscale = torch.sqrt(sizes).to(dev)[:, None]
        self.l1_ratio = float(l1_ratio)
        self.exact_prox = exact_prox

    def prox(self, x, step, alpha):
        """Prox of the sparse group lasso penalty, column-wise. x: (d, K); step, alpha: (K,).
        Soft-threshold by l1_ratio*alpha*step, then shrink each group towards 0 by
        (1-l1_ratio)*alpha*step*sqrt(|g|) in l2 norm (exact_prox=True), or groupyr's variant."""
        thr1 = self.l1_ratio * alpha * step
        if self.feat_w is not None:
            thr1 = thr1 * self.feat_w
        u = torch.sign(x) * torch.clamp(x.abs() - thr1, min=0)
        sq = torch.zeros((self.G, x.shape[1]), dtype=x.dtype, device=x.device).index_add_(0, self.g_of, u * u)
        norms = torch.sqrt(sq) / self.gscale
        thr2 = (1.0 - self.l1_ratio) * alpha * step
        mask = norms > thr2
        shrink = 1.0 - thr2 / torch.where(mask, norms, torch.ones_like(norms))
        if self.exact_prox == 'hard':
            # group hard-thresholding: a group is kept unchanged (after the l1 step) or set to zero
            fac = mask.to(norms.dtype)
        elif self.exact_prox:
            fac = torch.where(mask, shrink, torch.zeros_like(norms))
        else:
            # groupyr applies the group shrinkage only when *every* group survives;
            # otherwise surviving groups keep the l1-thresholded values unchanged.
            allm = mask.all(0, keepdim=True)
            fac = torch.where(mask, torch.where(allm, shrink, torch.ones_like(norms)), torch.zeros_like(norms))
        return u * fac[self.g_of]

    def _z(self, B, c, cols=None):
        return self.X @ B + c[None, :]

    @staticmethod
    def _f(z, b, W):
        return (W * ((1.0 - b)[:, None] * z - torch.nn.functional.logsigmoid(z))).sum(0)

    def _grad(self, z, b, W):
        R = W * (torch.sigmoid(z) - b[:, None])
        return self.X.T @ R, R.sum(0)

    def solve(self, b, W, alphas, max_iter=1000, tol=1e-7, bt_factor=0.6, max_bt=1000):
        X = self.X
        d = X.shape[1]
        K = W.shape[1]
        dt, dev = X.dtype, X.device
        alphas = _to_t(alphas, dev, dt)
        B = torch.zeros((d, K), dtype=dt, device=dev)
        c = torch.full((K,), 0.5, dtype=dt, device=dev)
        z = self._z(B, c)
        f = self._f(z, b, W)
        gB, gc = self._grad(z, b, W)

        # copt.utils.init_lipschitz
        L0 = torch.full((K,), 1e-3, dtype=dt, device=dev)
        for _ in range(100):
            ft = self._f(self._z(B - gB / L0, c - gc / L0), b, W)
            bad = ft > f
            if not bool(bad.any()):
                break
            L0 = torch.where(bad, L0 * 10, L0)
        step = 1.8 / L0

        active = torch.ones(K, dtype=torch.bool, device=dev)
        n_iter = torch.zeros(K, dtype=torch.long, device=dev)
        for it in range(max_iter + 1):
            act = torch.nonzero(active).flatten()
            if act.numel() == 0:
                break
            Ba, ca, fa, gBa, gca = B[:, act], c[act], f[act], gB[:, act], gc[act]
            Wa, al, st = W[:, act], alphas[act], step[act]
            Bn = self.prox(Ba - st * gBa, st, al)
            cn = ca - st * gca
            st = st * 1.1
            zn = self._z(Bn, cn)
            fn = self._f(zn, b, Wa)
            pend = torch.ones(act.numel(), dtype=torch.bool, device=dev)
            for _ in range(max_bt):
                dB, dc = Bn - Ba, cn - ca
                rhs = fa + (gBa * dB).sum(0) + gca * dc + ((dB * dB).sum(0) + dc * dc) / (2.0 * st)
                pend = pend & ~(fn <= rhs)
                if not bool(pend.any()):
                    break
                pi = torch.nonzero(pend).flatten()
                st = torch.where(pend, st * bt_factor, st)
                Bp = self.prox(Ba[:, pi] - st[pi] * gBa[:, pi], st[pi], al[pi])
                cp = ca[pi] - st[pi] * gca[pi]
                zp = self._z(Bp, cp)
                Bn[:, pi], cn[pi], zn[:, pi] = Bp, cp, zp
                fn[pi] = self._f(zp, b, Wa[:, pi])
            gBn, gcn = self._grad(zn, b, Wa)
            dB, dc = Ba - Bn, ca - cn
            cert = torch.sqrt((dB * dB).sum(0) + dc * dc) / st
            B[:, act], c[act], f[act], gB[:, act], gc[act], step[act] = Bn, cn, fn, gBn, gcn, st
            n_iter[act] = it
            active[act] = ~(cert < tol)
        return B, c, n_iter


# ----------------------------------------------------------------------------------------
# scKnockPath on GPU
# ----------------------------------------------------------------------------------------
class scKnockPathGPU(scKnockPath):
    def __init__(self, device='cuda', dtype=torch.float64, exact_prox=True) -> None:
        super().__init__()
        self.device = torch.device(device)
        self.dtype = dtype
        self.exact_prox = exact_prox
        self.timings = {}

    # ---------------- 1. pathway unfolding ----------------
    def duplicate_genes(self, X):
        if not hasattr(self, '_unfold_idx') or getattr(self, '_unfold_key', None) is not id(self.pathway_dict):
            idx, groups_list, group_gene_names = [], [], []
            num_col = 0
            for group in self.pathway_dict.values():
                cols = np.where(np.isin(self.gene_names, group))[0]
                idx.append(cols)
                groups_list.append(np.arange(num_col, num_col + len(cols)))
                group_gene_names.append(self.gene_names[cols])
                num_col += len(cols)
            counts = {}
            for names in group_gene_names:
                for g in names:
                    counts[g] = counts.get(g, 0) + 1
            self.group_gene_names = group_gene_names
            self.unique_cols = [cols[np.array([counts[g] == 1 for g in names], dtype=bool)]
                                for cols, names in zip(groups_list, group_gene_names)]
            self._groups_list = groups_list
            self._unfold_idx = torch.as_tensor(np.concatenate(idx), dtype=torch.long, device=self.device)
            self._unfold_key = id(self.pathway_dict)
        Xt = X.to(self.device) if isinstance(X, torch.Tensor) else _to_t(X, self.device, self.dtype)
        return Xt[:, self._unfold_idx], self._groups_list

    def prepare_data(self, *args, **kwargs):
        """Same as the CPU version, but X and Xh are returned as tensors on the device."""
        t0 = time.time()
        X, Xh, y, groups_list = super().prepare_data(*args, **kwargs)
        X = _to_t(X, self.device, self.dtype)
        self.timings['prepare_sec'] = time.time() - t0
        return X, Xh, y, groups_list

    # ---------------- 2. knockoffs ----------------
    def estimate_sigma_t(self, X, estimator='ledoitwolf', ridge=1e-2, glasso_alpha=0.1, Sigma_ref=None):
        info = {'cov_estimator': estimator}
        if estimator == 'ledoitwolf':
            Sigma, invSigma, sh = ledoit_wolf_t(X)
            info['lw_shrinkage'] = sh
        elif estimator == 'ridge':
            Xc = X - X.mean(0)
            Sigma = Xc.T @ Xc / (X.shape[0] - 1)
            Sigma.diagonal().add_(ridge)
            invSigma = torch.cholesky_inverse(torch.linalg.cholesky(Sigma))
            info['ridge'] = ridge
        elif estimator in ('graphicallasso', 'reference'):
            Xn = X.cpu().numpy()
            keep = Xn.std(0) > 0 if estimator == 'graphicallasso' else np.ones(Xn.shape[1], bool)
            S_k, inv_k, info = scKnockPath.estimate_sigma(Xn[:, keep], estimator=estimator,
                                                          glasso_alpha=glasso_alpha, Sigma_ref=Sigma_ref)
            if keep.all():
                S_np, inv_np = S_k, inv_k
            else:
                # constant genes (zero columns after standardisation) make the graphical lasso ill-conditioned:
                # fit it on the other genes and give the constant ones unit variance, independent of the rest
                S_np = np.eye(Xn.shape[1]); inv_np = np.eye(Xn.shape[1])
                ix = np.ix_(keep, keep); S_np[ix] = S_k; inv_np[ix] = inv_k
                info['n_constant_genes'] = int((~keep).sum())
            return _to_t(S_np, self.device, X.dtype), _to_t(inv_np, self.device, X.dtype), info
        else:
            raise ValueError(estimator)
        info['min_eig'] = float(torch.linalg.eigvalsh(Sigma)[0])
        return Sigma, invSigma, info

    def knockoff_sampler(self, X, Xh=None, method='sdp', cov_estimator='ledoitwolf', num_processes=None,
                         Sigma_ref=None, S=None, sdp_kwargs=None, **cov_kwargs):
        """Returns Xc = [Xh_std, Xkh] (tensor on device), X_std, Xk.
        S: optional precomputed diagonal of the S-matrix (e.g. to reuse across fits)."""
        print('start generating knockoffs (GPU)')
        dev = self.device
        t0 = time.time()
        # the knockoff construction always runs in float64 (the SDP does not converge in float32)
        X = standardize_t(_to_t(X, dev, torch.float64))
        if method == 'sdp':
            Sigma, invSigma, info = self.estimate_sigma_t(X, estimator=cov_estimator, Sigma_ref=Sigma_ref, **cov_kwargs)
            _sync(dev); t1 = time.time()
            if S is None:
                s, sinfo = compute_s_sdp_t(Sigma, **(sdp_kwargs or {}))
                info.update(sinfo)
            else:
                s = _to_t(S, dev, X.dtype)
            _sync(dev); t2 = time.time()
            Z = _to_t(np.random.randn(X.shape[0], X.shape[1], 1)[:, :, 0], dev, X.dtype)
            Xk = gaussian_knockoffs_t(X, X.mean(0), invSigma, s, Z)
            _sync(dev); t3 = time.time()
            s_np = s.cpu().numpy()
            info['s_mean'] = float(np.mean(s_np)); info['s_min'] = float(np.min(s_np))
            info['corr_orig_knock_mean'] = float(np.mean(1 - s_np))
            self.S_diag = s_np
            self.timings.update(cov_sec=t1 - t0, sdp_sec=t2 - t1, sample_sec=t3 - t2)
        elif method == 'permutation':
            perm = torch.as_tensor(np.random.permutation(X.shape[0]), device=dev)
            Xk = X[perm]
            info = {'method': 'permutation'}
        else:
            raise ValueError(method)
        self.knockoff_info = info
        Xh_s, _ = self.duplicate_genes(X)          # = StandardScaler(Xh): columns are copies of X's
        Xkh, _ = self.duplicate_genes(Xk)
        Xc = torch.cat([Xh_s, Xkh], dim=1).to(self.dtype)
        _sync(dev)
        self.timings['knockoff_sec'] = time.time() - t0
        print(f'shape of knock off={tuple(Xk.shape)}; {info}')
        sys.stdout.flush()
        return Xc, X, Xk

    # ---------------- 3. batched SGL fit ----------------
    def fit(self, Xc, y, groups_list, cv=None, n_alphas=10, alphas=None, l1_ratio=0.2, merge_groups=False,
            n_jobs=None, verbose=True, max_iter=1000, tol=1e-7, store_path=False, cv_groups=None, shared_l1_power=0.0,
            cv_gain=1.0):
        """cv_gain: lambda rule. 1.0 = lambda with the smallest CV loss. x < 1 = the LARGEST lambda whose CV-loss
        improvement over the null model (largest lambda of the grid) reaches x times the best improvement, e.g. 0.9.
        The CV loss is flat near its minimum, and the minimum lies at a lambda so small that redundant pathways enter.
        cv_groups: optional group label per sample (e.g. donor). The CV folds then keep each group together
        (StratifiedGroupKFold), so lambda is chosen by prediction on unseen groups.
        shared_l1_power = gamma > 0: the l1 penalty of every copy of a gene that lies in m pathways is multiplied by
        m**gamma, so that a pathway enters the model less easily through genes it shares with other pathways.
        store_path=True additionally fits the full data at every lambda of the path (one batch) and keeps
        Z, Z~, Z^U, Z~^U per lambda in self.path, so that other lambda rules can be evaluated afterwards."""
        Xc = _to_t(Xc, self.device, self.dtype)
        y = np.asarray(y)
        self.unfold_length = int(Xc.shape[1] // 2)
        self.l1_ratio = l1_ratio
        self.merge_groups = merge_groups
        if merge_groups:
            total_groups_list = [np.concatenate([g, g + self.unfold_length]) for g in groups_list]
        else:
            total_groups_list = list(groups_list) + [g + self.unfold_length for g in groups_list]

        alphas = np.logspace(-1, -3, num=n_alphas) if alphas is None else alphas
        alphas = np.array([float(f'{a:.2g}') for a in alphas])
        print(f'Path Range: {alphas[0]} -> {alphas[-1]}, n_alphas={len(alphas)}')
        self.alphas = alphas

        classes = np.unique(y)
        if len(classes) != 2:
            raise ValueError('binary labels required')
        b = _to_t((y != classes[0]).astype(float), self.device, self.dtype)
        n = len(y)
        feat_w = None
        if shared_l1_power:
            gi = self._unfold_idx.cpu().numpy()
            mult = np.bincount(gi)[gi].astype(float) ** shared_l1_power
            feat_w = np.concatenate([mult, mult])
        self.shared_l1_power = shared_l1_power
        solver = BatchedLogisticSGL(Xc, total_groups_list, l1_ratio, exact_prox=self.exact_prox, feat_w=feat_w)
        start = time.time()

        def weights(rows_list):
            W = torch.zeros((n, len(rows_list)), dtype=self.dtype, device=self.device)
            for k, rows in enumerate(rows_list):
                W[torch.as_tensor(rows, device=self.device), k] = 1.0 / len(rows)
            return W

        if cv is None:
            print(f'no CV, processing alpha={alphas[0]}')
            best_alpha = float(alphas[0]); self.cv_results = None
        else:
            if cv_groups is None:
                splits = list(KFold(n_splits=cv, shuffle=True, random_state=42).split(np.zeros(n)))
            else:
                from sklearn.model_selection import StratifiedGroupKFold
                splits = list(StratifiedGroupKFold(n_splits=cv, shuffle=True, random_state=42).split(np.zeros(n), y, np.asarray(cv_groups)))
            for tr, _ in splits:
                if len(np.unique(y[tr])) != 2:
                    raise ValueError('a training fold contains one class only')
            rows, al = [], []
            for a in alphas:
                for tr, _ in splits:
                    rows.append(tr); al.append(a)
            Bk, ck, it = solver.solve(b, weights(rows), np.array(al), max_iter=max_iter, tol=tol)
            self.cv_n_iter = it.cpu().numpy()
            k = 0
            results = []
            for a in alphas:
                fs = []
                for _, te in splits:
                    te_t = torch.as_tensor(te, device=self.device)
                    p1 = torch.sigmoid(Xc[te_t] @ Bk[:, k] + ck[k]).cpu().numpy()
                    fs.append(-log_loss(y[te], np.column_stack([1 - p1, p1]), labels=[0, 1]))
                    k += 1
                results.append({'alpha': a, 'mean_score': float(np.mean(fs)), 'fold_scores': fs})
            best = max(results, key=lambda r: r['mean_score'])
            self.best_alpha_cvmin = float(best['alpha']); self.cv_gain = cv_gain
            if cv_gain < 1.0:
                sc_ = np.array([r['mean_score'] for r in results]); gain = sc_ - sc_[0]
                if gain.max() > 0:
                    best = results[int(np.where(gain >= cv_gain * gain.max())[0][0])]
            best_alpha = float(best['alpha'])
            self.cv_results = [{'alpha': float(r['alpha']), 'mean_score': float(r['mean_score']),
                                'fold_scores': [float(x) for x in r['fold_scores']]} for r in results]
            print(f"Best l1_ratio: {l1_ratio}, Best Alpha: {best_alpha}, CV Mean Score: {best['mean_score']:.4f}")

        B, c, it = solver.solve(b, weights([np.arange(n)]), np.array([best_alpha]), max_iter=max_iter, tol=tol)
        self.n_iter_ = int(it[0])
        model = _FittedSGL(B[:, 0].cpu().numpy(), float(c[0]), classes)
        self.coef_ = model.coef_
        self.best_alpha = best_alpha
        self.path = None
        if store_path:
            Bp, cp_, itp = solver.solve(b, weights([np.arange(n)] * len(alphas)), np.asarray(alphas), max_iter=max_iter, tol=tol)
            self.path = dict(alphas=[float(a) for a in alphas], n_iter=itp.cpu().numpy().tolist(), Z=[], Zk=[], ZU=[], ZkU=[],
                             ZC=[], ZkC=[], df_coef=[], df_gene=[], n_active_groups=[])
            # mean negative log-likelihood of the full-data fit at every lambda (for BIC-type criteria)
            self.path['nll'] = torch.nn.functional.binary_cross_entropy_with_logits(
                Xc @ Bp + cp_[None, :], b[:, None].expand(-1, len(alphas)), reduction='none').mean(0).cpu().numpy().tolist()
            Bp = Bp.cpu().numpy(); gi2 = self._unfold_idx.cpu().numpy(); gi2 = np.concatenate([gi2, gi2 + gi2.max() + 1])
            for k in range(len(alphas)):
                self.update_Wg(_FittedSGL(Bp[:, k], float(cp_[k]), classes))
                for key in ('Z', 'Zk', 'ZU', 'ZkU'):
                    self.path[key].append(getattr(self, key).tolist())
                zc, zkc = self.claim_stats(Bp[:, k])
                self.path['ZC'].append(zc.tolist()); self.path['ZkC'].append(zkc.tolist())
                nz = Bp[:, k] != 0
                self.path['df_coef'].append(int(nz.sum()))                    # nonzero coefficients (unfolded copies)
                self.path['df_gene'].append(int(len(np.unique(gi2[nz]))))     # distinct columns (gene or knockoff gene) in the model
                self.path['n_active_groups'].append(int(sum(nz[g].any() for g in total_groups_list)))
        self.update_Wg(model)
        if verbose:
            print('scKnockPath selects', self.select_pathway(fdr=0.2))
        _sync(self.device)
        self.fit_time_sec = time.time() - start
        self.timings['fit_sec'] = self.fit_time_sec
        print(f'Total time taken: {self.fit_time_sec / 60:.2f} minutes')

    def claim_stats(self, coef=None):
        """Z^C, Z~^C: as Z, Z~ (mean |coef| over the pathway's genes), but every gene that lies in several pathways is
        credited to ONE of them only: the pathway h with the largest T_h = Z_h + Z~_h (symmetric in original/knockoff).
        A pathway that is in the model only through genes it shares with a stronger pathway gets Z^C ~ 0."""
        co = np.abs(self.coef_ if coef is None else coef); L = self.unfold_length; gi = self._unfold_idx.cpu().numpy()
        T = np.asarray(self.Z) + np.asarray(self.Zk)
        g_of = np.empty(L, dtype=int)
        for i, g in enumerate(self._groups_list):
            g_of[g] = i
        best = np.full(gi.max() + 1, -np.inf); owner = np.full(gi.max() + 1, -1)
        for c in np.argsort(-T[g_of], kind='stable'):
            if owner[gi[c]] < 0:
                owner[gi[c]] = g_of[c]
        own = owner[gi] == g_of
        ZC = np.array([co[g][own[g]].sum() / len(g) for g in self._groups_list])
        ZkC = np.array([co[g + L][own[g]].sum() / len(g) for g in self._groups_list])
        return ZC, ZkC

    def save_model(self, file_name):
        # tensors/devices are not needed after fitting; keep the pickle loadable on CPU-only machines
        for k in ('_unfold_idx',):
            self.__dict__.pop(k, None)
        self.device = str(self.device)
        super().save_model(file_name)


class _FittedSGL:
    def __init__(self, coef, intercept, classes):
        self.coef_ = coef
        self.intercept_ = intercept
        self.classes_ = classes
