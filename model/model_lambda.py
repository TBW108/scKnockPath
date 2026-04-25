import knockpy
import numpy as np
import pandas as pd
import json
import warnings
import pickle
import sys
from groupyr import LogisticSGL
from collections import OrderedDict
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import gseapy as gp
from sklearn.model_selection import KFold
from joblib import Parallel, delayed
from sklearn.base import clone
from sklearn.metrics import log_loss
import time
from sklearn.model_selection import GridSearchCV


class scKnockPath:
    def __init__(self,) -> None:
        pass
    
    def get_pathways_dict(self,file_name):
        '''
        read geneset file and return a dict
        '''
        if file_name.endswith('.json'):
            geneset_dict=dict()
            # get the genes for each gene set
            with open(file_name, 'r') as json_file:
                data=json.load(json_file)
                for pathway in data.keys():
                    geneset_dict[pathway]=np.array(data[pathway]['geneSymbols'])
                json_file.close()
            print(f'Original file has {len(geneset_dict)} pathways.')
            
        elif file_name.endswith('.gmt'):
            geneset_dict=gp.read_gmt(file_name)
            # Convert each value (gene list) to numpy array
            for key in geneset_dict:
                geneset_dict[key] = np.array(geneset_dict[key])
            print(f'Original file has {len(geneset_dict)} pathways.')

        return geneset_dict

    def filter_genes(self,geneset_dict,X_genenames):
        '''
        filter genesets based on gene_thresh and genes in X
        '''
        filter_genesets=dict()
        for pathway,genes in geneset_dict.items():
            # if the size of the geneset is larger than thresh and genes of it are in adata's genes
            index=np.isin(genes.astype(str),X_genenames.astype(str))
            if index.sum()>=self.gene_thresh:
                filter_genesets[pathway]=genes[index]
        print(f'After filtering, there are {len(filter_genesets)} pathways in total.')
        
        return filter_genesets
     
    def duplicate_genes(self,X):
        '''
        construct the unfolded matrix Xh and the groups_list of Xh
        '''
        # rearange X to get Xh
        groups_list=[]# the index of cols for each group in X_hat
        n_total_col=sum([genes.shape[0] for genes in self.pathway_dict.values()])
        Xh=np.zeros((X.shape[0],n_total_col)) # prepare the duplicated matrix X_hat
        num_col=0 # count for the loop below
        
        for group in self.pathway_dict.values():
            gene_indices=np.isin(self.gene_names,group)
            Xh[:,num_col:num_col+group.size]=X[:,gene_indices]
            groups_list.append(np.arange(num_col,num_col+group.size))
            num_col+=group.size

        return Xh,groups_list
    
    def knockoff_sampler(self,X,Xh,method='sdp'):
        '''
        get knock off matrix Xk for Xkh
        '''
        print('start generating knock offs')
        # duplicate overlapping genes to generate Xk_h, data structure of the X should be changed
        if method=='sdp':
            X = StandardScaler().fit_transform(X)
            Xh = StandardScaler().fit_transform(Xh)
            sampler=knockpy.knockoffs.GaussianSampler(X,method=method,max_block=3000, num_processes=100,)
            Xk=sampler.sample_knockoffs()
            print(f'shape of knock off={Xk.shape}')
            Xkh,_=self.duplicate_genes(Xk,)
            sys.stdout.flush()
            
        elif method=='permutation':
            print("Using permutation method to generate knockoff data...")
            X = StandardScaler().fit_transform(X)
            Xh = StandardScaler().fit_transform(Xh)
            
            n_samples = X.shape[0]
            random_labels = np.random.permutation(n_samples)
            #  生成 knockoff 数据
            Xk = X[random_labels, :].copy()
            print(f'shape of knock off={Xk.shape}')
            Xkh,_=self.duplicate_genes(Xk,)
            sys.stdout.flush()
            
        Xc=np.hstack([Xh,Xkh])
        
        return Xc,X,Xk
    
    def prepare_data(self,adata,obs_y,genesets,layer=None,gene_thresh=15,gene_names=None,class1=None,class2=None):
        '''
        return the X,Xh,y,groups_list to input into the algorithm
        '''
        self.gene_thresh=gene_thresh
        
        # get gene names
        if gene_names is None:
            self.gene_names=np.array(adata.var_names)
        else:
            self.gene_names=np.array(adata.var[gene_names])
        
        # get gene set dict
        if isinstance(genesets,str):
            geneset_dict=self.get_pathways_dict(genesets)
            self.pathway_dict=self.filter_genes(geneset_dict,self.gene_names)
        else:
            self.pathway_dict=self.filter_genes(genesets,self.gene_names)
        self.pathway_names=np.array(list(self.pathway_dict.keys()))
        
        # get X, y used for calculation
        adata=adata[(adata.obs[obs_y]==class2)|(adata.obs[obs_y]==class1),:]
        if layer is None:
            if hasattr(adata.X, 'toarray'):
                X=np.array(adata.X.toarray())
            else:
                X=np.array(adata.X)
        else:
            if hasattr(adata.layers[layer], 'toarray'):
                X=np.array(adata.layers[layer].toarray())
            else:
                X=np.array(adata.layers[layer])
            
        Xh,groups_list=self.duplicate_genes(X)
        y=np.array(adata.obs[obs_y])
        self.groups_list=groups_list
        
        print(f"(After sampling)n_cells={X.shape[0]},n_features={X.shape[1]},n_unfold_features={Xh.shape[1]}, n_pathways={len(self.pathway_dict)}")
        sys.stdout.flush()
        
        return X,Xh,y,groups_list
    
    def get_robust_lambda_path(self,Xc, y, n_alphas=10, eps=0.001):
        """
        计算稳健的 Lambda 路径，专门解决 FDR 高的问题。
        
        参数:
        X : (n_samples, n_features) 设计矩阵 (包含原始 + Knockoff)
        y : (n_samples, ) 0/1 标签
        eps : lambda_min / lambda_max 的比例。
            建议设为 0.05 而不是默认的 0.001。
            eps 越大，路径截断得越早，越能排除尾部噪音。
        """
        n_samples = Xc.shape[0]
        
        # 2. 计算 Null Model (只含截距) 的预测概率
        # 对于 Logistic，这就是正样本的比例
        y = y.astype(float)
        p_null = np.mean(y)
        
        # 3. 计算残差 (Residuals)
        residuals = y - p_null
        
        # 4. 计算 Lambda Max (理论上限)
        # 公式：特征与残差内积的最大绝对值 / 样本数
        # 这代表了“最强变量”的初始信号强度
        grad = np.abs(Xc.T @ residuals)
        lambda_max = np.max(grad) / n_samples
        
        # 5. 计算 Lambda Min (截断点)
        # 如果您的 FDR 很高，请尝试调大 eps (如 0.05 或 0.1)
        # 这相当于告诉模型：“太微弱的信号我不要了，只要头部的”
        lambda_min = lambda_max * eps
        
        # 6. 生成对数尺度的路径 (100个点)
        alphas = np.geomspace(lambda_max, lambda_min, n_alphas)
        
        return alphas, lambda_max

    
    def fit(self,Xc,y,groups_list,cv=None,n_alphas=10,alphas=None,l1_ratio=0.2):
        '''
            repeat the selection
        '''
        # total group list
        self.unfold_length=int(Xc.shape[1]//2)
        total_groups_list=groups_list.copy()
        for group in groups_list:
            total_groups_list.append(group+self.unfold_length)
        
        warnings.filterwarnings("ignore", category=DeprecationWarning)

        best_score = -np.inf
        best_alpha = None
        self.sgl_models = []
        self.cv_scores = []
        
        # if alphas is None:
        #     alphas, max_val = self.get_robust_lambda_path(Xc, y, n_alphas=n_alphas, eps=eps)
            
        # 为保证数值稳定，取每个 alpha 的前两位非 0 的小数
        alphas = np.logspace(-1, -3, num=n_alphas) if alphas is None else alphas
        alphas = np.array([float(f"{alpha:.2g}") for alpha in alphas])
        print(f"Path Range: {alphas[0]} -> {alphas[-1]}, n_alphas={len(alphas)}")
        self.alphas = alphas
        
        # grid—search 在评估时使用的是模型的 score 方法，而 LogisticSGL 的 score 方法默认是 accuracy，这可能不适合不平衡数据集。建议在 GridSearchCV 中使用 log_loss 作为评估指标，以更好地反映模型在不平衡数据上的性能。
        start_time = time.time()
        if cv is None:
            
            print(f'no CV, processing alpha={alphas[0]}')
            # create the model
            sgl_model = LogisticSGL(l1_ratio=l1_ratio, groups=total_groups_list, alpha=alphas[0])
            sgl_model.fit(Xc, y)
            
            self.sgl_models.append(sgl_model)
            self.update_Wg(sgl_model)
            print('scKnockPath selects',self.select_pathway(fdr=0.2,offset=0))

            sys.stdout.flush()
            end_time = time.time()
            print(f"Total time taken: {(end_time - start_time)/60:.2f} minutes")
            
        else:
            # param_grid={'alpha':alphas,
            #             'l1_ratio':[0,0.1,0.2]
            # }
            # gs_cv=GridSearchCV(LogisticSGL(groups=total_groups_list,),param_grid=param_grid,cv=cv,n_jobs=-1,scoring='neg_log_loss')
            # gs_cv.fit(Xc,y)
            # Define a custom scorer for negative log loss (since GridSearchCV maximizes the score)
            # l1_ratios = [0.1, 0.2, 0.3]
            # neg_log_loss_scorer = make_scorer(log_loss, greater_is_better=False, needs_proba=True)

            # param_grid = {
            #     'alpha': alphas,
            #     'l1_ratio': [0.1, 0.2, 0.3]
            # }
            # gs_cv = GridSearchCV(
            #     LogisticSGL(groups=total_groups_list),
            #     param_grid=param_grid,
            #     cv=cv,
            #     n_jobs=-1,
            #     scoring=neg_log_loss_scorer
            # )
            # gs_cv.fit(Xc, y)
            
            # best_model=gs_cv.best_estimator_

            # find the best alpha and l1_ratio
            results = []
            # Prepare KFold splits once for reproducibility
            kf = KFold(n_splits=cv, shuffle=True, random_state=42)
            splits = list(kf.split(Xc))
            # l1_ratio=0.05
            def fit_cv(alpha):
                # For each fold, fit and score
                fold_scores = []
                for train_index, test_index in splits:
                    model = LogisticSGL(l1_ratio=l1_ratio, groups=total_groups_list, alpha=alpha, warm_start=True)
                    X_train, X_test = Xc[train_index], Xc[test_index]
                    y_train, y_test = y[train_index], y[test_index]
                    model.fit(X_train, y_train)
                    if not hasattr(model, 'classes_'):
                        model.classes_ = np.unique(y_train)
                    y_pred_proba = model.predict_proba(X_test)
                    fold_scores.append(-log_loss(y_test, y_pred_proba))
                mean_score = np.mean(fold_scores)
                return {
                    'alpha': alpha,
                    'mean_score': mean_score,
                    'fold_scores': fold_scores
                }

            # Run all (l1_ratio, alpha) combinations in parallel
            # param_grid = [(l1, a) for l1 in l1_ratios for a in alphas]
            parallel_results = Parallel(n_jobs=-1)(
                delayed(fit_cv)(alpha) for alpha in alphas
            )

            # Find the best combination
            best_result = max(parallel_results, key=lambda x: x['mean_score'])
            best_l1_ratio = l1_ratio  # Fixed l1_ratio since it's not part of the grid search
            best_alpha = best_result['alpha']
            best_score = best_result['mean_score']

            print(f'Best l1_ratio: {best_l1_ratio}, Best Alpha: {best_alpha}, CV Mean Score: {best_score:.4f}')

            # Fit the best model on the whole data
            best_model = LogisticSGL(l1_ratio=best_l1_ratio, groups=total_groups_list, alpha=best_alpha, warm_start=True)
            best_model.fit(Xc, y)

            self.update_Wg(best_model)
            print('scKnockPath selects',self.select_pathway(fdr=0.2,offset=0))
            end_time = time.time()
            print(f"Total time taken: {(end_time - start_time)/60:.2f} minutes")

    def fit_W(self,Xc,y,groups_list):
        Xh=Xc[:,:Xc.shape[1]//2]
        Xhk=Xc[:,Xc.shape[1]//2:]
        Wh=np.abs(Xh.T.dot(y - y.mean())) - np.abs(Xhk.T.dot(y - y.mean()))# 每个基因的 W 统计量
        
        self.Wg=np.zeros(len(groups_list))
        # 每一个 pathway都有一个 Wg，计算方式是 pathway 内基因的 W 的平均值
        for i,group in enumerate(groups_list):
            self.Wg[i]=Wh[group].mean()

    def update_Wg(self,sgl_model):
        # get betas, Wg
        n_groups = len(self.groups_list)
        Z=np.zeros(n_groups)
        Zk=np.zeros(n_groups)
        Wg=np.zeros(n_groups)
        for i,group in enumerate(self.groups_list):
            Z[i]=np.abs(sgl_model.coef_[group]).sum()/group.shape[0]
            Zk[i]=np.abs(sgl_model.coef_[group+self.unfold_length]).sum()/group.shape[0]
            Wg[i] = np.sign(Z[i] - Zk[i]) * (Z[i]-Zk[i])
            sys.stdout.flush()
        self.Z=Z
        self.Zk=Zk
        self.Wg=Wg
        self.beta=sgl_model.coef_
    
    def save_model(self,file_name):
        with open(file_name, 'wb') as f:
            pickle.dump(self, f)
        print('save success')
        
    def rank_pathways(self):
        "rank pathways based on mean of Wgs"
        rank_scores=self.Wg
        df=pd.DataFrame({
            'pathway_names':self.pathway_names,
            'scores':rank_scores
        })
        df.sort_values(by='scores',inplace=True,ascending=True)
        
        self.pathway_rank=df
        
        return df

    def select_pathway(self,fdr,offset=0):
        "select pathways based on the T"
        self.Tg=knockpy.knockoff_stats.data_dependent_threshhold(self.Wg,fdr=fdr,offset=offset)
        self.selected_id = self.Wg>=self.Tg
        
        # print(f'Select {np.sum(self.selected_id)} pathways')
        self.pathway_names=np.array(self.pathway_names)
        self.selected_pathways=self.pathway_names[self.selected_id]

        return self.selected_pathways
    
    def transform_geneset_dict(self,):
        '''
        make geneset_dict into a dict where keys are genes and values represents the pathways the gene belongs to.
        '''
        overlap_gene_dict=OrderedDict()
        for gene in self.gene_names:
            overlap_gene_dict[gene]=[]
            for pathway,pathway_genes in self.pathway_dict.items():
                if gene in pathway_genes:
                    overlap_gene_dict[gene].append(pathway)
        
        self.overlap_gene_dict=overlap_gene_dict
        return self.overlap_gene_dict
    
    
    def find_gene_effect(self,gene=None,pathway=None):
        '''
        Given a gene, the function offers the effect of the gene in different pathways.
        Given a pathway, the function offers the effect of genes in the pathway.
        Parameters:
        - gene: str, the name of the gene.
        - pathway: str, the name of the pathway.

        Returns:
        - df: pd.DataFrame, a data frame containing the effect of genes in different pathways.
        '''
        # get effect for every gene
        n_all_features=self.beta.shape[0]
        
        self.beta=np.array(self.beta)
        W=np.abs(self.beta[:n_all_features//2])-np.abs(self.beta[n_all_features//2:])
        self.genes_effect=W
        
        if (pathway is not None) & (gene is None):
            pathway_ind=np.where(self.pathway_names==pathway)[0][0]
            genes_ind=self.groups_list[pathway_ind]
            target_effect=self.genes_effect[genes_ind]
            df=pd.DataFrame({
                'genes_names':self.pathway_dict[pathway],
                'genes_effect':target_effect})
        
        if gene is not None:
            
            pathway_effect=[]
            overlap_pathway=self.overlap_gene_dict[gene]
            for pathway in overlap_pathway:
                gene_ind=np.where(self.pathway_dict[pathway]==gene)[0][0]
                pathway_ind=np.where(self.pathway_names==pathway)[0][0]
                pathway_effect.append(self.genes_effect[self.groups_list[pathway_ind][gene_ind]])
        
            df=pd.DataFrame({
                'pathway_names':overlap_pathway,
                'pathway_effect':pathway_effect})
        
        return df
    
    def load_model(self,file_name):
        with open(file_name, 'rb') as f:
            loaded_model = pickle.load(f)
        self.__dict__.update(loaded_model.__dict__)
        print('load success')
    
    def plot_JImat(self,sig_pathways,geneset_dict=None,show_text=False,show_axes=False,axes_size=10,plot=False,cmap='Reds',JI_th=-0.1,x_ticks=None):
        
        if geneset_dict is None:
            geneset_dict=self.pathway_dict
        
        def jaccard_index(arr1, arr2):
            intersection = np.intersect1d(arr1, arr2)
            union = np.union1d(arr1, arr2)
            return intersection.size / union.size
        
        sig_pathways_num=len(sig_pathways)
        # 如果一个pathway和其他的pathway的JI只要有一个大于阈值，就保留下来画图
        JI_mat=np.zeros((sig_pathways_num,sig_pathways_num))
        for i,pathway1 in enumerate(sig_pathways):
            for j,pathway2 in enumerate(sig_pathways):
                JI_mat[i,j]=jaccard_index(geneset_dict[pathway1],geneset_dict[pathway2],)

        np.fill_diagonal(JI_mat, 0)
        target_pathway_idx = np.any(JI_mat > JI_th, axis=0)
        sig_pathways = np.array(sig_pathways)[target_pathway_idx]
        sig_pathways_num=len(sig_pathways)
        JI_mat = JI_mat[target_pathway_idx][:, target_pathway_idx]
        
        if plot:
            plt.figure(dpi=300)
            plt.imshow(JI_mat,cmap=cmap)
            plt.clim(0, 1)
            plt.colorbar()
            
            if show_axes:
                if x_ticks is None:
                    plt.xticks(list(range(sig_pathways_num)),sig_pathways,rotation=90,size=axes_size)
                    plt.yticks(list(range(sig_pathways_num)),sig_pathways,size=axes_size)
                else:
                    plt.xticks(list(range(sig_pathways_num)),[],rotation=90,size=axes_size)
                    plt.yticks(list(range(sig_pathways_num)),x_ticks,size=axes_size)
            else:
                plt.xticks(list(range(sig_pathways_num)),[],size=axes_size)
                plt.yticks(list(range(sig_pathways_num)),list(range(1,sig_pathways_num+1)),size=axes_size)

            if show_text:
                for i in range(JI_mat.shape[1]):
                    for j in range(JI_mat.shape[0]):
                        plt.text(j,i,np.round(JI_mat[i,j],2),ha='center',va='center')
        
        return JI_mat,sig_pathways
        
        


            
# if __name__ == "__main__":

    # preprocess data
    # old_adata=sc.read_h5ad('./data/HCA/lung_disease_filtered.h5ad')
    # adata=util.lung_disease_preprocess(old_adata,cell_type='capillary endothelial cell',n_sample=3000,seed=42,n_top_genes=2000)
    # print(adata)
    # adata.write_h5ad(f'./data/HCA/processed_lung_disease_filtered.h5ad')
    
    # adata=sc.read_h5ad('./data/HCA/processed_lung_disease_filtered.h5ad')
    # geneset_path='./data/pathway_db/c2.cp.reactome.v2025.1.Hs.symbols.gmt'

    # # Initialize model
    # model=scKnockPath()
    # X,Xh,y,groups_list=model.prepare_data(adata=adata,obs_y='disease',
    #                                     layer='lognorm',genesets=geneset_path,
    #                                     gene_names=adata.var['feature_name'],
    #                                     gene_thresh=15)

    # # generate knockoff data
    # file_name=f'norm_GMMknockoff3.npz'
    # data_file = os.path.join('./data/HCA', file_name)
    # load=True
    # np.random.seed(42)
    # if load and os.path.exists(data_file):
    #     print('loading data')
    #     data = np.load(data_file,allow_pickle=True)
    #     Xc = data['Xc']
    #     y = pd.Series(data['y'])
    # else:
    #     print('generating data')
    #     Xc,X,Xk = model.knockoff_sampler(X,Xh,method='GMM')
    #     np.savez(data_file, Xc=Xc, y=y)
    

    # lambdas = np.logspace(-1, -5, num=5)  # Generates log-spaced numbers 
    
    # out=model.fit_lambda(Xc,y,groups_list,lambdas,fdr=0.2,use_knockoff_plus=False)
    
    # print(out)

