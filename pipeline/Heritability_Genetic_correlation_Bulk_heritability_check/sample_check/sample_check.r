library(glmnet)
library(matrixStats)
library(abind)
library(pracma)
library(matrixcalc)
library(TCA)
library(MCMCpack)

summary_statistics <- function(mdl,train_X,test_X,train_c1,test_c1,
                               train_G,test_G,train_Z,test_Z,beta,
                               gamma_c1,gamma_c2,model=1){
    # TCA estimation & parameters
    Z_hat_tca = tensor(train_G,mdl,log_file=NULL,verbose=FALSE,debug=FALSE)
    rmse_tca = t(do.call(rbind,lapply(1:length(Z_hat_tca),
                                    function(x) sqrt(sum((Z_hat_tca[[x]]-train_Z[,x])^2)/length(train_Z[,x])))))
    Z_hat_cor_tca = sapply(1:length(Z_hat_tca), function(x) cor(unlist(Z_hat_tca[x]),train_Z[,x]))
    cell_type = ncol(mdl$W)
    if(model==1){
        tca_beta = t(matrix(mdl$gammas_hat[,grepl("_SNP",colnames(mdl$gammas_hat))],ncol=cell_type))
        beta_hat_cor_tca = lapply(1:cell_type,function(x) 
            cor(beta[x,],mdl$gammas_hat[,grepl(paste("Cell_type",x,".G",sep=''),colnames(mdl$gammas_hat))]))        
    }else{ # model 2 does not predict beta for snps
        tca_beta = NULL
        beta_hat_cor_tca = NULL
    }
    gamma_hat_cor_tca = cor(array(gamma_c1),
                        mdl$gammas_hat[,grepl('male|smok|age', colnames(mdl$gammas_hat))])
    if(model==1){
        pred = cbind(train_X,train_c1)
        test_pred = cbind(test_X,test_c1)        
    }else{ # model 2 does not use train_X/test_X as C1
        pred = train_c1
        test_pred = test_c1
    }
    # TCA parameter direct estimation <----> C1 & X
    Z_hat_train_second_eq_tca = sapply(1:cell_type,function(x) pred %*% mdl$gammas_hat[,((x-1)*dim(pred)[2]+1):(x*dim(pred)[2])])
    Z_hat_train_second_eq_cor_tca = diag(cor(train_Z,Z_hat_train_second_eq_tca))
    Z_hat_test_second_eq_tca = sapply(1:cell_type,function(x) test_pred %*% mdl$gammas_hat[,((x-1)*dim(test_pred)[2]+1):(x*dim(test_pred)[2])])
    Z_hat_test_second_eq_cor_tca = diag(cor(test_Z,Z_hat_test_second_eq_tca))
    
    # Lasso bulk data regression <----> X only
    glmnet.mdl.X.cv <- cv.glmnet(x=train_X,y=t(train_G),nfolds=5)
    glmnet.mdl.X <- glmnet(x=train_X,y=t(train_G),lambda=glmnet.mdl.X.cv$lambda.min)
    beta_full_X_bulk <- as.numeric(glmnet.mdl.X$beta)
    # extract non zero predictors and recorrelate
    predictors.X <- colnames(train_X)[which(beta_full_X_bulk!=0)]
    beta_X <- as.matrix(c(glmnet.mdl.X$a0,as.matrix(glmnet.mdl.X$beta[predictors.X,])))
    bias_one <- numeric(nrow(train_X))+1
    G_hat_train_lasso <- cbind(bias_one,train_X[,predictors.X]) %*% beta_X
    G_hat_train_cor_bulk <- cor(t(train_G),G_hat_train_lasso)
    G_hat_test_lasso <- cbind(numeric(nrow(test_X))+1,test_X[,predictors.X]) %*% beta_X
    G_hat_test_cor_bulk <- cor(t(test_G),G_hat_test_lasso)
    
    # cell type specific lasso
    corrs = numeric(ncol(mdl$W))
    corrs.real = numeric(ncol(mdl$W))
    corrs.beta = numeric(ncol(mdl$W))
    corrs.test.real = numeric(ncol(mdl$W))
    dev_ratio = numeric(ncol(mdl$W))
    rmse_lasso = numeric(ncol(mdl$W))
    beta_full_cell = matrix(0,nrow=ncol(mdl$W),ncol=dim(train_X)[2])
    Z_hat_test_lasso = matrix(0,nrow=nrow(test_Z),ncol=ncol(test_Z))
    for (h in 1:ncol(mdl$W)){
        
        glmnet.mdl.cv <- cv.glmnet(x=train_X,y=Z_hat_tca[[h]],standardize=FALSE,alpha=1,nfolds=5)
        rmse_lasso[h] = sqrt(glmnet.mdl.cv$cvm[glmnet.mdl.cv$lambda == glmnet.mdl.cv$lambda.min])
        glmnet.mdl <- glmnet(x=train_X,y=Z_hat_tca[[h]],standardize=FALSE,alpha=1,lambda=glmnet.mdl.cv$lambda.min)
        dev_ratio[h] <- glmnet.mdl$dev.ratio
        beta.full <- as.numeric(glmnet.mdl$beta)
        beta_full_cell[h,] <- beta.full
        predictors <- colnames(train_X)[which(beta.full != 0)]
        beta_lasso <- as.matrix(c(glmnet.mdl$a0,as.matrix(glmnet.mdl$beta[predictors,])))
        Z_hat_train_lasso <- cbind(numeric(nrow(train_X))+1,train_X[,predictors]) %*% beta_lasso
        Z_hat_test_lasso[,h] <- cbind(numeric(nrow(test_X))+1,test_X[,predictors]) %*% beta_lasso
        Z_hat_test_lasso[,h] = signif(Z_hat_test_lasso[,h], digits = 5)
        if(sum(beta_lasso)==0 | sd(Z_hat_train_lasso)==0){
            # model 2 lasso forces SNPs effect to become zero.
            # If there were no predictors, then cor would be NaN. Which also means
            # no snps is correlated with Z. cor just set to 0.
            corrs[h] = 0
            corrs.real[h] = 0
            corrs.beta[h] = 0
            corrs.test.real[h] = 0
        }else{
            corrs[h] <- cor(t(Z_hat_tca[[h]]),Z_hat_train_lasso)
            corrs.real[h] <- cor(train_Z[,h],Z_hat_train_lasso)
            corrs.beta[h] <- cor(beta[h,],beta.full)
            corrs.test.real[h] <- cor(test_Z[,h],Z_hat_test_lasso[,h])
        }
    }
    # Precision & Recall for lasso
    Binary_True=(beta!=0)
    Binary_Pred_Lasso=(beta_full_cell!=0)
    TP = do.call(rbind,lapply(1:cell_type,function(x) sum(as.integer(Binary_True[x,]&Binary_Pred_Lasso[x,]))))
    FP = do.call(rbind,lapply(1:cell_type,function(x) sum(as.integer((!Binary_True[x,])&Binary_Pred_Lasso[x,]))))
    FN = do.call(rbind,lapply(1:cell_type,function(x) sum(as.integer(Binary_True[x,]&(!Binary_Pred_Lasso[x,])))))
    precision_lasso = TP/(TP+FP)
    recall_lasso = TP/(TP+FN)
    
    return(list('Z_hat_test_lasso'=Z_hat_test_lasso,'Z_hat_tca'=Z_hat_tca,'Z_hat_cor_tca'=Z_hat_cor_tca,'beta_hat_cor_tca'=beta_hat_cor_tca,'gamma_hat_cor_tca'=gamma_hat_cor_tca,
               'Z_hat_train_second_eq_cor_tca'=Z_hat_train_second_eq_cor_tca,'Z_hat_test_second_eq_cor_tca'=Z_hat_test_second_eq_cor_tca,
               'G_hat_train_cor_bulk'=G_hat_train_cor_bulk,'G_hat_test_cor_bulk'=G_hat_test_cor_bulk,
               'cor_lasso_tca'=corrs,'cor_lasso_real_train'=corrs.real,'cor_lasso_real_test'=corrs.test.real,
               'cor_beta_lasso_real'=corrs.beta,'beta_hat_lasso'=beta_full_cell,'dev_rat'=dev_ratio,
               'precision_lasso'=precision_lasso,'recall_lasso'=recall_lasso,'rmse_lasso'=rmse_lasso,'rmse_tca'=rmse_tca))
}


generate_params <- function(cell_her=TRUE,gene_cor=FALSE,bulk_her=FALSE,seed=1,N=5000,M=100,herr_arr=c(),D=20,pslab=0.5,her=0.05,her_bulk=0.8){
    set.seed(1)
    if(length(herr_arr)!=0){
        M = length(herr_arr)
    }else{
        M = M
    }
    K = 4
    N = N
    pc_num = 2
    # cis_snps_nums = floor(runif(M, min = 150, max = 400))
    D = D
    cis_snps_nums = array(D,M)
    sigma_g = 0.01
    sigma_z = 0.1
    
    if(gene_cor){ # varing genetic correlation across genes
        corr_seq = seq(from=0,to=pslab,length.out=M)
    }else{
        corr_seq = seq(from=0,to=0,length.out=M)
    }
    pslab = matrix(pslab,nrow=M,ncol=K)
    corr_matrix = lapply(1:M,function(x) matrix(corr_seq[x],nrow=K,ncol=K))
    for(i in 1:M){
        diag(corr_matrix[[i]])=1
    }
    
    if(cell_her){ # varing heribility across genes
        heritibility_cell_specific = do.call(rbind,lapply(1:K,function(x) herr_arr))
    }else{
        heritibility_cell_specific = matrix(her,nrow=K,ncol=M)
    }
    
    if(bulk_her){ # varing bulk heritibility across genes
        heritibility_bulk = seq(from=0.1,to=her_bulk,length.out=M)
    }else{
        heritibility_bulk = array(her_bulk,M)
    }
        
    MAF = lapply(1:length(cis_snps_nums),function(x) runif(cis_snps_nums[x],min=0.1,max=0.5))
    id_name = sapply(1:N, function(x) paste('SAMPLE',toString(x),sep = ''))
    cell_type_name=sapply(1:K, function(x) paste('Cell_type',x,sep = ''))
    
    ### generate C1 and C2
    male = matrix(rbinom(N, 1, 0.5),nrow=N)
    smoking = matrix(rbinom(N, 2, 0.2),nrow=N)
    smoking = (smoking - min(smoking))/(max(smoking))
    age = matrix(sapply(sapply(rnorm(N,50,20), function(x) floor(x)), function(x)  if (x<0){x=20} else{x=x}),nrow=N) 
    age = (age-min(age))/(max(age)-min(age))
    c1 = scale(cbind(male, smoking, age))
    c2 = scale(matrix(rnorm(N*pc_num), nrow=N))
                            
    rownames(c1) = id_name
    colnames(c1) = c('male','smoking','age')
    rownames(c2) = id_name
    colnames(c2) = sapply(1:pc_num, function(x) paste('PC',toString(x),sep = ''))
                          
    p1 = dim(c1)[2]
    sigma_gamma = sqrt((1-heritibility_cell_specific[1,]-sigma_z^2)/p1)
    var_beta = heritibility_cell_specific[1,]*(p1*sigma_gamma^2+sigma_z^2)/(1-heritibility_cell_specific[1,])/D
    
    # cell type proportion, from the real data estimated
    W_alpha = c(26.553683792256,17.6621467979005,4.48671525658667,1.56874856517803)#,0.178797663350679)
    #W_alpha = 50.4500920752719
    W_xsi = c(0.526335685426257,0.350091468050216,0.088933737720289,0.0310950585152043,0.00354405028803339)
    if(K <= length(W_alpha)){
        W = rdirichlet(N, W_alpha[1:K])
    }else{ #more cell type prop
        W = rdirichlet(N, runif(K, 0, 1))
    }
    colnames(W) = cell_type_name
    rownames(W) = id_name
    
    alpha = W_alpha
    alpha_0 = sum(alpha)
    alpha_tilde = alpha/alpha_0
    m2_alpha = alpha_tilde %*% t(alpha_tilde)*alpha_0/(alpha_0+1)
    diag(m2_alpha) = alpha_tilde*(1-alpha_tilde)/(alpha_0+1)+alpha_tilde^2
    
    return(list(
        'M'=M,'K'=K,'N'=N,'pc_num'=pc_num,'D'=D,'id_name'=id_name,'cell_name'=cell_type_name,
        'her'=heritibility_cell_specific,'her_bulk'=heritibility_bulk,'pslab'=pslab,
        'MAF'=MAF,'c1'=c1,'c2'=c2,'W'=W,'sigma_gamma'=sigma_gamma,'sigma_g'=sigma_g,'m2_alpha'=m2_alpha,
        'sigma_z'=sigma_z,'beta_cor'=corr_matrix,'var_beta'=var_beta,'W_alpha'=W_alpha[1:K]))
}

one_gene <- function(par,g,seed=1){
    set.seed(seed)
    cell_type = par$K
    n_snps = par$D
    her = par$her[,g]
    cor = par$beta_cor[[g]]
    pslab = par$pslab[g,]
    var_beta = array(par$var_beta[g],par$K)
    
    maf = par$MAF[[g]]
    dummy_var = 10
    
    var_matrix=sqrt(var_beta%*%t(var_beta))
    Sigma_beta = cor*var_matrix/(1-pslab)^2
    diag(Sigma_beta) = diag(Sigma_beta) *(1-pslab)
    beta_ = mvrnorm(n_snps,mu=rep(0,par$K),Sigma=Sigma_beta,tol=1e-4,empirical=TRUE)
    rmask = do.call(rbind,lapply(1:cell_type,function(x) rbinom(n_snps,1,1-pslab[x])))
    beta = t(beta_)*rmask

    X = t(do.call(rbind,lapply(1:n_snps,function(x) rbinom(par$N,2,maf[x]))))
    colnames(X) = sapply(1:n_snps, function(x) paste('G',g,'_SNP',toString(x),sep = ''))
    rownames(X) = par$id_name
    X = scale(X)
    # enforce variance of samples by adjusting sd of beta
    xbeta_var = n_snps*var_beta[1]
    beta = do.call(rbind,
                      lapply(1:cell_type,function(x) beta[x,]*(sqrt(xbeta_var)/sd(X %*% beta[x,]))))

    mu_z = X%*%t(beta)
    epsilon_z = array(rnorm(par$N*cell_type,mean=0,sd=par$sigma_z), c(par$N,cell_type))
    epsilon_z = apply(epsilon_z,2,function(x) x*par$sigma_z/sd(x))
    gamma_c1 = do.call(rbind,lapply(1:cell_type,function(x) rnorm(dim(par$c1)[2],mean=0,sd=par$sigma_gamma)))     # enforce variance of samples by adjusting sd of gamma_c1
    # enforce the c1_gamma variance
    c1gamma_var = dim(par$c1)[2]*(par$sigma_gamma[g])^2
    gamma_c1 = do.call(rbind,
                       lapply(1:cell_type,function(x) gamma_c1[x,]*(sqrt(c1gamma_var)/sd(par$c1 %*% gamma_c1[x,]))))
    c1_gamma = par$c1 %*% t(gamma_c1)

    Z = epsilon_z + mu_z + c1_gamma
    rownames(Z) = par$id_name
    colnames(Z) = par$cell_name
    
    ### Generate Gene expressions
    G = rowSums((par$W)*Z)
    epsilon_G = rnorm(par$N,mean=0,sd=par$sigma_g)
    epsilon_G = epsilon_G*par$sigma_g/sd(epsilon_G)
    
    # use bulk level heritability to calculate sigma_gamma
    bulk_her_nom = sum(hadamard.prod(par$m2_alpha,var(mu_z)))
    bulk_her_c1_gamma = sum(hadamard.prod(par$m2_alpha,var(c1_gamma)))
    bulk_her_epsilon_z = sum(hadamard.prod(par$m2_alpha,var(epsilon_z)))
    bulk_her_z = sum(hadamard.prod(par$m2_alpha,var(Z)))
    
    bulk_her_denom = sum(hadamard.prod(par$m2_alpha,var(Z)))+var(epsilon_G)
    bulk_her_val = par$her_bulk[g]*her[1]
    sd_gamma_c2 = sqrt((bulk_her_nom/bulk_her_val-bulk_her_denom)/par$pc_num)
    
    gamma_c2 = matrix(rnorm(par$pc_num,mean=0,sd=sd_gamma_c2),nrow=par$pc_num,ncol=1)
    c2gamma_var = par$pc_num * (sd_gamma_c2)^2
    gamma_c2 = gamma_c2*(sqrt(c2gamma_var)/sd(par$c2%*%gamma_c2))
    c2_gamma = par$c2 %*% gamma_c2
    
    G = t(G + epsilon_G + c2_gamma)
    real_denom = var(t(G))
    colnames(G) <- par$id_name
    rownames(G) <- paste('gene',g,sep='_')
    G = as.data.frame(G)
    
    real_bulk_her = bulk_her_nom/(bulk_her_denom+var(c2_gamma))
                              
    return(list('X'=X,'beta'=beta,'gamma_c1'=t(gamma_c1),'bulk_her'=real_bulk_her,
               'gamma_c2'=gamma_c2,'Z'=Z,'G'=G,'epsilon_z'=epsilon_z))
}
#data = one_gene(par=params,1)

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}
summary_statistics_z_hat_test_lasso <- function(mdl,train_X,test_X,train_c1,test_c1,
                               train_G,test_G,train_Z,test_Z,beta,
                               gamma_c1,gamma_c2,model=1){
    
    Z_hat_tca = tensor(train_G,mdl,log_file=NULL,verbose=FALSE,debug=FALSE)
    Z_hat_test_lasso = matrix(0,nrow=nrow(test_Z),ncol=ncol(test_Z))
    beta_full_cell = matrix(0,nrow=ncol(mdl$W),ncol=dim(train_X)[2])
    for (h in 1:ncol(mdl$W)){
        glmnet.mdl.cv <- cv.glmnet(x=train_X,y=Z_hat_tca[[h]],standardize=FALSE,alpha=1,nfolds=5)
        glmnet.mdl <- glmnet(x=train_X,y=Z_hat_tca[[h]],standardize=FALSE,alpha=1,lambda=glmnet.mdl.cv$lambda.min)
        beta.full <- as.numeric(glmnet.mdl$beta)
        beta_full_cell[h,] <- beta.full
        predictors <- colnames(train_X)[which(beta.full != 0)]
        beta_lasso <- as.matrix(c(glmnet.mdl$a0,as.matrix(glmnet.mdl$beta[predictors,])))
        Z_hat_test_lasso[,h] <- cbind(numeric(nrow(test_X))+1,test_X[,predictors]) %*% beta_lasso
        Z_hat_test_lasso[,h] = signif(Z_hat_test_lasso[,h], digits = 5)
    }
    
    return(list('Z_hat_test_lasso'=Z_hat_test_lasso,'beta_hat_lasso'=beta_full_cell))
}

#params = generate_params(cont_her='other',gene_cor=FALSE,seed=1,N=105000,M=1000,D=20,pslab=0.5,her=0.05)
params_sample_check = generate_params(cell_her=FALSE,gene_cor=FALSE,bulk_her=FALSE,
                         seed=1,N=105000,M=1000,
                         D=20,pslab=0.5,her=0.05,
                         her_bulk=0.8)
data_sample_check = lapply(1:params_sample_check$M,function(x) one_gene(par=params_sample_check,x))

sample_size_arr = c(100,200,300,400,500)
GWAS_size = 100000
li_sample_check_power = list()
p_value_matrix = array(0,c(length(sample_size_arr),length(data_sample_check),params_sample_check$K))
for(sample_size_i in 1:length(sample_size_arr)){
    sample_size = sample_size_arr[sample_size_i]
    for(g in 1:length(data_sample_check)){
        train_X = (data_sample_check[[g]]$X)[1:sample_size,]
        test_X = data_sample_check[[g]]$X[(params_sample_check$N-GWAS_size+1):params_sample_check$N,]
        train_W = params_sample_check$W[1:sample_size,]
        test_W = params_sample_check$W[(params_sample_check$N-GWAS_size+1):params_sample_check$N,]
        train_G = data_sample_check[[g]]$G[,1:sample_size]
        test_G = data_sample_check[[g]]$G[,(params_sample_check$N-GWAS_size+1):params_sample_check$N]
        train_c1 = params_sample_check$c1[1:sample_size,]
        test_c1 = params_sample_check$c1[(params_sample_check$N-GWAS_size+1):params_sample_check$N,]
        train_c2 = params_sample_check$c2[1:sample_size,]
        test_c2 = params_sample_check$c2[(params_sample_check$N-GWAS_size+1):params_sample_check$N,]
        train_Z = data_sample_check[[g]]$Z[1:sample_size,]
        test_Z = data_sample_check[[g]]$Z[(params_sample_check$N-GWAS_size+1):params_sample_check$N,]
        beta = data_sample_check[[g]]$beta
        gamma_c1 = data_sample_check[[g]]$gamma_c1
        gamma_c2 = data_sample_check[[g]]$gamma_c2
        tca.mdl1 = tca(X=train_G,W=train_W,C1=cbind(train_X,train_c1),verbose=FALSE)
        tca.mdl1.summary = summary_statistics_z_hat_test_lasso(mdl=tca.mdl1,train_X=train_X,test_X=test_X,
                                                  train_c1=train_c1,test_c1=test_c1,train_G=train_G,
                                                  test_G=test_G,train_Z=train_Z,test_Z=test_Z,
                                                  beta=beta,gamma_c1=gamma_c1,gamma_c2=gamma_c2,model=1)
        li_sample_check_power[[paste(g,'.',sample_size,sep='')]] = tca.mdl1.summary
        r_testz = dim(test_Z)[1]
        c_testz = dim(test_Z)[2]
        Y = test_Z + matrix(rnorm(r_testz*c_testz,mean=0,sd=sqrt(499)),nrow=r_testz)
        for(cell_i in 1:params_sample_check$K){
            if(length(unique(tca.mdl1.summary$Z_hat_test_lasso[,cell_i]))==1){
                p_value_matrix[sample_size_i,g,cell_i] = 1
                next #If lasso predict everything to be the same. p value breaks
            }
            lm.model <- lm(Y[,cell_i]~tca.mdl1.summary$Z_hat_test_lasso[,cell_i])
            p_value_matrix[sample_size_i,g,cell_i] = lmp(lm.model)            
        }
    }
    print(paste('sample size',sample_size,'done'))
}
save(li_sample_check_power,sample_size_arr,params_sample_check,data_sample_check,p_value_matrix,file = "power_check.RData")

