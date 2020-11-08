"
=================================================
This file is used to generate data from specific |
distribution to simulate the TCA+TWAS pipepine. |
=================================================
"

  "
  Generate SNPs and Gene Expression Data.
  @Params:
      N:          #samples
      D:          #SNP sequence length
      MAF:        #minor allele frequency
      M:          #gene expressions
      K:          #cell types
      p_slap:     sparsity. Proportion of gene not related to GE
      sigma_beta:         std of Beta (mu = X*beta+epsilon_{mu},beta gaussian)
      sigma_mu:           std of mu_z (mu = X*beta+epsilon_{mu},epsilon_mu gaussian)
      sigma_z:            std of Z (Z^{i}_{km}|mu^{i}_{km},epsilon^{i}_{km}~N{mu^{i}_{km},epsilon^{i}_{km})
      sigma_g:            std of Gene Expression
      alpha:              alpha for Cell Type Proportion, which follows dirichlet distribution
      Q                   number of traits
      train_proportion:   the proportion of training dataset
  @Shape information:
      This part tell us the shape information before we split them into training data and test data
      X:          N*D Matrix SNPs
      Z:          N*K*M Array cell type specific gene expression
      G:          M*N Matrix mixed gene expressions
      beta:       D*K*M Array effect size of SNPs on Z
      W:          N*K Matrix cell type proportion
      T:          Q Traits
      c1:         N*3
  @Returns:
      traing_data: ['X'=X, 'Z'=Z, 'G'=G, 'beta'=beta, 'W'=W, 'T'=T, 'c1'=c1]
      test_data:   ['X'=X, 'Z'=Z, 'G'=G, 'beta'=beta, 'W'=W, 'T'=T, 'c1'=c1]
  "


generate_data <- function(N=100,
                          D=2000,
                          MAF=0.1,
                          M=10,
                          K=6,
                          p_slab=0.5,
                          sigma_beta=0.3,
                          sigma_mu=0.3,
                          sigma_z=0.3,
                          sigma_g=0.3,
                          alpha=0.3,
                          Q=5,
                          train_proportion=0.7,
                          male_proportion=0.5,
                          smoking_proportion=0.2,
                          pc_num=2
                          ) {
    library('MCMCpack')
    
    #### names
    id_name = sapply(1:N, function(x) paste('SAMPLE',toString(x),sep = ''))
    cell_type=sapply(1:K, function(x) paste('Cell_type',x,sep = ''))
    snp_name=sapply(1:D, function(x) paste('SNP',toString(x),sep = ''))
    ge_name=sapply(1:M, function(x) paste('GE',toString(x),sep=''))
                   
    #### generate X
    X = sapply(1:D,function(x) rbinom(N, 2,MAF))
    
    colnames(X, do.NULL = TRUE, prefix = "col")
    colnames(X) <- snp_name
    rownames(X, do.NULL = TRUE, prefix = "col")
    rownames(X) <- id_name
    
    # mask a normal to form spikeslab
    beta = rnorm(D*K*M, mean=0, sd=sigma_beta)
    mask = rbinom(D*K*M, 1, 1-p_slab)
    beta = array(beta*mask, c(D,K,M))
    
    # \mu = X*\beta+\epsilon_mu
    epsilon_mu = array(rnorm(N*K*M, mean = 0, sd=sigma_mu), c(N,K,M))
    mu_z = lapply(c(1:M),function(x) X%*%beta[,,x])
    mu_z = array(as.numeric(unlist(mu_z)), dim=c(N,K,M)) + epsilon_mu
    # Simulate mu_jh that is shared across individuals
    mu <- array(runif(M*K, min=0, max=10), c(M,K))
    mu_stacked <- aperm(abind(lapply(1:N, function(x) mu), along = 3), c(3,2,1))
    mu_z <- mu_z + mu_stacked
    
    #W = rdirichlet(N, rep(alpha,K))
    W <- rdirichlet(N, runif(K, 0, 1))
    colnames(W, do.NULL = TRUE, prefix = "col")
    colnames(W) <- cell_type
    rownames(W, do.NULL = TRUE, prefix = "col")
    rownames(W) <- id_name
  
    # Beta :QxM
    # T: NxQ
    # epsilon: NxQ
    
    ### generate C1 and C2
    male = matrix(rbinom(N, 1, 0.5),nrow=N)
    smoking = matrix(rbinom(N, 2, 0.2),nrow=N)
    smoking = (smoking - min(smoking))/(max(smoking))
    age =matrix(sapply(sapply(rnorm(N,50, 20), function(x) floor(x)), function(x)  if (x<0){x=20} else{x=x}),nrow=N) 
    age = (age-min(age))/(max(age)-min(age))
    c1=cbind(male, smoking, age)
    c2=matrix(rnorm(N*pc_num), nrow=N)
    rownames(c1, do.NULL = TRUE, prefix = "col")
    rownames(c1) <- id_name
    colnames(c1, do.NULL = TRUE, prefix = "col")
    colnames(c1) <- list('male','smoking','age')
    rownames(c2, do.NULL = TRUE, prefix = "col")
    rownames(c2) <- id_name
    colnames(c2, do.NULL = TRUE, prefix = "col")
    colnames(c2) <- sapply(1:pc_num, function(x) paste('PC',toString(x),sep = ''))
                           
    ### Generating Z
    # sample z from mu_z and epsilon_za and C1
    epsilon_z = array(rnorm(N*K*M,mean=0,sd=sigma_z), c(N,K,M))
    Z = mu_z + epsilon_z
    dimnames(Z)[[1]] <- id_name
    dimnames(Z)[[2]] <- cell_type
    dimnames(Z)[[3]] <- ge_name
    # sample gama
    gama = array(sapply(1:K*M, function(x) runif(3, 0, 1)), dim = c(K, 3, M))
    c_gama = array(rep(0, N*K*M), c(N,K,M))
    for (i in 1:K) {
      c_gama[,i,] = c1 %*% gama[i,,]
    }
    Z = Z + c_gama
                        
    ### Generate Gene expressions
    G = sapply(c(1:M),function(x) rowSums(W*Z[,,c(x)]))
    epsilon_G = replicate(M,rnorm(N,mean=0,sd=sigma_g))
    G = G+epsilon_G
    T_t = runif(1)
    epsilon_G = replicate(Q,rnorm(N))*T_t*T_t
    beta_t = matrix(runif(Q*M),ncol=M)
    T = G%*%t(beta_t)+epsilon_G
    #generate delta
    delta = matrix(sapply(1:M*pc_num, function(x) runif(pc_num,0,1)), nrow = pc_num)
    c_delta = c2%*%delta
    G = G + c_delta
    
    rownames(G, do.NULL = TRUE, prefix = "col")
    rownames(G) <- id_name
    colnames(G, do.NULL = TRUE, prefix = "col")
    colnames(G) <- ge_name
    G = t(G)                      
    ### split into train and test data
    prop = floor(train_proportion*dim(X)[1])
    train_data = list('X'=X[c(1:prop),],
                      'Z'=Z[c(1:prop),,],
                      'G'=G[,c(1:prop)],
                      'beta'=beta,
                      'W'=W[c(1:prop),],
                      'T'=T[c(1:prop),],
                      'C1'=c1[c(1:prop),],
                      'C2'=c2[c(1:prop),],
                      'MU'=mu)
    test_data = list('X'=X[-c(1:prop),],
                     'Z'=Z[-c(1:prop),,],
                     'G'=G[,-c(1:prop)],
                     'beta'=beta,
                     'W'=W[-c(1:prop),],
                     'T'=T[-c(1:prop),],
                     'C1'=c1[-c(1:prop),],
                     'C2'=c2[-c(1:prop),],
                     'MU'=mu)
    
    return(list('train_data'=train_data,'test_data'=test_data))
}
