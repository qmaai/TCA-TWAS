<<<<<<< HEAD

library("glmnet")
lasso_on_one_gene <- function(Z,SNPs,cell_names){
    glmnet.cv = lapply(1:length(cell_names),function(x) cv.glmnet(x=SNPs,y=t(Z[x,]),nfold=5))
    betas = lapply(1:length(cell_names),function(x) as.matrix(t(coef(glmnet.cv[[x]],s="lambda.min"))))
    betas = as.data.frame(do.call(rbind, betas))
    rownames(betas) = cell_names
=======
.libPaths(c(.libPaths(), "/u/home/q/qmaai/.R"))
shhh <- suppressPackageStartupMessages
shhh(library(glmnet))
shhh(library(optparse))

report_statistics <- function(glmnet.cv,pred_num){
    lambda = c(glmnet.cv$lambda.1se,glmnet.cv$lambda.min)
    ind = which(glmnet.cv$lambda==lambda[1] | glmnet.cv$lambda==lambda[2])
    df = glmnet.cv$nzero[ind]
    result = cbind(lambda,df,pred_num)
    rownames(result) = c('1se','min')
    return(result)
}

lasso_on_one_gene <- function(Z,SNPs,cell_names,i,log_info=FALSE){
    glmnet.cv = lapply(1:length(cell_names),function(x) cv.glmnet(x=SNPs,y=t(Z[x,]),nfold=5))
    betas = lapply(1:length(cell_names),function(x) as.matrix(t(coef(glmnet.cv[[x]],s="lambda.1se"))))
    betas = as.data.frame(do.call(rbind, betas))
    rownames(betas) = cell_names
    if(log_info){
        print(paste('gene',i,'statistics',sep=' '))
        info = lapply(1:length(cell_names),
                      function(x) report_statistics(glmnet.cv[[x]],dim(SNPs)[2]))
        print(info)
    }
>>>>>>> bd8fe80edc5f588f001e3a0cafad85a1a93b8d14
    return(betas)
}

# each gene has its own SNPs. So do lasso separately
main <- function(opt){
    load(opt$z_path)
    load(opt$mdl_path)
    num_of_gene = length(z_result)
<<<<<<< HEAD
    beta_chunk = lapply(1:num_of_gene,function(x) lasso_on_one_gene(z_result[[x]],
                                                       model_result[[x]]$C1,
                                                       colnames(model_result[[x]]$mus_hat)))
=======
    beta_chunk = lapply(1:num_of_gene,
                        function(x) lasso_on_one_gene(z_result[[x]],
                                                       model_result[[x]]$C1,
                                                       colnames(model_result[[x]]$mus_hat),
                                                       x,
                                                       log_info=opt$verbose))
>>>>>>> bd8fe80edc5f588f001e3a0cafad85a1a93b8d14
    
    if(opt$output_path!='NULL'){
        dir.create(opt$output_path,showWarnings = FALSE)
        opt$end = opt$start + num_of_gene
        if(opt$start==0){
<<<<<<< HEAD
            matches = unlist(regmatches(path, gregexpr("[[:digit:]]+", path)))
            opt$start = matches[1]
            opt$end = matches[2]
        }
        save(beta_chunk,file=paste(opt$output_path,'lasso_beta',opt$start,'_',opt$end,'.RData'))
=======
            matches = unlist(regmatches(opt$z_path, gregexpr("[[:digit:]]+", opt$z_path)))
            opt$start = matches[1]
            opt$end = matches[2]
        }
        file_path = paste(opt$output_path,'lasso_beta_',opt$start,'_',opt$end,'.RData',sep='')
        save(beta_chunk,file=file_path)
>>>>>>> bd8fe80edc5f588f001e3a0cafad85a1a93b8d14
    }
}

option_list = list(
    make_option("--z_path", type="character", 
              help="TCA prediction file path for a chunk of genes", metavar="character"),
    make_option("--mdl_path",type="character",
               help="TCA model file path for a chunk of genes",metavar="character"),
    make_option("--start", type="integer",default=0, 
              help="first gene in the chunk's order", metavar="character"),
<<<<<<< HEAD
    make_option("--output_path",type="character",
                default="/u/home/q/qmaai/project-sriram/tca-twas/data/tca_result/lasso_result/",
=======
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                help="Print extra output [default]"),
    make_option("--output_path",type="character",
                default="/u/home/q/qmaai/project-sriram/tca-twas/data/lasso_result/",
>>>>>>> bd8fe80edc5f588f001e3a0cafad85a1a93b8d14
               help='lasso output save position',metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (!interactive()){
    main(opt);
}
