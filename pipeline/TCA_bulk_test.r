
library(snpStats)
library(TCA)
library(pracma)
library(matrixcalc)
library(optparse)
library(pbapply)
library(tidyr)

find_cis_snps_gene <-function(S_geno,G_locs,S_map,G,flank = 1000000){
    G_locs = G_locs[G,]
    S_map = S_map[S_map$chromosome==G_locs$CHROM,]
    cis_snps = rownames(S_map[((S_map$position>G_locs$START-flank) & (S_map$position<G_locs$START))
                              | ((S_map$position<G_locs$END+flank) & (S_map$position>G_locs$END)),])
    return(S_geno[,colnames(S_geno) %in% cis_snps])
}

tca_per_gene <- function(snps,geno_locs,ge,prop,i){
    k = ncol(prop)
    n = nrow(prop)
    W_norms <- rowSums(prop**2)**0.5
    gene = rownames(ge)[i]
    X_raw = find_cis_snps_gene(snps$genotype,geno_locs,snps$map,gene)
    X = as(X_raw,'numeric')
    W = prop
    G = ge[i,]
    
    # find bad individual snps
    bad_rows = array(0,ncol(X))
    for (j in 1:ncol(X)){
        X_ = matrix(X[,j],nrow=nrow(X),ncol=1)
        C1_ = hadamard.prod(Reshape(Reshape(apply(W, 2, function(v) repmat(v,1,1)), n*k,1), n,k), repmat(X_, 1, k))
        D = as.matrix(cbind(W,C1_))
        if(qr(D)$rank<dim(D)[2]){
            bad_rows[j] = 1
        }
    }
    bad_rows = which(bad_rows==1)
    bad_snps = colnames(X[,c(bad_rows)])                                               
    X = X[,!colnames(X) %in% bad_snps]
    
    # remove across column colinearity
    rnm = matrix(rnorm(ncol(X)*nrow(X),mean=0,sd=0.1),ncol=ncol(X),nrow=nrow(X))
    
    tryCatch({
        tca.mdl = tca(X=G,W=W,C1=X+rnm,verbose=FALSE)
        z = tensor(G,tca.mdl,log_file=NULL,verbose=FALSE)
        
        # parse gammas into matrix
        gammas = tca.mdl$gammas_hat
        rn = data.frame(colnames(gammas))
        colnames(rn) = 'mixed_name'
        rn = rn %>% separate('mixed_name',c('cell','gene'))
        gammas = data.frame(matrix(gammas,nrow=dim(tca.mdl$mus_hat)[2],byrow=TRUE),
                            row.names=unique(rn[,1]))
        colnames(gammas) = unique(rn[,2])
        
        return(list('z'=do.call(rbind.data.frame,z),
                    'mdl'=list('mus_hat'=tca.mdl$mus_hat,
                              'sigmas_hat'=tca.mdl$sigmas_hat,
                              'tau_hat'=tca.mdl$tau_hat,
                              'deltas_hat'=tca.mdl$deltas_hat,
                              'gammas_hat'=gammas)))
    },error=function(e){
        return(list('z'=0,'mdl'=0))})
}
                                                  
# snp_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_filted_individual_ld_pruned'
# geno_locs_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_gene_locs.csv'
# GE_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_expression_without_rbc.csv'
# prop_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_props.csv'

main <- function(opt){   
    # loading data
    snps = read.plink(opt$snp_path)
    geno_locs = read.csv(opt$geno_locs_path,row.names=1)
    ge = read.csv(opt$gene_expression_path,row.names=1)
    colnames(ge) = substr(colnames(ge),2,10)
    prop = read.csv(opt$cell_type_proportion_path,row.names=1)
    
    # perform tca
    start = opt$start_gene
    end = opt$start_gene + opt$bulk_size
    tca_result = pblapply(start:end,
            function(x) tca_per_gene(snps,geno_locs,ge,prop,x))
    z_result = list()
    model_result = list()
    for(i in 1:length(tca_result)){
        z_result[i] = tca_result[[i]]$z
        model_result[i] = tca_result[[i]]$mdl
    }
    
    # saving result
    if(!opt$tca_mdl_path=='NULL'){
        dir.create(opt$tca_mdl_path,showWarnings = FALSE)
        save(z_result,file=paste(opt$tca_mdl_path,'z_result',start,'_',end,'.RData'))
        save(model_result,file=paste(opt$tca_mdl_path,'mdl_result',start,'_',end,'.RData'))
    }
}
                          
option_list = list(
    make_option("--snp_path", type="character",
              default="/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_filted_individual_ld_pruned", 
              help="snps file path", metavar="character"),
    make_option("--geno_locs_path",type="character",
               default="/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_gene_locs.csv",
               help="gene location file path",metavar="character"),
    make_option("--gene_expression_path", type="character",
              default="/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_expression_without_rbc.csv", 
              help="gene expression path", metavar="character"),
    make_option("--cell_type_proportion_path", type="character",
              default="/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_props.csv", 
              help="cell type proportion path", metavar="character"),
    make_option("--start_gene_pose",type="integer",default=1,
               help="the start position of the gene to analyse",metavar="number"),
    make_option("--bulk_size", type="integer", default=100, 
              help="number of genes to analyze per bulk", metavar="number"),
    make_option("--tca_mdl_path",type="character",
                default="/u/home/q/qmaai/project-sriram/tca-twas/data/tca_result/",
               help='tca output save position',metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (!interactive()){
    main(opt);
}

# snp_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_filted_individual_ld_pruned'
# geno_locs_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_gene_locs.csv'
# GE_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_expression_without_rbc.csv'
# prop_path = '/u/home/q/qmaai/project-sriram/tca-twas/data/pruned_data/dutch_props.csv'

# snps = read.plink(snp_path)
# geno_locs = read.csv(geno_locs_path,row.names=1)
# ge = read.csv(GE_path,row.names=1)
# colnames(ge) = substr(colnames(ge),2,10)
# prop = read.csv(prop_path,row.names=1)
# tca_result = pblapply(1:2,function(x) tca_per_gene(snps,geno_locs,ge,prop,x))
