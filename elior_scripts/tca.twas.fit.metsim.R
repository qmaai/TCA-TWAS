## This script extracts for each gene that was found associated with a phenotype (and reported by the output of tcareg.fit.metsim.sum.R) the set of genotypes in the gene and within 100KB upstream/downstream distance.
## Then, the script fits an elastic net model for predicting the cell-type-specific expression levels (estimated by TCA) of each of these genes based on their cis-SNPs.


library("data.table")
library("TCA")

datadir1 <- "/u/project/eeskin/halperin/ehsafe/METSIM/"
datadir2 <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/"

load_data <- function(filepath){
	X <- data.frame(fread(file  = filepath), row.names=1)
	col_names <- X[1,]
	X <- X[2:nrow(X),]
	colnames(X) <- col_names
	return (X)	
}


flank <- 1000000

# Extract the chromosome numbers and start, end positions in bp for each of the significant genes
pheno_names <- c("BMI","S_tottg","matsuda","B_GHbA1C")
genes <- character(0)
for (i in 1:length(pheno_names)){
	genes <- append(genes,rownames(read.table(file = paste("tcareg.metsim.",pheno_names[i],".marginal.significant_genes.pvals.txt",sep=""), header = TRUE, sep=",", row.names = 1)))
}
genes <- unique(genes)
info <- read.table(file = paste(datadir1,"gencode_v25lift37_gene_info.txt",sep=""), header = TRUE)
rownames(info) <- info$gene_id
genes.info <- info[genes,c("seqname","start","end")]


# Read the genotypes data
X <- fread(file = "/u/project/eeskin/halperin/ehsafe/FromPaivi/METSIM.FULL/334_subset/metsim_omni.334_subset.raw")
ids <- as.character(X$IID)
X <- as.matrix(X[,7:ncol(X)])
rownames(X) <- ids
X <- t(X)

# remove snps with any na values
X <- X[(rowSums(is.na(X)) == 0),]
snps <- rownames(X)

# extract the chr, position information of each snp in X and find the list of SNPs that are within 'flank' distance from at least one of the genes in 'genes'
check_snp <- function(G, snps.map, snp, flank){
	j <- gregexpr(pattern ='_',snp)
	s <- substr(snp, 1, j[[1]][1]-1)
	chr <- snps.map[s,1]
	pos <- snps.map[s,2]
	if (sum(G[,1] == chr & (G[,2] < pos & G[,3] > pos | abs(G[,2]-pos) < flank | abs(G[,3]-pos) < flank ))) TRUE else FALSE	
}

bim <- fread("/u/project/eeskin/halperin/ehsafe/FromPaivi/METSIM.FULL/metsim_omni.bim")
snps.map <- cbind(bim$V1,bim$V4)
rownames(snps.map) <- bim$V2
colnames(snps.map) <- c("chr","pos")


library("tictoc")
library("pbapply")
library("parallel")
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("flank","genes.info","snps.map","snps","check_snp"), envir=environment())
op <- pboptions(nout = 10, type = "none")
tic()
res <- pblapply(1:length(snps),function(j) check_snp(genes.info, snps.map, snps[j], flank), cl = cl)
# without parallel:
# res <- lapply(1:length(snps),function(j) check_snp(genes.info, snps.map, snps[j], flank))
toc()


# extract only the relevant SNPs
a <- c(as.logical(res))
Z <- X[a,]
	



library("glmnet")
library("matrixStats")
load("tca.mdl.metsim.Rdata")

X <- load_data(paste(datadir1,"Kal_v25_TPM_gene_filtered.csv",sep=""))
phenos <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335.txt",sep="")), row.names=1)
# Remove sample 2812 for which we have NA value in the Matsuda index
keep <- setdiff(1:nrow(phenos),which(rownames(phenos) == "2812"))
X <- X[,keep]

corrs <- matrix(0,nrow(genes.info),ncol(tca.mdl$W))
corrs.bulk <- numeric(nrow(genes.info))

s <- rownames(Z)
# Extract the genotypes to use for each gene
for (i in 1:nrow(genes.info)){
	print(sprintf("----- %d out of %s",i,nrow(genes.info)))
	tic()
	# get the cis-SNPs for the current gene
	clusterExport(cl, c("flank","genes.info","snps.map","s","check_snp","i"), envir=environment())
	r <- pblapply(1:length(s),function(j) check_snp(genes.info[i,], snps.map, s[j], flank), cl = cl)
	G <- Z[as.logical(r),]
	if (!nrow(G)){
		print(sprintf("num of SNPs for the glmnet model: %d",nrow(G))) 
		toc()
		next
	}
	# remove constant sites
	G <- G[setdiff(1:nrow(G),which(rowVars(G) == 0)),]
	print(sprintf("num of SNPs for the glmnet model: %d",nrow(G))) 
	#write.csv(Y, file = paste(datadir1,"cis_genotypes_per_gene/for_elastic_net_models/",rownames(genes.info[i,]),".cis_genotypes.for_elastic_net.txt",sep=""), quote = FALSE)

	# Fit a lasso model for each gene in the set of hits
	gene_name <- rownames(genes.info)[i]
	tca.mdl.gene <- tcasub(tca.mdl, gene_name, log_file = NULL, verbose = FALSE)
	Z_hat <- tensor(X[gene_name,], tca.mdl.gene, log_file = NULL, verbose = FALSE)
	# standardize G and save mean and std of this gene across all genotypes
	sds <- as.matrix(rowVars(as.matrix(G))**0.5)
	means <- as.matrix(rowMeans(as.matrix(G)))
	rownames(sds) <- rownames(G)
	rownames(means) <- rownames(G)
	G.scaled <- t(scale(t(G), center = TRUE, scale = TRUE))
	# prediction using the bulk data only - as a baseline
	glmnet.mdl.X.cv <- cv.glmnet(x = t(G.scaled), y = t(X[gene_name,]), standardize = FALSE, alpha=1, nfolds = 10)
	glmnet.mdl.X <- glmnet(x = t(G.scaled), y = t(X[gene_name,]), standardize = FALSE, alpha=1, lambda = glmnet.mdl.X.cv$lambda.min)
	beta.full.X <- as.numeric(glmnet.mdl.X$beta)
	predictors.X <- rownames(G)[which(beta.full.X != 0)]
	beta.X <- as.matrix(c(glmnet.mdl.X$a0,as.matrix(glmnet.mdl.X$beta[predictors.X,])))
	rownames(beta.X) <- c("Intercept",predictors.X)
	colnames(beta.X) <- gene_name
	write.csv(beta.X, file = paste(datadir2,"elastic_net_models/",rownames(genes.info[i,]),".elastic_net.bulk.predictors.txt",sep=""), quote = FALSE)
	X.pred <- cbind(numeric(ncol(G))+1,t(G.scaled[predictors.X,,drop=FALSE])) %*% beta.X
	if (sd(X.pred) > 0.00001){
		corrs.bulk[i] <- cor(t(X[gene_name,]),X.pred)
		print(sprintf("in-smaple correlation for gene %s with bulk: %f", gene_name, corrs.bulk[i]))
	}
	write.csv(means[predictors.X,], file = paste(datadir2,"elastic_net_models/",rownames(genes.info[i,]),".elastic_net.bulk.predictors.means.txt",sep=""), quote = FALSE)
	write.csv(sds[predictors.X,], file = paste(datadir2,"elastic_net_models/",rownames(genes.info[i,]),".elastic_net.bulk.predictors.sds.txt",sep=""), quote = FALSE)
	for (h in 1:ncol(tca.mdl$W)){
		cell_type <- colnames(tca.mdl$W)[h]
		glmnet.mdl.cv <- cv.glmnet(x = t(G.scaled), y = Z_hat[[h]], standardize = FALSE, alpha=1, nfolds = 10)
		glmnet.mdl <- glmnet(x = t(G.scaled), y = Z_hat[[h]], standardize = FALSE, alpha=1, lambda = glmnet.mdl.cv$lambda.min)
		beta.full <- as.numeric(glmnet.mdl$beta)
		predictors <- rownames(G)[which(beta.full != 0)]
		beta <- as.matrix(c(glmnet.mdl$a0,as.matrix(glmnet.mdl$beta[predictors,])))
		rownames(beta) <- c("Intercept",predictors)
		colnames(beta) <- paste(gene_name,".",cell_type,sep="")
		Z_hat.h.pred <- cbind(numeric(ncol(G))+1,t(G.scaled[predictors,,drop=FALSE])) %*% beta
		if (sd(Z_hat.h.pred) > 0.00001){
			corrs[i,h] <- cor(t(Z_hat[[h]]),Z_hat.h.pred)
			print(sprintf("%d) in-smaple correlation for gene %s, cell type %s: %f", h, gene_name, cell_type, corrs[i,h]))
		}
		write.csv(means[predictors,], file = paste(datadir2,"elastic_net_models/",rownames(genes.info[i,]),".elastic_net.",cell_type,".predictors.means.txt",sep=""), quote = FALSE)
		write.csv(sds[predictors,], file = paste(datadir2,"elastic_net_models/",rownames(genes.info[i,]),".elastic_net.",cell_type,".predictors.sds.txt",sep=""), quote = FALSE)
		write.csv(beta, file = paste(datadir2,"elastic_net_models/",rownames(genes.info[i,]),".elastic_net.",cell_type,".predictors.txt",sep=""), quote = FALSE)
	}
	tic()
}

colnames(corrs) <- colnames(tca.mdl$W)
write.csv(corrs, file = paste(datadir2,"elastic_net_models/elastic_net.in_sample_correlation.txt",sep=""), quote = FALSE)
write.csv(corrs.bulk, file = paste(datadir2,"elastic_net_models/elastic_net.in_sample_correlation.bulk.txt",sep=""), quote = FALSE)


