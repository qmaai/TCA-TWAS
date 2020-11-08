## Pool the gene found in all phenotypes of interest and look for cell type specific eqtls for these genes -- only in cell types that were found as candidate cell types in these genes in at least on phenotype.
## Determine the significant cell-type-specific eqtls under FDR
## Then, for each phenotype, consider its candidate (gene,cell-type) pairs, and for each cell-type-specific eqtl that was found for this gene and cell type, test the phenotype for an effect of the cell-type-specific eqtl through changes in cell-type-specific expression.


library(TCA)
library("data.table")
require("stats")
require("rlist")
require("pracma")

datadir1 <- "/u/project/eeskin/halperin/ehsafe/METSIM/"
datadir2 <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/"

num_sva_comp <- 10
pheno_names <- c("BMI","S_tottg","matsuda","B_GHbA1C")
#pheno_names <- c("BMI","S_tottg","P_ffa0","matsuda","fmass","WHR","P_IL1RA","P_Adipon","B_GHbA1C")

assert <- function (expr, error) {
  if (! expr) stop(error, call. = FALSE)
}

load_data <- function(filepath){
	X <- data.frame(fread(file  = filepath), row.names=1)
	col_names <- X[1,]
	X <- X[2:nrow(X),]
	colnames(X) <- col_names
	return (X)	
}
# Returns the (minus) log likelihood of the entire model in a list together with the derivative of tau
# Input:
#  U <- (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
#  const <- -n*log(2*pi)
#  W_squared <- W**2
#  tau
minus_log_likelihood_tau <- function(U,W_squared,sigmas,const,tau){
  m <- ncol(U)
  res <- matrix(0,m,2)
  tmp <- lapply(1:m, function(j)
    {V <- tcrossprod(W_squared,t(sigmas[j,]**2))+tau**2;
    V_squared <- V**2;
    return (c(-0.5*(const-sum(log(V))-sum(U[,j]/V)), -(tau*(sum(U[,j]/V_squared)-sum(1./V))))) } )
  for (j in 1:m){
    res[j,] = tmp[[j]]
  }
  res <- colSums(res)
  return(res[1])
}
require("matrixcalc")
calc_C1_W_interactions <- function(W,C1){
  n <- nrow(W)
  k <- ncol(W)
  p1 <- ncol(C1)
  if (p1){
    return( hadamard.prod(Reshape(Reshape(apply(W, 2, function(v) repmat(v,1,p1)), n*p1*k,1), n,p1*k), repmat(C1, 1, k)) )
  }else{
    return(matrix(0,n,0))
  }
}


X <- load_data(paste(datadir1,"Kal_v25_TPM_gene_filtered.csv",sep=""))

# Remove sample 2812 for which we have NA value in the Matsuda index
phenos <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335.txt",sep="")), row.names=1)
HbA1C <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335_HbA1C.txt",sep="")), row.names=1)
phenos <- cbind(phenos,HbA1C)
keep <- setdiff(1:nrow(phenos),which(rownames(phenos) == "2812"))
X <- X[,keep]
phenos <- phenos[keep,]


load(file = "tca.mdl.metsim.Rdata")

## Test each gene for cell-type-specific eqtls
genes <- c()
for (i in 1:length(pheno_names)){
 	genes <- append(genes, rownames(read.table(file = paste("tcareg.metsim.",pheno_names[i],".marginal.significant_genes.pvals.txt",sep=""), header = TRUE, sep=",", row.names = 1)))
}
genes <- unique(genes)
hits <- matrix(0,length(genes),ncol(tca.mdl$W))
rownames(hits) <- genes
for (i in 1:length(pheno_names)){
	genes.info <- read.table(file = paste(datadir2,"../tcareg.metsim.",pheno_names[i],".marginal.significant_genes.pvals.txt",sep=""), header = TRUE, sep=",", row.names = 1)
 	hits[rownames(genes.info),] <- as.matrix(hits[rownames(genes.info),] + genes.info)
}
colnames(hits) <- colnames(genes.info)

n <- nrow(tca.mdl$W)
const <- -n*log(2*pi)
W_squared <- tca.mdl$W**2
C1_ <- calc_C1_W_interactions(tca.mdl$W,tca.mdl$C1)
C1.map <- matrix(1,1+ncol(tca.mdl$C1),ncol(tca.mdl$W))
results.pvals <- c()
results.ct <- c()
results.genotypes <- c()
results.genes <- c()
for (g in 1:nrow(hits)){
	print(g)
	gene <- rownames(hits)[g]
	G <- read.table(file = paste(datadir1,"cis_genotypes_per_gene/",gene,".cis_genotypes.txt",sep=""), header = TRUE, sep = ",", check.names=FALSE, row.names = 1)
	# calculate the null likelihood - that is, the likelihood of the gene without modeling the effect of a snp
	U <- (tcrossprod(tca.mdl$W,tca.mdl$mus_hat[gene,,drop=FALSE]) + tcrossprod(tca.mdl$C2,tca.mdl$deltas_hat[gene,,drop=FALSE]) + tcrossprod(C1_,tca.mdl$gammas_hat[gene,,drop=FALSE]) - t(X[gene,,drop=FALSE]) )**2
	ll0 <- -minus_log_likelihood_tau(U,W_squared,tca.mdl$sigmas_hat[gene,,drop=FALSE],const,tca.mdl$tau_hat)
	counter <- 1
	num_tests <- sum(hits[gene,]>=1)*nrow(G)
	tmp.pvals <- numeric(num_tests)
	tmp.ct <- character(num_tests)
	tmp.genotypes <- character(num_tests)
	tmp.genes = rep(gene, num_tests)
	for (s in 1:nrow(G)){
		C1.s <- cbind(tca.mdl$C1,t(G[s,]))
		for (h in 1:ncol(tca.mdl$W)){
			if (hits[gene,h]>=1){
				# use C1.map to make sure that the effects of snp s on all cell types except for h are zero (i.e. a marginal test)
				C1.map[nrow(C1.map),h] <- 1
				C1.map[nrow(C1.map),setdiff(1:ncol(tca.mdl$W),h)] <- 0
				tca.mdl.g.1 <- tca(X = X[genes[g],], W = tca.mdl$W, C1 = C1.s, C1.map = C1.map, C2 = tca.mdl$C2, tau = tca.mdl$tau_hat, parallel = FALSE, log_file = NULL, verbose = FALSE)
				C1.s_ <- calc_C1_W_interactions(tca.mdl$W,C1.s)
				U <- (tcrossprod(tca.mdl.g.1$W,tca.mdl.g.1$mus_hat) + tcrossprod(tca.mdl.g.1$C2,tca.mdl.g.1$deltas_hat) + tcrossprod(C1.s_,tca.mdl.g.1$gammas_hat) - t(X[gene,,drop=FALSE]) )**2
				ll1 <- -minus_log_likelihood_tau(U,W_squared,tca.mdl.g.1$sigmas_hat,const,tca.mdl.g.1$tau_hat)
				# null model:
				tca.mdl.g.0 <- tca(X = X[genes[g],], W = tca.mdl$W, C1 = tca.mdl$C1, C2 = tca.mdl$C2, tau = tca.mdl$tau_hat, parallel = FALSE, log_file = NULL, verbose = FALSE)
				U <- (tcrossprod(tca.mdl.g.0$W,tca.mdl.g.0$mus_hat) + tcrossprod(tca.mdl.g.0$C2,tca.mdl.g.0$deltas_hat) + tcrossprod(C1_,tca.mdl.g.0$gammas_hat) - t(X[gene,,drop=FALSE]) )**2
				ll0 <- -minus_log_likelihood_tau(U,W_squared,tca.mdl.g.0$sigmas_hat,const,tca.mdl.g.0$tau_hat)
				#print(pchisq(-2*(ll0-ll1), df = 1, lower.tail=FALSE))
				tmp.pvals[counter] <- pchisq(-2*(ll0-ll1), df = 1, lower.tail=FALSE)
				tmp.ct[counter] <- colnames(hits)[h]
				tmp.genotypes[counter] <- rownames(G)[s]
				counter <- counter + 1
			}
		}
	}
	assert(num_tests == counter-1, "num_tests must equal counter-1")
	results.pvals <- c(results.pvals,tmp.pvals)
	results.ct <- c(results.ct,tmp.ct)
	results.genotypes <- c(results.genotypes,tmp.genotypes)
	results.genes <- c(results.genes,tmp.genes)
}

# calculate qvals, sort, and save into file
df <- data.frame(results.genes,results.genotypes,results.ct,results.pvals,p.adjust(results.pvals, method="BY"))
colnames(df) <- c("gene_id","snp","cell_type","pval","qval")
ord <- order(df$pval)
df <- data.frame(df[ord,])
write.csv(df, file = "cts_eqtl.metsim.txt", quote = FALSE, row.names = FALSE)
# df <- read.table("cts_eqtl.metsim.txt", header = TRUE, sep =",")




## For each of the significant cell-type-specific eqtls, for each of the phenotypes that were found associated with the gene, test the phenotype for a cell-type-specific effect that goes through cell type speicfic expression.
candidates <- df[df$pval < 0.05/nrow(df),]
write.csv(candidates, file = "cts_eqtl.metsim.pval_significant.txt", quote = FALSE, row.names = FALSE)
head(candidates)
candidates.genes <- unique(candidates[,1])

final.pvals <- c()
final.phenos <- c()
final.genes <- c()
final.genotypes <- c()
final.cts <- c()



cell_type_names <- c("Adipocytes","Fibroblasts","Endothelial.cells","Macrophages","T.memory.cells")

for (pheno in 1:length(pheno_names)){
	print(pheno_names[pheno])
	genes.info <- read.table(file = paste("tcareg.metsim.",pheno_names[pheno],".marginal.significant_genes.txt",sep=""), header = TRUE, sep=",", row.names = 1)
 	genes <- rownames(genes.info)
	if(!length(genes)) next

	y <- phenos[colnames(X),pheno_names[pheno],drop=FALSE]
	if (pheno_names[pheno] == "WHR"){
		sva <- t(load_data(paste(datadir1,"sva/METSIM_TPM_SVA_335_WHRAdjBMI.csv",sep="")))	
	}else{
		sva <- t(load_data(paste(datadir1,"sva/METSIM_TPM_SVA_335_",pheno_names[pheno],".csv",sep="")))
	}
	sva <- sva[colnames(X),]
	sva_comp_names <- character(num_sva_comp)
	for (i in 1:num_sva_comp) sva_comp_names[i] <- paste("sva",i,sep="")
	C3 <- cbind(tca.mdl$C1, tca.mdl$C2,tca.mdl$W[,1:(ncol(tca.mdl$W)-1)],sva[,1:num_sva_comp])

	for (g in 1:length(candidates.genes)){
		candidate <- as.character(candidates.genes[g])
		if (sum(genes == candidate) > 0){
			G <- read.table(file = paste(datadir1,"cis_genotypes_per_gene/",candidate,".cis_genotypes.txt",sep=""), header = TRUE, sep = ",", check.names=FALSE, row.names = 1)
			for (h in 1:ncol(tca.mdl$W)){
				if (genes.info[candidate,h]){
					# test each snp that was found associated with the gene 'candidate' in cell type h (if such SNPs exist)
					for (s in 1:nrow(candidates)){
						if (candidates[s,1] == candidate & candidates[s,3] == cell_type_names[h] ){
							print(s)
							# model with cell type specific effect of snp s on cell type h
							C1.s <- cbind(tca.mdl$C1,t(G[as.character(candidates[s,2]),]))
							C1.map[nrow(C1.map),h] <- 1
							C1.map[nrow(C1.map),setdiff(1:ncol(tca.mdl$W),h)] <- 0
							tca.mdl.g.1 <- tca(X = X[genes[g],], W = tca.mdl$W, C1 = C1.s, C1.map = C1.map, C2 = tca.mdl$C2, tau = tca.mdl$tau, parallel = FALSE, log_file = NULL, verbose = FALSE)
							# null model
							tca.mdl.g.0 <- tca(X = X[genes[g],], W = tca.mdl$W, C1 = tca.mdl$C1, C2 = tca.mdl$C2, tau = tca.mdl$tau, parallel = FALSE, log_file = NULL, verbose = FALSE)
							# fit null and alternative for y
							C3_ <- cbind(C3,t(G[as.character(candidates[s,2]),]))
							mdl0 <- tcareg(X = X[genes[g],], tca.mdl = tca.mdl.g.0, y = y, C3 = C3_, test = "custom", alternative_model = c(colnames(tca.mdl$W)[h]), save_results = FALSE, log_file = NULL, features_metadata = NULL, parallel = FALSE, verbose = FALSE)
							mdl1 <- tcareg(X = X[genes[g],], tca.mdl = tca.mdl.g.1, y = y, C3 = C3_, test = "custom", alternative_model = c(colnames(tca.mdl$W)[h]), save_results = FALSE, log_file = NULL, features_metadata = NULL, parallel = FALSE, verbose = FALSE)

							final.pvals <- c(final.pvals,pchisq(-2*(mdl0$alternative_ll-mdl1$alternative_ll), df = 1, lower.tail=FALSE))
							final.cts <- c(final.cts,colnames(tca.mdl$W)[h])
							final.genes <- c(final.genes,candidate)
							final.genotypes <- c(final.genotypes,as.character(candidates[s,2]))
							final.phenos <- c(final.phenos,pheno_names[pheno])

						}
					}
				}
			}
		}
	}
}

df.final <- data.frame(final.phenos,final.genes,final.cts,final.genotypes,final.pvals,p.adjust(final.pvals, method="BY"))
colnames(df.final) <- c("pheno","gene_id","snp","cell_type","pval","qval")
ord <- order(df.final$pval)
df.final <- data.frame(df.final[ord,])

write.csv(df.final, file = "cts_eqtl_pheno.metsim.txt", quote = FALSE, row.names = FALSE)

	