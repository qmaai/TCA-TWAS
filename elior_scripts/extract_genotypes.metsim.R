
## Extract all the genotypes that are in cis to at least one of the genes in the list of candidate genes that were found by tcareg.fit.metsim.sum.R
## a cis snp with respect to a gene is defined as a snp that is within the gene or withing a 'flank' distance from the gene (upstream/downstread)

flank <- 100000
r_squared_th <- 0.8 # prunning of SNPs will be done based on this threshold


library("data.table")

datadir1 <- "/u/project/eeskin/halperin/ehsafe/METSIM/"
datadir2 <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/"

# Extract the chromosome numbers and start, end positions in bp for each of the significant genes
#pheno_names <- c("BMI","S_tottg","P_ffa0","matsuda","fmass","WHR","P_IL1RA","P_Adipon","B_GHbA1C")
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
X <- fread(file = paste(datadir1,"METSIM_variant.dosage",sep=""))


# extract the chr, position information of each snp in X and find the list of SNPs that are within 'flank' distance from at least one of the genes in 'genes'
check_snp <- function(G, str, flank){
	j <- gregexpr(pattern ='_',str)
	chr <- as.numeric(substr(str, 1, j[[1]][1]-1))
	pos <- as.numeric(substr(str, j[[1]][1]+1, j[[1]][2]-1))
	if (sum(G[,1] == chr & (G[,2] < pos & G[,3] > pos | abs(G[,2]-pos) < flank | abs(G[,3]-pos) < flank ))) TRUE else FALSE	
}

snps <- X[,1][2:nrow(X)]
library("tictoc")
tic()
res <- lapply(1:nrow(snps),function(j) check_snp(genes.info, snps[j], flank))
toc()

# extract only the relevant SNPs
a <- logical(nrow(X))
a[1:(length(res)+1)] <- c(TRUE,as.logical(res))
X_prime <- X[a,]
col_names <- X_prime[1,2:ncol(X_prime)]
row_names <- X_prime[2:nrow(X_prime),1]
X_prime <- X_prime[2:nrow(X_prime),2:ncol(X_prime)]
X_prime <- as.matrix(X_prime)
colnames(X_prime) <- col_names
rownames(X_prime) <- as.character(as.matrix(row_names))

# keep only samples that were used in the TCA analysis
phenos <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335.txt",sep="")), row.names=1)
keep <- setdiff(1:nrow(phenos),which(rownames(phenos) == "2812"))
samples <- rownames(phenos)[keep]
X_final <- X_prime[,samples]


## Can run the following to endure low maf and variance snps are excluded; in the metsim data that's already the case so no need to run this

## get the distributions of maf and var (for imputed genotypes) across the selected genotypes
#mafs <- c()
#vars <- c()
#s <- ncol(X_final)
#for (i in 1:nrow(X_final)){
#	if (sum(X_final[i,] == 0 | X_final[i,] == 1 | X_final[i,] == 2) == s){
#		mafs <- append(mafs, min(1-sum(X_final[i,]) / (s*2), sum(X_final[i,]) / (s*2)))
#	}else{
#		vars <- append(vars, var(X_final[i,]))
#	}
#}

## Get the var_th - will be used to filter out imputed genotype with variaiton lower than that
#maf_th <- 0.03
#b <- 1
#while (quantile(mafs, b) > maf_th & b > 0.01) b <- b-0.01
#min_var <- quantile(vars, b)

#variants_keep <- logical(nrow(X_final))
#for (i in 1:nrow(X_final)){
#	if (sum(X_final[i,] == 0 | X_final[i,] == 1 | X_final[i,] == 2) == s){
#		if (min(1-sum(X_final[i,]) / (s*2), sum(X_final[i,]) / (s*2)) >= maf_th) variants_keep[i] <- TRUE
#	}else{
#		if (var(X_final[i,]) >= min_var) variants_keep[i] <- TRUE
#	}
#}
#Z = X_final[variants_keep,]

Z <- X_final
prune_snps <- function(Y, r_squared_th){
	Y_scaled <- t(scale(t(Y), center = TRUE, scale = apply(t(Y), 2, sd)))
	i = 1
	to_remove <- c()
	counter <- 0
	while (i <= nrow(Y_scaled)){
		counter <- counter + 1
		if (max((tcrossprod(Y_scaled[setdiff(1:nrow(Y_scaled),i),],t(Y_scaled[i,]))/ncol(Y_scaled))**2) > r_squared_th){
			Y_scaled <- Y_scaled[setdiff(1:nrow(Y_scaled),i),]
			 append(to_remove, counter)
			next
		}
		i <- i + 1
	}
	return(Y[setdiff(1:nrow(Y_scaled),to_remove),])
}

s <- rownames(Z)
total_num_snps <- 0
# Extract the genotypes to use for each gene
for (i in 1:nrow(genes.info)){
	print(i)
	# get the cis-SNPs for the current gene
	r <- lapply(1:length(s),function(j) check_snp(genes.info[i,], s[j], flank))
	Y <- Z[as.logical(r),]
	Y_pruned <- prune_snps(Y, r_squared_th)
	total_num_snps <- total_num_snps + nrow(Y_pruned)
	# prune SNPs based on correlation; save the genotypes	
	write.csv(Y_pruned, file = paste(datadir1,"cis_genotypes_per_gene/",rownames(genes.info[i,]),".cis_genotypes.txt",sep=""), quote = FALSE)
}


