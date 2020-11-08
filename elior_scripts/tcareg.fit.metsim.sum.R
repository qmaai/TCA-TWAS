## Sum the results of tcareg.fit.metsim.R
## For each phenotype create and save a matrix with the associated genes found and 0/1 value for each cell type, indicating whther an associaiton with the gene though the corresponding cell type was reported.


refit_W = 0
library("data.table")
library("stats")
#pheno_names <- c("BMI","S_tottg","P_ffa0","matsuda","fmass","WHR","P_IL1RA","P_Adipon")
pheno_names <- c("BMI","S_tottg","matsuda","B_GHbA1C")
cell_types <- c("Adipocytes","Fibroblasts","Endothelial cells","Macrophages","T memory cells")
qval_th <- 0.05
pval_th <- 0.05
test_names <- c("regression","joint","single_effect","marginal_ct1","marginal_ct2","marginal_ct3","marginal_ct4","marginal_ct5")
results.qvals <- matrix(0,length(pheno_names),length(test_names))
refit_W_str <- ""
if (refit_W) refit_W_str <- "refit_W."

# qvals
for (i in 1:length(pheno_names)){
	print(i)
	pheno_name <- pheno_names[i]
	res <- data.frame(fread(file = paste("regression.", refit_W_str, "metsim.",pheno_name,".txt",sep="")), row.names=1)
	results.qvals[i,1] <- sum(p.adjust(res$reg.pval, method = "BY") < qval_th)
	res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".joint.txt",sep="")), row.names=1)
	results.qvals[i,2] <- sum(p.adjust(res$pval, method = "BY") < qval_th)
	res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".single_effect.txt",sep="")), row.names=1)
	results.qvals[i,3] <- sum(p.adjust(res$pval, method = "BY") < qval_th)
	genes <- c()
	for (h in 1:length(cell_types)){
		res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".marginal.",cell_types[h],".txt",sep="")), row.names=1)
		qvals <- p.adjust(res$pval, method = "BY")
		results.qvals[i,h+3] <- sum(qvals < qval_th)
		genes <- append(genes, rownames(res[qvals < qval_th,]))
	}
	genes <- unique(genes)
	Y <- matrix(0,length(genes),length(cell_types))
	rownames(Y) <- genes
	colnames(Y) <- cell_types
	for (h in 1:length(cell_types)){
		res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".marginal.",cell_types[h],".txt",sep="")), row.names=1)
		qvals <- p.adjust(res$pval, method = "BY")
		Y[rownames(res[qvals < qval_th,]),h] <- 1
	}
	write.table(Y, file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".marginal.significant_genes.qvals.txt",sep=""), quote = FALSE, row.names = TRUE, col.names = NA, sep=",")
}
rownames(results.qvals) <- pheno_names
colnames(results.qvals) <- test_names


# pvals
results.pvals <- matrix(0,length(pheno_names),length(test_names))
for (i in 1:length(pheno_names)){
	print(i)
	pheno_name <- pheno_names[i]
	res <- data.frame(fread(file = paste("regression.", refit_W_str, "metsim.",pheno_name,".txt",sep="")), row.names=1)
	results.pvals[i,1] <- sum(res$reg.pval < pval_th/length(res$reg.pvals))
	res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".joint.txt",sep="")), row.names=1)
	results.pvals[i,2] <- sum(res$pval < pval_th/length(res$pval))
	res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".single_effect.txt",sep="")), row.names=1)
	results.pvals[i,3] <- sum(res$pval < pval_th/length(res$pval))
	genes <- c()
	for (h in 1:length(cell_types)){
		res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".marginal.",cell_types[h],".txt",sep="")), row.names=1)
		results.pvals[i,h+3] <- sum(res$pval < pval_th/length(res$pval))
		genes <- append(genes, rownames(res[res$pval < pval_th/length(res$pval),]))
	}
	genes <- unique(genes)
	Y <- matrix(0,length(genes),length(cell_types))
	rownames(Y) <- genes
	colnames(Y) <- cell_types
	for (h in 1:length(cell_types)){
		res <- data.frame(fread(file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".marginal.",cell_types[h],".txt",sep="")), row.names=1)
		Y[rownames(res[res$pval < pval_th/length(res$pval),]),h] <- 1
	}
	write.table(Y, file = paste("tcareg.",refit_W_str,"metsim.",pheno_name,".marginal.significant_genes.pvals.txt",sep=""), quote = FALSE, row.names = TRUE, col.names = NA, sep=",")
}
rownames(results.pvals) <- pheno_names
colnames(results.pvals) <- test_names


