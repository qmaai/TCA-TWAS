# Requires R >=3.5; use:
# module load R/3.5.1

#require("devtools")
#devtools::install_local("/u/project/eeskin/halperin/erahmani/projects/TCA_expression/TCA_1.0.1.tar.gz")

#pheno_name <- "BMI"
#pheno_name <- "S_tottg"
#pheno_name <- "P_ffa0"
#pheno_name <- "matsuda"
#pheno_name <- "fmass"
#pheno_name <- "WHR"
#pheno_name <- "P_IL1RA"
#pheno_name <- "P_Adipon"
#pheno_name <- "B_GHbA1C"


refit_W <-0
#refit_W <-1

# The number of SVA components to account for
num_sva_comp <- 10 
if (refit_W){
	load("tca.mdl.refit_W.metsim.Rdata")
}else{
	load("tca.mdl.metsim.Rdata")
}

library(TCA)
library("data.table")
datadir1 <- "/u/project/eeskin/halperin/ehsafe/METSIM/"
datadir2 <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/"

load_data <- function(filepath){
	X <- data.frame(fread(file  = filepath), row.names=1)
	col_names <- X[1,]
	X <- X[2:nrow(X),]
	colnames(X) <- col_names
	return (X)	
}

X <- load_data(paste(datadir1,"Kal_v25_TPM_gene_filtered.csv",sep=""))
phenos <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335.txt",sep="")), row.names=1)
HbA1C <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335_HbA1C.txt",sep="")), row.names=1)
phenos <- cbind(phenos,HbA1C)
covars <- data.frame(fread(file = paste(datadir1,"METSIM_cov_335.txt",sep="")), row.names=1)


# Remove sample 2812 for which we have NA value in the Matsuda index
keep <- setdiff(1:nrow(phenos),which(rownames(phenos) == "2812"))
X <- X[,keep]
phenos <- phenos[keep,]
covars <- covars[keep,]

if (pheno_name == "WHR"){
	sva <- t(load_data(paste(datadir1,"sva/METSIM_TPM_SVA_335_WHRAdjBMI.csv",sep="")))	
}else{
	sva <- t(load_data(paste(datadir1,"sva/METSIM_TPM_SVA_335_",pheno_name,".csv",sep="")))
}
sva <- sva[colnames(X),]
sva_comp_names <- character(num_sva_comp)
for (i in 1:num_sva_comp) sva_comp_names[i] <- paste("sva",i,sep="")
C3 <- cbind(tca.mdl$C1,tca.mdl$C2,covars[,c("Age"),drop=FALSE]**2,tca.mdl$W[,1:(ncol(tca.mdl$W)-1)],sva[,1:num_sva_comp])
colnames(C3) <- c(colnames(tca.mdl$C1),colnames(tca.mdl$C2),"Age_squared",colnames(C3)[(ncol(tca.mdl$C1)+ncol(tca.mdl$C2)+2):ncol(C3)])

y <- as.matrix(phenos[,pheno_name])
rownames(y) <- rownames(phenos)
colnames(y) <- pheno_name

# run TCA regression
if (refit_W){
	outfile <- paste("tcareg.refit_W.metsim.",pheno_name,sep="")
}else{
	outfile <- paste("tcareg.metsim.",pheno_name,sep="")
}
res <- tcareg(X = X, tca.mdl = tca.mdl, y = y, C3 = C3, test = "marginal", save_results = TRUE, log_file = paste("tcareg.",pheno_name,".log",sep=""), features_metadata = paste(datadir2,"metsim.features_metadata.txt",sep=""), output = outfile, parallel = TRUE, debug = TRUE)
res <- tcareg(X = X, tca.mdl = tca.mdl, y = y, C3 = C3, test = "single_effect", save_results = TRUE, log_file = paste("tcareg.",pheno_name,".log",sep=""), features_metadata = paste(datadir2,"metsim.features_metadata.txt",sep=""), output = outfile, parallel = TRUE, debug = TRUE)
res <- tcareg(X = X, tca.mdl = tca.mdl, y = y, C3 = C3, test = "joint", save_results = TRUE, log_file = paste("tcareg.",pheno_name,".log",sep=""), features_metadata = paste(datadir2,"metsim.features_metadata.txt",sep=""), output = outfile, parallel = TRUE, debug = TRUE)


# run standard regression analysis
require("pbapply")
reg.pvals <- data.frame(matrix(0,nrow(X),1)+1)
for (i in 1:nrow(X)){
	r <- lm(y~.,data = data.frame(cbind(C3,t(X[i,]))))
	reg.pvals[i,1] <- summary(r)[["coefficients"]][ncol(C3)+2,4]
}
ord <- order(reg.pvals[,1])
reg.pvals <- data.frame(reg.pvals[ord,1])
rownames(reg.pvals) <- rownames(X)[ord]
colnames(reg.pvals) <- c("reg.pvals")

if (refit_W){
	outfile <- paste("regression.refit_W.metsim.",pheno_name,".txt",sep="")
}else{
	outfile <- paste("regression.metsim.",pheno_name,".txt",sep="")
}
write.csv(reg.pvals, file = outfile, quote = FALSE)

