# This script fits the TCA model to the metseim data.
# At the end, the script saves the estimates parameters of the model.
# Requires R >=3.5; use:
# module load R/3.5.1


## Install the latest local version of TCA
#require("devtools")
#devtools::install_local("/u/project/eeskin/halperin/erahmani/projects/TCA_expression/TCA_1.0.1.tar.gz")

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
W <- t(load_data(paste(datadir2,"metsim.bisque.props.csv",sep="")))
W <- W[colnames(X),]
phenos <- data.frame(fread(file = paste(datadir1,"METSIM_pheno_335.txt",sep="")), row.names=1)
covars <- data.frame(fread(file = paste(datadir1,"METSIM_cov_335.txt",sep="")), row.names=1)

# Remove sample 2812 for which we have NA value in the Matsuda index
keep <- setdiff(1:nrow(phenos),which(rownames(phenos) == "2812"))
X <- X[,keep]
W <- W[keep,]
phenos <- phenos[keep,]
covars <- covars[keep,]
stopifnot(rownames(covars) == rownames(W))
stopifnot(rownames(phenos) == rownames(W))
stopifnot(rownames(phenos) == colnames(X))

#C1 <- cbind(covars[,c("Age")],as.numeric(phenos[,c("smoke")] == 0), as.numeric(phenos[,c("smoke")] == 1))
#rownames(C1) <- rownames(covars)
#colnames(C1) <- c("Age","smoke_yes","smoke_no")	# not that there are also individuals that are either previous smokers or intermittent smokers

# there are 4 batches; create indicators for batches 1-3
batches <- as.numeric(as.factor(covars[,"Batch"]))
C2 <- cbind(as.numeric(phenos[,c("smoke")] == 0), as.numeric(phenos[,c("smoke")] == 1),covars[,c("RIN", "Map_per","MT_per","Bias_3prime")],as.numeric(batches == 1), as.numeric(batches == 2), as.numeric(batches == 3))
rownames(C2) <- rownames(covars)
colnames(C2) <- c("smoke_yes","smoke_no","RIN","Map_per","MT_per","Bias_3prime","Batch_1","Batch_2","Batch_3")

C1 <- covars[,c("Age"),drop=FALSE]

# fit the tca model
#tca.mdl <- tca(X = X, W = W, C1 = C1, C2 = C2, parallel = TRUE, log_file = "tca.mesim.log", debug = TRUE)
tca.mdl <- tca(X = X, W = W, C1 = C1, C2 = C2, parallel = TRUE, log_file = "tca.mesim.log", debug = TRUE)

save('tca.mdl', file = "tca.mdl.metsim.Rdata")

# refit W
tca.mdl <- tca(X = X, W = W, C1 = C1, C2 = C2, refit_W = TRUE, parallel = TRUE, log_file = "tca.refit_W.mesim.log", debug = TRUE)
tca.mdl <- tca(X = X, W = tca.mdl$W, C1 = C1, C2 = C2, parallel = TRUE, log_file = "tca.refitted_W.mesim.log", debug = TRUE)
save('tca.mdl', file = "tca.mdl.refit_W.metsim.Rdata")

