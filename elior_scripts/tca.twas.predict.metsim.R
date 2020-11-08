# Predict on the full metsim data (exclude the samples that were used for learning the model)

library("BEDMatrix")
library("data.table")
library("pracma")

plink_file <- "/u/project/eeskin/halperin/ehsafe/FromPaivi/METSIM.FULL/metsim_omni"
models_dir <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/elastic_net_models/"
#models_dir <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/10000000/elastic_net_models/"
datadir1 <- "/u/project/eeskin/halperin/ehsafe/METSIM/"
datadir2 <- "/u/project/eeskin/halperin/erahmani/projects/TCA_expression/metsim/data/"

load("tca.mdl.metsim.Rdata")

INT <- function(x,c = 3/8){
  ranks <- rank(x,na.last = "keep")
  N <- sum(!is.na(x))
  result <- qnorm((ranks-c)/(N - 2*c + 1))
  return(result)
}

G <- BEDMatrix(plink_file)

pcs <- read.table(paste(datadir2,"metsim_omni.pca.eigenvec",sep=""))
pcs.ids_map <- pcs[,1,drop = FALSE]
rownames(pcs.ids_map) <- apply( pcs[,1:2], 1, paste, collapse = "_")
rownames(pcs) <- pcs[,1]
pcs <- pcs[,3:ncol(pcs)]

# Remove the samples for which we don't have all phenotypes and samples which we usedx for fitting the TWAS models; also, consider only smaples for which we have pcs
phenos <- read.csv("/u/project/eeskin/halperin/ehsafe/FromPaivi/METSIM.FULL/F10k_0719.txt", header=T, na.strings=c("","NA"),sep="\t",row.names = 1)
samples <- as.character(setdiff(rownames(phenos)[rowSums(is.na(phenos)) == 0], as.numeric(rownames(tca.mdl$C2))))
samples <- apply( cbind(samples,samples), 1, paste, collapse = "_")
# make sure to consider only sample with genotypes.
samples <- intersect(intersect(samples, rownames(G)),rownames(pcs.ids_map))

samples.phenos <- character(length(samples))
for (i in 1:length(samples)){
	j <- gregexpr(pattern ='_',samples[i])
	samples.phenos[i] <- substr(samples[i], 1, j[[1]][1]-1)
}

covars <- cbind(phenos[samples.phenos,"Age",drop=FALSE], phenos[samples.phenos,"Age",drop=FALSE]**2, as.numeric(phenos[samples.phenos,c("smoke")] == 0), as.numeric(phenos[samples.phenos,c("smoke")] == 1), pcs[as.character(pcs.ids_map[samples,]),1:10])
colnames(covars) <- c("Age","Age_squared","smoke_yes","smoke_no","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")
phenos <- phenos[samples.phenos,3:6]

# recode phenotypes
for (i in 1:4){
	phenos[,i] <- INT(phenos[,i])
}


# TODO: Check for outliers


pheno_names <- c("BMI","S_tottg","matsuda","B_GHbA1C")
G.colnames <- colnames(G)
cell_types <- colnames(tca.mdl$W)
for (i in 1:length(pheno_names)){
	pheno <- pheno_names[i]
	print(pheno)
	candidates <- read.table(file = paste("tcareg.metsim.",pheno_names[i],".marginal.significant_genes.pvals.txt",sep=""), header = TRUE, sep=",", row.names = 1)
	genes <- rownames(candidates)
	y <- INT(phenos[samples.phenos,pheno])
	for (gene in genes){
		if (!file.exists(paste(models_dir,gene,".elastic_net.bulk.predictors.txt",sep=""))) next
		twas.mdl.bulk <- read.table(file = paste(models_dir,gene,".elastic_net.bulk.predictors.txt",sep=""), header = TRUE, sep=",")
		if (nrow(twas.mdl.bulk) == 1) next
		twas.mdl.bulk.means <- read.table(file = paste(models_dir,gene,".elastic_net.bulk.predictors.means.txt",sep=""), header = TRUE, sep=",")[,2]
		twas.mdl.bulk.sds <- read.table(file = paste(models_dir,gene,".elastic_net.bulk.predictors.sds.txt",sep=""), header = TRUE, sep=",")[,2]

		twas.mdl.bulk[,1] <- as.character(twas.mdl.bulk[,1])
		for (a in 2:nrow(twas.mdl.bulk)){
			# Tale care of the case where the suffix of an rs id was change because of a change in the MAF after extracting a subset of the samples in plink
			if (!sum(G.colnames == twas.mdl.bulk[a,1])){
				j <- gregexpr(pattern ='_',twas.mdl.bulk[a,1])
				s <- substr(twas.mdl.bulk[a,1], j[[1]][1]+1, j[[1]][1]+1)
				if (sum(G.colnames == paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_T", sep=""))) twas.mdl.bulk[a,1] <- paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_T", sep="")
				if (sum(G.colnames == paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_A", sep=""))) twas.mdl.bulk[a,1] <- paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_A", sep="")
				if (sum(G.colnames == paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_G", sep=""))) twas.mdl.bulk[a,1] <- paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_G", sep="")
				if (sum(G.colnames == paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_C", sep=""))) twas.mdl.bulk[a,1] <- paste(substr(twas.mdl.bulk[a,1], 1, j[[1]][1]-1), "_C", sep="")			
			}
		}
		Z <- G[samples,as.character(twas.mdl.bulk[2:nrow(twas.mdl.bulk),1])]
		Z.scaled <- (Z - repmat(twas.mdl.bulk.means,length(samples),1)) / repmat(twas.mdl.bulk.sds,length(samples),1)

		x <- data.frame(covars[samples.phenos,], cbind(cbind(numeric(length(samples)),Z.scaled) %*% twas.mdl.bulk[,2]))
		r <- lm(y~.,x)
		pval <- summary(r)[["coefficients"]][ncol(x)+1,4]
		print(sprintf("Phenotype: %s, gene: %s, bulk, pval = %f",pheno,gene,pval))

		for (h in 1:ncol(tca.mdl$W)){
			if (candidates[gene,h]){
				if (!file.exists(paste(models_dir,gene,".elastic_net.",cell_types[h],".predictors.txt",sep=""))) next
				twas.mdl <- read.table(file = paste(models_dir,gene,".elastic_net.",cell_types[h],".predictors.txt",sep=""), header = TRUE, sep=",")
				if (nrow(twas.mdl) == 1) next
				twas.mdl.means <- read.table(file = paste(models_dir,gene,".elastic_net.",cell_types[h],".predictors.means.txt",sep=""), header = TRUE, sep=",")[,2]
				twas.mdl.sds <- read.table(file = paste(models_dir,gene,".elastic_net.",cell_types[h],".predictors.sds.txt",sep=""), header = TRUE, sep=",")[,2]

				twas.mdl[,1] <- as.character(twas.mdl[,1])
				for (a in 2:nrow(twas.mdl)){
					# Tale care of the case where the suffix of an rs id was change because of a change in the MAF after extracting a subset of the samples in plink
					if (!sum(G.colnames == twas.mdl[a,1])){
						j <- gregexpr(pattern ='_',twas.mdl[a,1])
						s <- substr(twas.mdl[a,1], j[[1]][1]+1, j[[1]][1]+1)
						if (sum(G.colnames == paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_T", sep=""))) twas.mdl[a,1] <- paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_T", sep="")
						if (sum(G.colnames == paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_A", sep=""))) twas.mdl[a,1] <- paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_A", sep="")
						if (sum(G.colnames == paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_G", sep=""))) twas.mdl[a,1] <- paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_G", sep="")
						if (sum(G.colnames == paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_C", sep=""))) twas.mdl[a,1] <- paste(substr(twas.mdl[a,1], 1, j[[1]][1]-1), "_C", sep="")		
					}
				}
				Z <- G[samples,as.character(twas.mdl[2:nrow(twas.mdl),1])]
				#Z <- G[samples,as.character(twas.mdl[2:nrow(twas.mdl),1])]
				Z.scaled <- (Z - repmat(twas.mdl.means,length(samples),1)) / repmat(twas.mdl.sds,length(samples),1)
				x <- data.frame(covars[samples.phenos,], cbind(cbind(numeric(length(samples)),Z.scaled) %*% twas.mdl[,2]))
				r <- lm(y~.,x, na.action = "na.omit")
				pval <- summary(r)[["coefficients"]][ncol(x)+1,4]
				print(sprintf("Phenotype: %s, gene: %s, cell type: %s, pval = %f",pheno,gene,cell_types[h],pval))
			}
		}

	}
	
}

