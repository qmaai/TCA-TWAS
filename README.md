# TCA-TWAS
Identification of cell-type-specific genetic regulation of gene expression for transcriptome-wide association studies Genome-wide association studies (GWAS) linearly associate SNPS(DNA) with phenotypes, while Transcriptome-wide assocation studies (TWAS) linearly characterize the assocation with regulation of gene expression by SNPs. One remaining problem is that it remains unclear how SNPs affect phenotypes through different cell types, because traditionally cell-type-specific data are resource-intensive and expensive to acquire. In this TCA-TWAS study, we endevour to deconvolute builk-level gene expression into cell-type-specific gene expressions with SNPs and cell-type weights. We then associate specific gene expressions with phenotypes. Experiments and results are summarized into the poster and the presentation slides.

Two different submodules are contained: data_generation/data_generation.R contains SNPs simulation techniques that satisfies the SNPs heritability requirements on a buld level, with the associate pdf file to illustrate it; pipeline/TCA_gene_bulk_level_heritability contains the jupyternotebook covering the whole pipeline of TCA-TWAS. tca_gene_bulk.sh calls the original package of TCA, the result of which is saved for comparison with TCA-TWAS.
