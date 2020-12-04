library("dndscv")
load("signature_and_dndsCV.input.Rdata")

titles = c("Sample", "CHROM", "POS", "REF", "ALT")

all.snv = all.snv[, titles]
clonal.snv = clonal.snv[, titles]
subclonal.snv = subclonal.snv[, titles]
colnames(all.snv) = c("sampleID", "chr", "pos", "ref", "mut")
colnames(clonal.snv) = c("sampleID", "chr", "pos", "ref", "mut")
colnames(subclonal.snv) = c("sampleID", "chr", "pos", "ref", "mut")

## remove duplicate snv
all.snv$uid = paste(all.snv$chr, all.snv$pos, sep=":")
clonal.snv$uid = paste(clonal.snv$chr, clonal.snv$pos, sep=":")
subclonal.snv$uid = paste(subclonal.snv$chr, subclonal.snv$pos, sep=":")
all.snv = all.snv[-which(duplicated(all.snv$uid)),]
clonal.snv = clonal.snv[-which(duplicated(clonal.snv$uid)),]
subclonal.snv = subclonal.snv[-which(duplicated(subclonal.snv$uid)),]
all.snv$uid = NULL
clonal.snv$uid = NULL
subclonal.snv$uid = NULL

## run dndsCV program
all.dnds_res = dndscv(all.snv, outmats=T)
clonal.dnds_res = dndscv(clonal.snv, outmats=T)
subclonal.dnds_res = dndscv(subclonal.snv, outmats=T)

all.gene_ci = geneci(all.dnds_res)
clonal.gene_ci = geneci(clonal.dnds_res)
subclonal.gene_ci = geneci(subclonal.dnds_res)

## Find conserved/driver genes
all.sig_genes = subset(all.dnds_res$sel_cv, qglobal_cv<0.1)
clonal.sig_genes = subset(clonal.dnds_res$sel_cv, qglobal_cv<0.1)
subclonal.sig_genes = subset(subclonal.dnds_res$sel_cv, qglobal_cv<0.1)

## show top6 driver genes
head(all.sig_genes)

head(clonal.sig_genes)

head(subclonal.sig_genes)

## global dnds in the cohort
head(all.dnds_res$globaldnds)

head(clonal.dnds_res$globaldnds)

head(subclonal.dnds_res$globaldnds)

