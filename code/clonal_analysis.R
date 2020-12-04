library(foreach)

# Function to calculate cancer cell fraction
# Adapted from Mcgranahan et al, Science
absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number)
{
  f.function <- function (c,purity,local.copy.number)
  {
    return(min(c((purity*c) / (2*(1-purity) + purity*local.copy.number),1)))
  }
  x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
  if(is.na(min(x))){
    d = data.frame(
      ccf.lower.ci = NA, ccf.est = NA, ccf.upper.ci = NA,
      prob.subclonal=NA,prob.clonal=NA, 
      is.clonal=NA,is.clonal.byprob=NA
    )
    return (d)
  }
  if(min(x)==0)
  {
    x[length(x)] <- 1
  }
  names(x)       <- seq(0.01,1,length.out=100)
  sub.cint <- function(x, prob = 0.95,n.alt,depth) {
    xnorm   <- x/sum(x)
    xsort   <- sort(xnorm, decreasing = TRUE)
    xcumLik <- cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh <- xsort[n]
    cint  <- x[xnorm >= LikThresh]
    all   <- as.numeric(names(x))
    cellu <- as.numeric(names(cint))
    l.t   <- cellu[1]
    r.t   <- cellu[length(cellu)]
    m     <- cellu[which.max(cint)]
    
    prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative=‘less’)$p.val
    prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative=‘greater’)$p.val
    
    data.frame(
      ccf.lower.ci = l.t, ccf.est = m, ccf.upper.ci = r.t,
      prob.subclonal=prob.subclonal,prob.clonal=prob.clonal, 
      # Define clonality strictly as a confidence interval overlapping 1
      is.clonal=(l.t <= 1 & r.t >= 1),
      # Define clonality differently as prob.clonal > prob.subclonal (only for comparison purposes)
      is.clonal.byprob=(prob.clonal > prob.subclonal)
    )
  }
  return(sub.cint(x,n.alt=n.alt,depth=depth))
}

load("clonal_analysis.input.Rdata")
n = nrow(purity)

clonal_results = foreach (i=1:n, pat=purity$patient, reg=purity$region, .combine=rbind) %do% {
	cat(pat, reg, "...\n")
	## snv information
	snv=subset(somatic.snv,Tumor_Sample_Barcode==reg & Matched_Norm_Sample_Barcode==pat)
	rownames(snv) = paste(snv$Chromosome, snv$Start_Position, sep="_")
	snv. = data.frame(
		chrom = snv$Chromosome,
		position = snv$Start_Position,
		readcount = snv$t_ref_count + snv$t_alt_count,
		alt_count = snv$t_alt_count,
		vaf = snv$t_alt_count/(snv$t_ref_count + snv$t_alt_count)
	)
	if(nrow(snv.)==0){
		return (NULL)
	}
	if(!grepl("chr",snv$Chromosome[1])){
		snv.$chrom = paste("chr", snv$Chromosome,sep="")
	}
	## cnv information
	cnv = subset(ascat.cnv, patient==pat & region==reg)
	cnv$nTotal = cnv$nMajor+cnv$nMinor

	indices = sapply(1:nrow(snv.), function(x) 
		which(cnv$chr==snv.$chrom[x] & cnv$startpos<=snv.$position[x] & cnv$endpos>=snv.$position[x])[1])
	snv.$cnTotal = cnv[indices, "nTotal"]
	snv.$cnMajor = cnv[indices, "nMajor"]
	snv.$cnMinor = cnv[indices, "nMinor"]

	rho = purity$purity[i] ## tumor purity
	snv.$nmut = snv.$vaf * (rho * snv.$cnTotal + (1-rho) * 2) / rho ## multiplicity
	## calculate CCF
	ccf = lapply(1:nrow(snv.), function (x){
		absolute.cancer.cell.fraction(
			n.alt = snv.$alt_count[x],
			depth = snv.$readcount[x],
			purity = rho,
			local.copy.number=snv.$cnTotal[x]
		)
	})
	ccf <- data.table::rbindlist(ccf)
	return(cbind(snv[,1:(ncol(snv)-4)], snv.[,3:ncol(snv.)],ccf))

}
## final results
head(clonal_results)

## you can formated into tsv and run PyClone for clustering analsysis of multiple regions.
## use patient U_MIST1 for example

## adjust ref_counts and alt_counts according to CCF estimates
clonal_results = subset(clonal_results, !is.na(is.clonal))
clonal_results$adjust_alt_counts = round(clonal_results$readcount * clonal_results$ccf.est/2)
clonal_results$adjust_ref_counts = clonal_results$readcount - clonal_results$adjust_alt_counts

R1.indices = which(clonal_results$Matched_Norm_Sample_Barcode=="U_MIST1" & clonal_results$Tumor_Sample_Barcode=="R1")
R2.indices = which(clonal_results$Matched_Norm_Sample_Barcode=="U_MIST1" & clonal_results$Tumor_Sample_Barcode=="R2")
R3.indices = which(clonal_results$Matched_Norm_Sample_Barcode=="U_MIST1" & clonal_results$Tumor_Sample_Barcode=="R3")
R4.indices = which(clonal_results$Matched_Norm_Sample_Barcode=="U_MIST1" & clonal_results$Tumor_Sample_Barcode=="R4")
R5.indices = which(clonal_results$Matched_Norm_Sample_Barcode=="U_MIST1" & clonal_results$Tumor_Sample_Barcode=="R5")
R1.df = data.frame(
	mutation_id = paste(clonal_results$Chromosome[R1.indices],clonal_results$Start_Position[R1.indices],sep="_"),
	ref_counts = clonal_results$adjust_ref_counts[R1.indices],
	var_counts = clonal_results$adjust_alt_counts[R1.indices],
	normal_cn = 2,
	minor_cn = 0,
	major_cn = 2,
	variant_case = "R1",
	variant_freq = clonal_results$ccf.est[R1.indices]/2,
	ref_allele = clonal_results$Tumor_Seq_Allele1[R1.indices],
	alt_allele = clonal_results$Tumor_Seq_Allele2[R1.indices],
	func = clonal_results$Variant_Classification[R1.indices],
	gene = clonal_results$Hugo_Symbol[R1.indices]
)

R2.df = data.frame(
	mutation_id = paste(clonal_results$Chromosome[R2.indices],clonal_results$Start_Position[R2.indices],sep="_"),
	ref_counts = clonal_results$adjust_ref_counts[R2.indices],
	var_counts = clonal_results$adjust_alt_counts[R2.indices],
	normal_cn = 2,
	minor_cn = 0,
	major_cn = 2,
	variant_case = "R1",
	variant_freq = clonal_results$ccf.est[R2.indices]/2,
	ref_allele = clonal_results$Tumor_Seq_Allele1[R2.indices],
	alt_allele = clonal_results$Tumor_Seq_Allele2[R2.indices],
	func = clonal_results$Variant_Classification[R2.indices],
	gene = clonal_results$Hugo_Symbol[R2.indices]
)

R3.df = data.frame(
	mutation_id = paste(clonal_results$Chromosome[R3.indices],clonal_results$Start_Position[R3.indices],sep="_"),
	ref_counts = clonal_results$adjust_ref_counts[R3.indices],
	var_counts = clonal_results$adjust_alt_counts[R3.indices],
	normal_cn = 2,
	minor_cn = 0,
	major_cn = 2,
	variant_case = "R1",
	variant_freq = clonal_results$ccf.est[R3.indices]/2,
	ref_allele = clonal_results$Tumor_Seq_Allele1[R3.indices],
	alt_allele = clonal_results$Tumor_Seq_Allele2[R3.indices],
	func = clonal_results$Variant_Classification[R3.indices],
	gene = clonal_results$Hugo_Symbol[R3.indices]
)

R4.df = data.frame(
	mutation_id = paste(clonal_results$Chromosome[R4.indices],clonal_results$Start_Position[R4.indices],sep="_"),
	ref_counts = clonal_results$adjust_ref_counts[R4.indices],
	var_counts = clonal_results$adjust_alt_counts[R4.indices],
	normal_cn = 2,
	minor_cn = 0,
	major_cn = 2,
	variant_case = "R1",
	variant_freq = clonal_results$ccf.est[R4.indices]/2,
	ref_allele = clonal_results$Tumor_Seq_Allele1[R4.indices],
	alt_allele = clonal_results$Tumor_Seq_Allele2[R4.indices],
	func = clonal_results$Variant_Classification[R4.indices],
	gene = clonal_results$Hugo_Symbol[R4.indices]
)

R5.df = data.frame(
	mutation_id = paste(clonal_results$Chromosome[R5.indices],clonal_results$Start_Position[R5.indices],sep="_"),
	ref_counts = clonal_results$adjust_ref_counts[R5.indices],
	var_counts = clonal_results$adjust_alt_counts[R5.indices],
	normal_cn = 2,
	minor_cn = 0,
	major_cn = 2,
	variant_case = "R1",
	variant_freq = clonal_results$ccf.est[R5.indices]/2,
	ref_allele = clonal_results$Tumor_Seq_Allele1[R5.indices],
	alt_allele = clonal_results$Tumor_Seq_Allele2[R5.indices],
	func = clonal_results$Variant_Classification[R5.indices],
	gene = clonal_results$Hugo_Symbol[R5.indices]
)

write.table(R1.df, file="U_MIST1.R1.snv.tsv", quote=F, sep="\t", row.names=F)
write.table(R1.df, file="U_MIST1.R2.snv.tsv", quote=F, sep="\t", row.names=F)
write.table(R1.df, file="U_MIST1.R3.snv.tsv", quote=F, sep="\t", row.names=F)
write.table(R1.df, file="U_MIST1.R4.snv.tsv", quote=F, sep="\t", row.names=F)
write.table(R1.df, file="U_MIST1.R5.snv.tsv", quote=F, sep="\t", row.names=F)
