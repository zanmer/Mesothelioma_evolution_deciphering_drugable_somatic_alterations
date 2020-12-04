library(Palimpsest)
library(deconstructSigs)
library(foreach)
library(reshape2)

load("signature_and_dndsCV.input.Rdata")
## mutational matrix for 96 context
all.propMutsByCat <- palimpsestInput(vcf = all.snv,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
all.numMutsByCat <- palimpsestInput(vcf = all.snv,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = FALSE)
clonal.propMutsByCat <- palimpsestInput(vcf = clonal.snv,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
clonal.numMutsByCat <- palimpsestInput(vcf = clonal.snv,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = FALSE)
subclonal.propMutsByCat <- palimpsestInput(vcf = subclonal.snv,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
subclonal.numMutsByCat <- palimpsestInput(vcf = subclonal.snv,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = FALSE)

## deconstructSigs for all mutations
all.input = data.frame(t(all.numMutsByCat))
all.samples = rownames(all.input)
all.contribution = foreach(i = 1:nrow(all.input),sample=all.samples,.combine=rbind) %do% {
    sig.decons = whichSignatures(
		tumor.ref = all.input,
		signatures.ref = cosmicv3_signatures, 
		sample.id = sample, 
		contexts.needed = TRUE, 
		tri.counts.method = 'default')
    return (sig.decons$weights)
}

## deconstructSigs for clonal mutations
clonal.input = data.frame(t(clonal.numMutsByCat))
clonal.samples = rownames(clonal.input)
clonal.contribution = foreach(i = 1:nrow(clonal.input),sample=clonal.samples,.combine=rbind) %do% {
    sig.decons = whichSignatures(
		tumor.ref = clonal.input,
		signatures.ref = cosmicv3_signatures, 
		sample.id = sample, 
		contexts.needed = TRUE, 
		tri.counts.method = 'default')
    return (sig.decons$weights)
}

## deconstructSigs for subclonal mutations
subclonal.input = data.frame(t(subclonal.numMutsByCat))
subclonal.samples = rownames(subclonal.input)
subclonal.contribution = foreach(i = 1:nrow(subclonal.input),sample=subclonal.samples,.combine=rbind) %do% {
    sig.decons = whichSignatures(
		tumor.ref = subclonal.input,
		signatures.ref = cosmicv3_signatures, 
		sample.id = sample, 
		contexts.needed = TRUE, 
		tri.counts.method = 'default')
    return (sig.decons$weights)
}

sigv3.valid = rownames(subset(sig.names, artefact=="N"))
sigv3.valid = sigv3.valid[sigv3.valid %in% colnames(all.contribution)]
sig_num = length(sigv3.valid)
sig_names = as.character(sig.names[sigv3.valid,1])

all.contribution = all.contribution[,sigv3.valid]
clonal.contribution = clonal.contribution[,sigv3.valid]
subclonal.contribution = subclonal.contribution[,sigv3.valid]


## boxplot
pdf("Signatures.cosmicv3.boxplot.pdf", width=10,height=6)
par(mfrow=(c(4,1)))
par(mar=c(0,8,1,2))
plot(1:sig_num,rep(1,sig_num),bty="n",xaxt="n",yaxt="n",col="#FFFFFF",xlim=c(1,sig_num),xlab="",ylab="",ylim=c(0,1))
text(1:sig_num,0,labels=sig_names,xpd=T,srt=60,adj=0)
par(mar=c(2,8,0,2))
boxplot(all.contribution, names=rep("",sig_num),ylim=c(0,0.8))
mtext(side=2,line=3,"ALL\ncontribution")
par(mar=c(2,8,0,2))
boxplot(clonal.contribution, names=rep("",sig_num),ylim=c(0,0.8))
mtext(side=2,line=3,"Clonal\ncontribution")
par(mar=c(2.5,8,0,2))
boxplot(subclonal.contribution, names=rep("",sig_num),ylim=c(0,0.8))
mtext(side=2,line=3,"Subclonal\ncontribution")
text(1: sig_num,-0.1,labels=sig.names[sigv3.valid,2],xpd=T,cex=1,srt=90,adj=1)
dev.off()

