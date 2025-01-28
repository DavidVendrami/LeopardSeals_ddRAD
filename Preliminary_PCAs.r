library(adegenet)

files = c("Leop_GQDP5_MD80_maxDP_MAF01_HWE001.raw", "Leop_GQDP5_MD90_maxDP_MAF01_HWE001.raw", "Leop_GQDP10_MD80_maxDP129_MAF01_HWE001.raw", "Leop_GQDP10_MD90_maxDP149_MAF01_HWE001.raw")
mains = c("DP = 5; Missingness = 80%", "DP = 5; Missingness = 90%","DP = 10; Missingness = 80%","DP = 10; Missingness = 90%")

pdf("Preliminary_PCAs.pdf",width=9, height=9)
par(mfrow=c(2,2))

for (i in 1:length(files)){
data<-read.PLINK(files[i],n.cores=1)
pca1<-glPca(data,parallel=FALSE, nf = 100)
plot(pca1$scores[,1], pca1$scores[,2], xlab="PC1", ylab="PC2", main=mains[i], pch = 20, cex = 2.5, col=transp("dodgerblue",.4))
}
dev.off()
