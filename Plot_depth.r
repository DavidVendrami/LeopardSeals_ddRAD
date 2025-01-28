files <- c("Leop_Seal_Biall_GQDP10_MD80.ldepth.mean", "Leop_Seal_Biall_GQDP10_MD90.ldepth.mean")
miss <- c("Missingness = 80%", "Missingness = 90%")

pdf("SNPs_Coverage.pdf",width=11, height=5.5)
par(mfrow=c(1,2))
for (i in 1:length(files)){
data <- read.table(files[i], h = T)
hist(data$MEAN_DEPTH, col = 'dodgerblue', xlab = "Mean SNP depth", main = miss[i], breaks = 100, axes = F, xaxt = 'n', xlim = c(0, 300), ylim = c(0, 1000))
axis(1, seq(0, 300, by = 50), labels = as.character(seq(0, 300, by = 50)), pos = 0)
axis(2, seq(0, 1000, by = 200), labels = as.character(seq(0, 1000, by = 200)), pos = 0)
abline(v = 2 * mean(data$MEAN_DEPTH), lty = 2, col = 'red')
}
dev.off()
