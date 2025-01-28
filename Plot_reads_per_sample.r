library(adegenet)

data<-read.table('Reads_perSample.txt', h=T)
data<-data[-90,] # Remove Control
data <- data[order(data$Total_Sequences),]
reads <- data$Total_Sequences
samples <- gsub("_1", "", data$Sample)

pdf("Reads_per_sample.pdf",width=4, height=8)
plot(reads, 1:length(data$Total_Sequences), axes = F, ylab = 'Samples', xlab = 'Millions of reads', 
     pch = 20, cex = 2, col = transp('black', .3))

box(bty="l")
axis(2, 1:length(samples), labels = samples, las = 2, cex.axis=0.4)
axis(1, seq(0, 7000000, by = 1000000), labels = as.character(seq(0, 7, by = 1)))

abline(v = 1000000, col= 'red', lty = 2)
dev.off()
