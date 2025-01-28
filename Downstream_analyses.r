####### Relatedness
# Bash
mkdir relatedness
cd relatedness
cp ../Fin_DP10_MD90_Missi80.* ./

# A bit of processing in R
R
data<-read.table("Fin_DP10_MD90_Missi80.fam",h=F)
data$V1<-rep('Leop',length(data$V1))
write.table(data,"Fin_DP10_MD90_Missi80.fam",quote=F,col.names=F,row.names=F)
q()
n

# Convert to ped/map and process map file
# bash
/prj/mar-in-gen/bin/plink/plink --bfile Fin_DP10_MD90_Missi80 --recode --aec --out useme

R
data<-read.table('useme.map',h=F)
data$V1<-rep(1,length(data$V1))
write.table(data,"useme.map",quote=F,col.names=F,row.names=F)
q()
n

# Convert to bed
/prj/mar-in-gen/bin/plink/plink --file useme --make-bed --aec --out useme

# Calculate king kinship coefficients and rxy values with gcta
/prj/mar-in-gen/bin/king -b useme.bed --related
/prj/mar-in-gen/bin/gcta_1.93.2beta/gcta64 --bfile useme --make-grm-gz --out relat
gzip -d relat.grm.gz

# Check relatedness patterns from gcta file in r
ids<-read.table("relat.grm.id",h=F)
grm<-read.table("relat.grm",h=F)
dim(grm[(grm$V4 > 0.8 & (grm$V1!=grm$V2)),])
# 0
dim(grm[(grm$V4 > 0.3 & grm < 0.7 & (grm$V1!=grm$V2)),])
# 1
dim(grm[(grm$V4 > 0.1 & grm < 0.3 & (grm$V1!=grm$V2)),])
# 1

# Plot relatedness results
data<-read.table('Data/king.kin',h=T)
data$cols<-data$Kinship
pos<-which(data$cols>=0.1)
rel<-which(data$cols<0.1 & data$cols>0.07)
unrel<-which(data$cols<0.07)
data$cols[pos]<-"#d95f02"
data$cols[rel]<-"#1b9e77"
data$cols[unrel]<-"#7570b3"

# Function to prepare grm
clean_grm <- function(grm, ids){
colnames(ids)<-c('FID','IID')
fid1<-grm$V1
fid2<-grm$V2
iid1<-grm$V1
iid2<-grm$V2

for (i in 1:dim(ids)[1]){
fid1[fid1==i]<-ids$FID[i]
fid2[fid2==i]<-ids$FID[i]
iid1[iid1==i]<-ids$IID[i]
iid2[iid2==i]<-ids$IID[i]
}
grm<-data.frame(FID1=fid1,IID1=iid1,FID2=fid2,IID2=iid2,rxy=grm$V4)
ind<-which(grm$IID1==grm$IID2)
grm<-grm[-ind,]
return(grm)
}

grm<-read.table('Data/relat.grm', h=F)
ids<-read.table('Data/relat.grm.id', h=F)
grm<-clean_grm(grm,ids)

# PO = 146 - 139391: 0.457
# 2nd/3rd = 127 - 139391: 0.21
# more distant: 146 127: 0.0957991

grm$cols<-grm$rxy
pos<-which(grm$cols>=0.4)
rel1<-which(grm$cols<0.4 & grm$cols>0.2)
rel2<-which(grm$cols<0.2 & grm$cols>0.05)
unrel<-which(grm$cols<0.05)
grm$cols[pos]<-"#d95f02"
grm$cols[rel1]<-"#1b9e77"
grm$cols[rel2]<-"#e6ab02"
grm$cols[unrel]<-"#7570b3"

pdf("Relatedness.pdf",width=11,height=6)
par(mfrow=c(1,2))
cols<-c(rep("#7570b3",21),rep("#1f78b4",31),rep("#1b9e77",40),rep("#d95f02",64))
hist(grm$rxy,breaks=100,col=cols,xaxt='n',xlab="",xlim=c(-0.1,0.5), main="",ylim=c(0,500))
arrows(0.095+0.0025,100,0.095+0.0025,25,length=0.1)
arrows(0.0025+0.210,100,0.0025+0.210,25,length=0.1)
arrows(0.0025+0.455,100,0.0025+0.455,25,length=0.1)
put.fig.letter(label="a)", location="topleft",cex=1.5)
abline(h=0)
axis(1,seq(-0.1,0.5,by=0.1),pos=0)
par(xpd=T)
text(0.2,-75,"Genetic relatedness")
par(xpd=F)

plot(data$IBS0,data$Kinship,pch=20,col=transp(data$cols,.6),cex=2,xlab="Proportion of 0 IBS", ylab="KING Kinship coefficient")
abline(h=c(0.354,0.177,0.0884,0.0442),lty=2,col='grey50')
put.fig.letter(label="b)", location="topleft",cex=1.5)
dev.off()

######### 13. Other analyses (PCA, F, MLH,g2)
# PCA
library(adegenet)
library(dplyr)
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.02,0.9),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, adj=0,...)
}


data<-read.PLINK("Data/Fin_DP10_MD90_Missi80.raw",n.cores=1)
lu<-read.table('Data/Look_up.txt',h=T)
data@pop<-factor(rep('LeopS',64))
pca1<-glPca(data,nf=50,parallel=FALSE)
cols<-rep("#7570b3",64)
cols[c(26,32)]<-"#d95f02"
cols[c(18)]<-"#1b9e77"
colors <- data.frame(ID = data@ind.names)
colors <- left_join(colors, lu, by = 'ID')
colors$Sex[colors$Sex == 'Female'] <- "#d95f02"
colors$Sex[colors$Sex == 'Male'] <- "#1b9e77"
colors$Sex[colors$Sex == 'Unknown'] <- "#7570b3"
colors$Residence_status[colors$Residence_status == 'Resident'] <- "#d95f02"
colors$Residence_status[colors$Residence_status == 'Transient'] <- "#1b9e77"
colors$Residence_status[colors$Residence_status == 'Unknown'] <- "#7570b3"

pdf("PCAs.pdf",width=7,height=10.5)
par(mfrow=c(3,2))
plot(pca1$scores[,1],pca1$scores[,2],pch=20,cex=3,col=transp(cols,.6),xlab='PC1 - 2.7%',ylab='PC2 - 1.9%',main='')
#legend(10,10,legend=c("Females","Males","Unknown"),fill=c("#d95f02","#1b9e77","#7570b3"),pch=20,col=c("#d95f02","#1b9e77","#7570b3"))
legend(10,10,legend=c("Parent-offspring","2nd/3rd Degree","Unknown"),pch=20,col=transp(c("#d95f02","#1b9e77","#7570b3"),.6),pt.cex=3)
put.fig.letter(label="a)", location="topleft",cex=1.5)
plot(pca1$scores[,2],pca1$scores[,3],pch=20,cex=3,col=transp(cols,.6),xlab='PC2 - 1.9%',ylab='PC3 - 1.9%',main='')
put.fig.letter(label="b)", location="topleft",cex=1.5)
plot(pca1$scores[,1],pca1$scores[,2],pch=20,cex=3,col=transp(colors$Residence_status,.6),xlab='PC1 - 2.7%',ylab='PC2 - 1.9%',main='')
#legend(10,10,legend=c("Females","Males","Unknown"),fill=c("#d95f02","#1b9e77","#7570b3"),pch=20,col=c("#d95f02","#1b9e77","#7570b3"))
legend(10,10,legend=c("Resident","Transient","Unknown"),pch=20,col=transp(c("#d95f02","#1b9e77","#7570b3"),.6),pt.cex=3)
put.fig.letter(label="c)", location="topleft",cex=1.5)
plot(pca1$scores[,2],pca1$scores[,3],pch=20,cex=3,col=transp(colors$Residence_status,.6),xlab='PC2 - 1.9%',ylab='PC3 - 1.9%',main='')
put.fig.letter(label="d)", location="topleft",cex=1.5)
plot(pca1$scores[,1],pca1$scores[,2],pch=20,cex=3,col=transp(colors$Sex,.6),xlab='PC1 - 2.7%',ylab='PC2 - 1.9%',main='')
#legend(10,10,legend=c("Females","Males","Unknown"),fill=c("#d95f02","#1b9e77","#7570b3"),pch=20,col=c("#d95f02","#1b9e77","#7570b3"))
legend(10,10,legend=c("Female","Male","Unknown"),pch=20,col=transp(c("#d95f02","#1b9e77","#7570b3"),.6),pt.cex=3)
put.fig.letter(label="e)", location="topleft",cex=1.5)
plot(pca1$scores[,2],pca1$scores[,3],pch=20,cex=3,col=transp(colors$Sex,.6),xlab='PC2 - 1.9%',ylab='PC3 - 1.9%',main='')
put.fig.letter(label="f)", location="topleft",cex=1.5)
dev.off()

# Fhat - sMLH
# Calulate Fhat stats
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile Fin_DP10_MD90_Missi80 --aec --ibc --out Leop_inbreeding

# Calculate MLH, sMLH and combine with Fhat
library(inbreedR)
library(dplyr)

ibc<-read.table('Data/Leop_inbreeding.ibc',h=T)
data<-read.table('Data/Fin_DP10_MD90_Missi80.raw',h=T)
geno<-data[,-c(1:6)]
geno[geno==2]<-0
smlh<-sMLH(geno)
ibc$sMLH<-smlh
mlh<-MLH(geno)
ibc$MLH<-mlh

summary(lm(ibc$Fhat1 ~ ibc$sMLH))

# Call:
# lm(formula = ibc$Fhat1 ~ ibc$sMLH)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.116460 -0.027583  0.001209  0.024715  0.086414 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   0.4569     0.1488   3.071  0.00317 **
# ibc$sMLH     -0.4565     0.1487  -3.069  0.00319 **
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.04095 on 62 degrees of freedom
# Multiple R-squared:  0.1319,    Adjusted R-squared:  0.1179 
# F-statistic: 9.418 on 1 and 62 DF,  p-value: 0.003185

summary(lm(ibc$Fhat2 ~ ibc$sMLH))

# Call:
# lm(formula = ibc$Fhat2 ~ ibc$sMLH)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.060791 -0.017859 -0.002798  0.017060  0.051670 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.95140    0.09236   10.30 4.67e-15 ***
# ibc$sMLH    -0.94543    0.09233  -10.24 5.93e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02542 on 62 degrees of freedom
# Multiple R-squared:  0.6284,    Adjusted R-squared:  0.6224 
# F-statistic: 104.8 on 1 and 62 DF,  p-value: 5.929e-15
 
summary(lm(ibc$Fhat3 ~ ibc$sMLH))
 
# Call:
# lm(formula = ibc$Fhat3 ~ ibc$sMLH)
# 
# Residuals:
#       Min        1Q    Median        3Q       Max 
# -0.032963 -0.009412 -0.001232  0.008756  0.036329 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.70752    0.04897   14.45   <2e-16 ***
# ibc$sMLH    -0.70135    0.04896  -14.33   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01348 on 62 degrees of freedom
# Multiple R-squared:  0.768,     Adjusted R-squared:  0.7642 
# F-statistic: 205.2 on 1 and 62 DF,  p-value: < 2.2e-16

data<-data.frame(IID = ibc$IID, Measure=c(rep("MLH",64),rep("sMLH",64),rep("Fhat1",64),rep("Fhat2",64),rep("Fhat3",64)),Value=c(ibc$MLH,ibc$sMLH,ibc$Fhat1,ibc$Fhat2,ibc$Fhat3))
data$Measure<-factor(data$Measure,levels=c("MLH","sMLH","Fhat1","Fhat2","Fhat3"))
pdf('Genetic_diversity.pdf',height=9,width=6)
boxplot(data$Value ~ data$Measure,lwd=1, lty=1,staplelty = 0, boxwex=0.55, outline=F, col=transp(c("#7570b3","#7570b3","#d95f02","#d95f02","#d95f02"),.6),xlab="Measures of genetic diversity", ylab="Value")
stripchart(data$Value ~ data$Measure, vertical = TRUE,method = "jitter", jitter=0.2, add = TRUE, pch = 20, col = transp("black",.2),cex=2)
abline(h=c(0,1),lty=2,col='grey55')
dev.off()

pdf('Genetic_diversity_2.pdf',height=7.5,width=5)
data <- data[data$Measure == 'MLH' | data$Measure == 'Fhat3',]
data$Measure <- as.character(data$Measure)
data$Measure <- factor(data$Measure, levels= c('MLH','Fhat3'))
boxplot(data$Value ~ data$Measure,lwd=1, lty=1,staplelty = 0, boxwex=0.55, outline=F, col=transp(c("#7570b3","#d95f02"),.6),xlab="Measures of genetic diversity", ylab="Value")
stripchart(data$Value ~ data$Measure, vertical = TRUE,method = "jitter", jitter=0.2, add = TRUE, pch = 20, col = transp("black",.2),cex=2)
abline(h=c(0,1),lty=2,col='grey55')
dev.off()

pdf('Genetic_diversity_3.pdf',height=7.5,width=5)
data <- data[data$Measure == 'sMLH' | data$Measure == 'Fhat3',]
data$Measure <- as.character(data$Measure)
data$Measure <- factor(data$Measure, levels= c('sMLH','Fhat3'))
boxplot(data$Value ~ data$Measure,lwd=1, lty=1,staplelty = 0, boxwex=0.55, outline=F, col=transp(c("#7570b3","#d95f02"),.6),xlab="Measures of genetic diversity", ylab="Value")
stripchart(data$Value ~ data$Measure, vertical = TRUE,method = "jitter", jitter=0.2, add = TRUE, pch = 20, col = transp("black",.2),cex=2)
abline(h=c(0,1),lty=2,col='grey55')
dev.off()

lku <- read.table('Data/Mappe1.txt', h = T) |>
	select(ID, Residence_status) |>
	rename(IID = ID)

data_st <- left_join(data, lku)
data_st <- data_st[data_st$Measure == 'MLH' | data_st$Measure == 'Fhat3',]
data_st <- data_st[data_st$Residence_status != 'Unknown',]
data_st$Group <- paste(data_st$Measure, data_st$Residence_status, sep = '_')
data_st$Group <- factor(data_st$Group, levels = c('MLH_Resident', 'MLH_Transient', 'Fhat3_Resident', 'Fhat3_Transient'))
data_st$Measure <- factor(as.character(data_st$Measure), levels = c('MLH', 'Fhat3'))

comp <- ibc |>
	select(IID, Fhat3, sMLH, MLH) |>
	left_join(lku) |>
	filter(Residence_status != 'Unknown')

t.test(comp$sMLH[comp$Residence_status == 'Transient'], comp$sMLH[comp$Residence_status == 'Resident'])
#         Welch Two Sample t-test
# 
# data:  comp$sMLH[comp$Residence_status == "Transient"] and comp$sMLH[comp$Residence_status == "Resident"]
# t = -1.6261, df = 51.531, p-value = 0.11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.031774671  0.003331685
# sample estimates:
# mean of x mean of y 
# 0.9934478 1.0076693

t.test(comp$MLH[comp$Residence_status == 'Transient'], comp$MLH[comp$Residence_status == 'Resident'])
#         Welch Two Sample t-test
# 
# data:  comp$MLH[comp$Residence_status == "Transient"] and comp$MLH[comp$Residence_status == "Resident"]
# t = -1.5792, df = 51.219, p-value = 0.1205
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.0058593481  0.0006995739
# sample estimates:
# mean of x mean of y 
# 0.1818241 0.1844040

t.test(comp$Fhat3[comp$Residence_status == 'Transient'], comp$Fhat3[comp$Residence_status == 'Resident'])
#         Welch Two Sample t-test
# 
# data:  comp$Fhat3[comp$Residence_status == "Transient"] and comp$Fhat3[comp$Residence_status == "Resident"]
# t = 0.74928, df = 44.379, p-value = 0.4576
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.008908011  0.019455791
# sample estimates:
#   mean of x   mean of y 
# 0.007843400 0.002569509

# g2:
data<-read.table('Data/Fin_DP10_MD90_Missi80.raw',h=T)
trans <- lku[which(lku$Residence_status == 'Transient'), 1]
res <- lku[which(lku$Residence_status == 'Resident'), 1]
trans_geno <- data[which(data$IID %in% trans), -c(1:6)]
trans_geno[trans_geno == 2] <- 0
res_geno <- data[which(data$IID %in% res), -c(1:6)]
res_geno[res_geno == 2] <- 0

g2_trans <- g2_snps(trans_geno, nperm = 100, nboot = 1000, CI = 0.95, parallel = FALSE, ncores = NULL)
g2_res <- g2_snps(res_geno, nperm = 100, nboot = 1000, CI = 0.95, parallel = FALSE, ncores = NULL)

g2_res
# Calculation of identity disequilibrium with g2 for SNP data
# -----------------------------------------------------------
# 
# Data: 33 observations at 10555 markers
# Function call = g2_snps(genotypes = res_geno, nperm = 100, nboot = 1000, CI = 0.95,     parallel = FALSE, ncores = NULL)
# 
# g2 = 0.0006320129, se = 0.0007065693
# 
# confidence interval 
#          2.5%         97.5% 
# -0.0002693657  0.0021637742 
# 
# p (g2 > 0) = 0.01 (based on 99 permutations)> 

g2_trans
# Calculation of identity disequilibrium with g2 for SNP data
# -----------------------------------------------------------
# 
# Data: 29 observations at 10555 markers
# Function call = g2_snps(genotypes = trans_geno, nperm = 100, nboot = 1000, CI = 0.95,     parallel = FALSE, ncores = NULL)
# 
# g2 = 0.001806035, se = 0.00165927
# 
# confidence interval 
#          2.5%         97.5% 
# -9.423416e-05  5.280208e-03 
# 
# p (g2 > 0) = 0.01 (based on 99 permutations)

tiff('Supplementary_Figure_B.tiff', width = 10, height = 6, unit = 'in', res = 600)
par(mfrow = c(1, 2))
plot(g2_res, main = 'Residents', col = transp("#d95f02", .6))
plot(g2_trans, main = 'Transients', col = transp("#1b9e77", .6))
dev.off()

g2r <- 0.0006320129
g2t <- 0.001806035
g2rl <- -0.0002693657
g2rh <- 0.0021637742
g2tl <- -9.423416e-05
g2th <- 5.280208e-03

tiff('Supplementary_Figure_B.tiff',height=9,width=7, unit= 'in', res = 600)
layout(matrix(c(1,2,3,3), 2, 2, byrow = T)) 
hist(g2_res$g2_boot, col = transp("#d95f02", .6), ylim = c(0,500), main = '', xlab = 'g2')
points(g2r,450,col = transp("#d95f02", .6), pch = 16, cex= 2)
segments(g2rl,450,g2rh,450, col = transp("#d95f02", .6), lwd = 2)
hist(g2_trans$g2_boot, col = transp("#1b9e77", .6), ylim = c(0,500), main = '', xlab = 'g2')
points(g2t,450,col = transp("#1b9e77", .6), pch = 16, cex= 2)
segments(g2tl,450,g2th,450, col = transp("#1b9e77", .6), lwd = 2)
boxplot(data_st$Value ~ data_st$Group,lwd=1, lty=1,staplelty = 0, boxwex=0.55, outline=F, 
	col=transp(c("#d95f02","#1b9e77","#d95f02","#1b9e77"),.6),xlab="", ylab="Value", at=c(1,2,4,5),
	names = c('Resident','Transient','Resident','Transient'))
stripchart(data_st$Value ~ data_st$Group, vertical = TRUE,method = "jitter", jitter=0.2, add = TRUE, pch = 20, col = transp(c("black"),.2),cex=2, at=c(1,2,4,5))
dev.off()
