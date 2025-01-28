# Check relatedness patterns

# Start by copying the dataset from a given filtering strategy to a temporary folder

mkdir temp
cd temp
cp ../Leop_GQDP5_MD90_maxDP_MAF01_HWE001.* ./

# A bit of processing in R
R
data<-read.table("Leop_GQDP5_MD90_maxDP_MAF01_HWE001.fam",h=F)
data$V1<-rep('Leop',length(data$V1))
write.table(data,"Leop_GQDP5_MD90_maxDP_MAF01_HWE001.fam",quote=F,col.names=F,row.names=F)
q()
n

# Convert to ped/map and process map file
/prj/mar-in-gen/bin/plink/plink --bfile Leop_GQDP5_MD90_maxDP_MAF01_HWE001 --recode --aec --out useme

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
View(grm[(grm$V4 > 0.8 & (grm$V1!=grm$V2)),])
dim(grm[(grm$V4 > 0.8 & (grm$V1!=grm$V2)),])
dim(grm[(grm$V4 > 0.3 & grm < 0.7 & (grm$V1!=grm$V2)),])

rm *
