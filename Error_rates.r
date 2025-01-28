# Expl filtering erorr rate
data<-read.table('Leop_GQDP5_MD90_maxDP_MAF01_HWE001.raw',h=T)
data<-data[c(30,44),]
data[,c(1:6)]
data<-data[,-c(1:6)]
one<-which(is.na(data[1,]))
two<-which(is.na(data[2,]))
out<-sort(c(one,two))
data<-data[,-out]

err<-c()
for (i in 1:dim(data)[2]){
if (data[1,i] == data[2,i]){
err[i]<-0
} else if (data[1,i] == 1 | data[2,i] == 1) {
err[i]<-0.5
} else {
err[i]<-1
}
}

sum(err)/dim(data)[2]



#Final error rate:
err_rate<-function(data){		
one<-which(is.na(data[1,]))
two<-which(is.na(data[2,]))
out<-sort(c(one,two))
data<-data[,-out]

err<-c()
for (i in 1:dim(data)[2]){
if (data[1,i] == data[2,i]){
err[i]<-0
} else if (data[1,i] == 1 | data[2,i] == 1) {
err[i]<-0.5
} else {
err[i]<-1
}
}
return(sum(err)/dim(data)[2])
}

input<-read.table('DP10_MD90_dups.raw',h=T)
uno<-input[c(1,24),-c(1:6)]
due<-input[c(2,15),-c(1:6)]
tre<-input[c(8,12),-c(1:6)]
qua<-input[c(9,5),-c(1:6)]
cin<-input[c(9,11),-c(1:6)]
sei<-input[c(10,18),-c(1:6)]
set<-input[c(19,23),-c(1:6)]
ott<-input[c(20,21),-c(1:6)]
nov<-input[c(22,16),-c(1:6)]
die<-input[c(25,27),-c(1:6)]
und<-input[c(25,14),-c(1:6)]
dod<-input[c(26,3),-c(1:6)]
trd<-input[c(26,7),-c(1:6)]
qtr<-input[c(27,14),-c(1:6)]
qnc<-input[c(3,7),-c(1:6)]
sdc<-input[c(4,13),-c(1:6)]
dcs<-input[c(5,11),-c(1:6)]
dco<-input[c(6,17),-c(1:6)]

out<-rbind(uno,
           due,
           tre,
           qua,
           cin,
           sei,
           set,
           ott,
           nov,
           die,
           und,
           dod,
           trd,
           qtr,
           qnc,
           sdc,
           dcs,
           dco)
		
errs<-c()
for (i in seq(1,18,by=2)){
errs[i]<-err_rate(out[c(i,i+1),])
}
errsi<-errs[!is.na(errs)]
mean(errsi)
sd(errsi)


