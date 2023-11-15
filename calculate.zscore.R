args = commandArgs(trailingOnly=TRUE)
flank = c(10,20,30,40,50,100,200,300,400,500)
motif = args[1]

sc<-c()
for (f in flank){
results <- read.table(paste0("results-",f,"/",motif,"_permutation_0.0001/results.txt"), quote="", comment.char="")
#View(results)

distribution<-as.integer(results[2:10001,])
m<-mean(distribution)
s<-sd(distribution)
z_score<-(as.integer(results[1,])-m)/s
print(f)
print(z_score)
sc<-c(sc,z_score)
}
sc<-format(round(sc, 2), nsmall = 2)
print(paste(as.character(sc), sep="''", collapse=","))
