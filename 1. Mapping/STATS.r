#===========Summary stats of mapping (R).
setwd("STATS_mem")

#-----------Check coverage and insert size (R).
files<-grep("\\.stats$",dir(),value=T)
dev.new(height=5,width=10)
par(mfrow=c(1,2))
for(i in 1:length(files))
{
	a<-scan(files[i],what="",quiet=T,sep="\n")
	
	write(grep("^COV",a,value=T),"temp")
	b<-read.table("temp",sep="\t",stringsAsFactors=F)
	b<-b[,3:4]
	plot(b[,1],b[,2],xlim=c(0,100),xlab="Depth",ylab="Frequency")
	
	write(grep("^IS",a,value=T),"temp")
	b<-read.table("temp",sep="\t",stringsAsFactors=F)
	b<-b[,2:3]
	plot(b[,1],b[,2],xlab="Insert size",ylab="Frequency")
	
	savePlot(sub("\\.stats$",".covis.tiff",files[i]),"png")
}
file.remove("temp")
dev.off()

files<-grep("\\.stats$",dir(),value=T)
cov<-numeric(0)
ins<-numeric(0)
for(i in 1:length(files))
{
	a<-scan(files[i],what="",quiet=T,sep="\n")
	
	write(grep("^COV",a,value=T),"temp")
	b<-read.table("temp",sep="\t",stringsAsFactors=F)
	b<-b[1:100,3:4]
	cov<-c(cov,sum(as.numeric(b[,1])*as.numeric(b[,2]))/sum(b[,2]))
	
	write(grep("^IS",a,value=T),"temp")
	b<-read.table("temp",sep="\t",stringsAsFactors=F)
	b<-b[,2:3]
	ins<-c(ins,sum(as.numeric(b[,1])*as.numeric(b[,2]))/sum(b[,2]))
}
file.remove("temp")
stats<-cbind(sub(".stats$","",files),cov,ins)
colnames(stats)<-c("Indi","COV","IS")
write.table(stats,"zstats",col.names=T,row.names=F,quote=F,sep="\t")

#-----------Summary of flagstat files (R).
files <- grep("\\.flagstat$", dir(), value = TRUE)
n <- length(files)
flagstat <- data.frame(Indi = character(n), Reads = integer(n), Mapped = numeric(n), Paired = numeric(n))
flagstat$Indi <- sub(".flagstat$", "", files)

for (i in 1:n) {
    a <- readLines(files[i])
    
    flagstat[i, "Reads"] <- as.integer(strsplit(a[1], " ")[[1]][1])
    
    #extract the mapped reads from the 5th line
    flagstat[i, "Mapped"] <- as.numeric(sub(".*mapped \\(([0-9.]+)%.*", "\\1", a[5])) / 100
    
    #extract the percentage of properly paired reads from the 9th line
    flagstat[i, "Paired"] <- as.numeric(sub(".* properly paired \\(([0-9.]+)%.*", "\\1", a[9])) / 100
}

write.table(flagstat, "zflagstat", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
