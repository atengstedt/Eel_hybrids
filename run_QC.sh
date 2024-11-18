## Quality control of sequencing data

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#===========QC
mkdir QC 

from="/faststorage/project/Eels/eel_combined_Aja/raw_data"
to="/faststorage/project/Eels/eel_combined_Aja/QC"

mkdir $to

for file in ${from}/*.fq.gz
do
sbatch -A Eels -t 12:00:00 --wrap "fastqc -o $to $file"
done

#-----------Info from FastQC result (R)
setwd("/faststorage/project/Eels/eel_combined_Aja/QC")

files<-grep("\\.html$",dir(),value=T)

n<-length(files)

l.read<-integer(n)
n.read<-l.read
encoding<-character(n)

for(i in 1:n)
{
	html<-scan(files[i],what="",quiet=T,sep="\n")
	html<-grep("<td>Sequence length</td>",html,value=T)
	a<-unlist(strsplit(html,"<[[:alpha:]/]{1,}>"))
	l.read[i]<-a[which(a=="Sequence length")+2]
	n.read[i]<-a[which(a=="Total Sequences")+2]
	encoding[i]<-a[which(a=="Encoding")+2]
}
info<-cbind(n.read,l.read,encoding)
rownames(info)<-files
write.table(info,"info",sep="\t",quote=F)

#===========Unzip
for file in raw_data/*.gz
do
sbatch -A Eels -t 12:00:00 --wrap "gunzip $file"
done

