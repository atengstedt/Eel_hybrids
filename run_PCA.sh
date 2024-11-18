#Principal component analysis on thinned VCF

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#============Thin VCF so no sites are within <1000 bp from each other
folder="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/VCF/"

sbatch -A Coregonus -t 12:00:00 --mem 8G --job-name vcftools --wrap\
 "vcftools --vcf ${folder}/chr01-19.filtered.ann.mac3.max2.miss1.recode.vcf --thin 1000 --recode --recode-INFO-all --out ${folder}/chr01-19.filtered.ann.mac3.max2.miss1.LD-thinned"
 
#835,600 SNPs (1 SNP per 1.000 bp)
 
#===========PCA (R)
f.vcf<-"/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/VCF/chr01-19.filtered.ann.mac3.max2.miss1.LD-thinned.recode.vcf" # Give the absolute pathname of the VCF subset.
f.popmap<-"/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/popmap_simple.txt"	# Give the absolute pathname of the popmap file.

#-----------Read the input files.
a<-read.table(f.vcf,sep="\t",stringsAsFactors=F)
popmap<-read.table(f.popmap,sep="\t",stringsAsFactors=F)

#-----------Parse the VCF file.
vcf<-as.matrix(a[,-c(1:9)])
nrow.vcf<-nrow(vcf)
ncol.vcf<-ncol(vcf)
chrom1<-as.integer(substring(vcf,1,1))
chrom2<-as.integer(substring(vcf,3,3))
chrom<-matrix(chrom1+chrom2,nrow.vcf,ncol.vcf)
colnames(chrom)<-popmap[,1]
rm(vcf,chrom1,chrom2)

#-----------Run PCA.
b<-prcomp(t(na.omit(chrom)))
e<-b$sdev^2
indis<-rownames(b$x)

#modify the dataframe before plotting

pca_data<-as.data.frame(b$x[,c(1,2,3,4,5,6)])  #save PC1-6
pca_data$ID <- row.names(pca_data)  #make row.numbers into a column

merged_data <- merge(pca_data, popmap, by.x = "ID", by.y = "V1", all = TRUE)  #merge PC information with popmap
names(merged_data) <- c("ID","PC1","PC2","PC3","PC4","PC5","PC6","GROUP")
write.table(merged_data, "/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/PCA_data.LD-thinned.txt", row.names=F)
#merged_data<-read.table("/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/PCA_data.txt", header=T)

#calculate percentage explained per PC
write.table(e, "/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/PCA_data.sdev.LD-thinned.txt", row.names=F)
#e<-read.table("/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/PCA_data.sdev.txt")

e[1]/sum(e)*100  #calculate percentage variance explained for PC1
#17.41719%
e[2]/sum(e)*100  #calculate percentage variance explained for PC2
#1.369538%

#make plot
library(ggplot2)

in2mm=25.4

pdf(file="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/PCAplot.LD-thinned.pdf", width = 80/in2mm, height = 80/in2mm)
par(fig=c(0,1,0,1), mar=c(3,3,2,2)+.1,mgp=c(2,0.5,0), tck = -0.01, cex=0.3)

ggplot(merged_data, aes(x = PC1/10, y = PC2/10, color = GROUP)) +
  geom_point(aes(shape = GROUP, fill = GROUP), size = 4) +
  theme_bw() +
  labs(x = "PC1 (17.4%)", y = "PC2 (1.4%)") +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  scale_shape_manual(values = c("AM" = 21, "EU" = 22, "Hybrid" = 24)) +  # Different shapes for each group
  scale_color_manual(values = c("black", "black", "black")) +
  scale_fill_manual(values = c("AM" = "#0072B2", "EU" = "#D55E00", "Hybrid" = "#999999")) +
  geom_text(data = subset(merged_data, ID == "Bur_05"), aes(label = "*"), size = 5, vjust = +0.8) +  #label one new hybrid
  guides(color = FALSE, fill = FALSE, shape = FALSE)
dev.off()

# Plot the eigenvalues.
dev.new() 
barplot(e,main="Eigenvalues")
savePlot("/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/PCA/eigenvalues.LD-thinned.tiff","png")
dev.off()

