## Make files for selection scans and summarise results

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#-----------Copy phased chromosome files to selection folder
folder="/faststorage/project/Eels/eel_combined_Aja/selection"

mkdir ${folder}
mkdir ${folder}/chromosomes

cp /faststorage/project/Eels/eel_combined_Aja/phasing/eu+am/Chr_*.phased.vcf ${folder}/chromosomes

#-------------Compress vcf files
folder="/faststorage/project/Eels/eel_combined_Aja/selection/chromosomes"

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --job-name compress --wrap "bgzip ${folder}/${line}.phased.vcf"
done

#-------------Index
folder="/faststorage/project/Eels/eel_combined_Aja/selection/chromosomes"

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --job-name index --wrap "tabix ${folder}/${line}.phased.vcf.gz"
done

#-------------Subsetting scaffold files
input="/faststorage/project/Eels/eel_combined_Aja/selection/chromosomes"
american="/faststorage/project/Eels/eel_combined_Aja/selection/am"
european="/faststorage/project/Eels/eel_combined_Aja/selection/eu"
european_list="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_european.txt"
american_list="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_american.txt"

mkdir $american
mkdir $american/VCF
mkdir $european
mkdir $european/VCF

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --job-name subset --wrap\
 "bcftools view -S ${european_list} -O v -o ${european}/VCF/${line}_eu.vcf ${input}/${line}.phased.vcf.gz"
sbatch -A Eels -t 12:00:00 --job-name subset --wrap\
 "bcftools view -S ${admixed_list} -O v -o ${admixed}/VCF/${line}_adm.vcf ${input}/${line}.phased.vcf.gz"
done

###Run XP-EHH, XP-nSL, XP-CLR and FST

###Proceed after running selection scans

#check for intersection between XP-EHH, XP-nSL, XP-CLR and FST, and keep only regions that are common to at least two methods for each population
folder="/faststorage/project/Eels/eel_combined_Aja"

multiIntersectBed -i ${folder}/selection/eu/cr.99.XP-CLR.25kb.EU.txt ${folder}/selection/eu/cr.99.XP-EHH.norm.10kb.EU.txt ${folder}/selection/eu/cr.99.XP-nSL.norm.10kb.EU.txt ${folder}/Fst.dxy/windowStats_w10000-s1000-m20/cr.99.Fst.txt -names XP-CLR XP-EHH XP-nSL FST | awk '$4 >= 2' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $5}' > ${folder}/selection/eu/cr.2.99.10kb.EU.XP-Fst.bed

#check for intersection between iHS and nSL for european eel, and keep only regions that are common to both methods
folder="/faststorage/project/Eels/eel_combined_Aja"

multiIntersectBed -i ${folder}/selection/eu/cr.99.iHS.norm.10kb.EU.txt ${folder}/selection/eu/cr.99.nSL.norm.10kb.EU.txt -names iHS nSL | awk '$4 >= 2' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $5}' > ${folder}/selection/eu/cr.2.99.10kb.EU.iHS-nSL.bed

#concatenate outlier regions from XP-EHH, XP-nSL, XP-CLR and Fst in European eels with those from iHS and nSL
cat "/faststorage/project/Eels/eel_combined_Aja/selection/eu/cr.2.99.10kb.EU.iHS-nSL.bed" "/faststorage/project/Eels/eel_combined_Aja/selection/eu/cr.2.99.10kb.EU.XP-Fst.bed" > "/faststorage/project/Eels/eel_combined_Aja/selection/eu/cr.2.99.10kb.EU.bed"

#merge overlapping and adjecent regions between runs for each population
folder="/faststorage/project/Eels/eel_combined_Aja/selection"

bedtools sort -i ${folder}/eu/cr.2.99.10kb.EU.bed | bedtools merge -d 1 -c 4 -o collapse > ${folder}/eu/cr.2.99.10kb.merged.bed #merge overlapping regions
bedtools sort -i ${folder}/am/cr.2.99.10kb.AM.bed | bedtools merge -d 1 -c 4 -o collapse > ${folder}/am/cr.2.99.10kb.merged.bed #merge overlapping regions

# extract unique values in column 4 (which methods detected a region) (R)
file.am<-"/faststorage/project/Eels/eel_combined_Aja/selection/am/cr.2.99.10kb.merged.bed"
file.eu<-"/faststorage/project/Eels/eel_combined_Aja/selection/eu/cr.2.99.10kb.merged.bed"

df.am<-read.table(file.am, colClasses=c('character','integer','integer','character'))
df.eu<-read.table(file.eu, colClasses=c('character','integer','integer','character'))

df.am$V4 <- sapply(strsplit(df.am$V4, ','), function(x) toString(unique(x)))
df.eu$V4 <- sapply(strsplit(df.eu$V4, ','), function(x) toString(unique(x)))

write.table(df.am, file.am, row.names=FALSE, col.names=FALSE,sep='\t',quote=FALSE)
write.table(df.eu, file.eu, row.names=FALSE, col.names=FALSE,sep='\t',quote=FALSE)

#check for intersection between candidate regions for american and european eel
folder="/faststorage/project/Eels/eel_combined_Aja/selection/"

bedtools intersect -a ${folder}/am/cr.2.99.10kb.merged.bed -b ${folder}/eu/cr.2.99.10kb.merged.bed -wb -wa -loj | sort -k 1,1n -k2,2n > ${folder}/common/cr.2.99.10kb_intersect.bed

#==============Prepare transcript file
tail -n +2 "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.txt" > "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"

sed -i 's/Chr_01/01/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_02/02/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_03/03/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_04/04/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_05/05/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_06/06/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_07/07/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_08/08/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_09/09/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_10/10/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_11/11/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_12/12/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_13/13/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_14/14/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_15/15/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_16/16/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_17/17/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_18/18/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"
sed -i 's/Chr_19/19/g' "/faststorage/project/Eels/eel_combined_Aja/selection/transcripts2.bed"

#==============Check if regions under selection intersect with genes
folder="/faststorage/project/Eels/eel_combined_Aja/selection/"

bedtools intersect -a  ${folder}/am/cr.2.99.10kb.merged.bed -b ${folder}/transcripts2.bed -wb -wa -loj | sort -k 1,1n -k2,2n > ${folder}/am/cr.2.99.10kb_intersect.bed
bedtools intersect -a  ${folder}/eu/cr.2.99.10kb.merged.bed -b ${folder}/transcripts2.bed -wb -wa -loj | sort -k 1,1n -k2,2n > ${folder}/eu/cr.2.99.10kb_intersect.bed

#--------------extract list of gene names (for topGO analysis), list of candidate regions with gene names (for individual population summary tables) and list of gene regions (for collective summary table) for each population
folder="/faststorage/project/Eels/eel_combined_Aja/selection"
GO="/faststorage/project/Eels/eel_combined_Aja/GOterm"

mkdir $GO
mkdir $GO/am
mkdir $GO/eu

cut -f 9 ${folder}/am/cr.2.99.10kb_intersect.bed | sort | uniq > ${GO}/am/am.2.99.10kb_genes.bed # gene list for topGO
cut -f 1-4,9 ${folder}/am/cr.2.99.10kb_intersect.bed > ${folder}/am/am.2.99.10kb_cr+genes.bed  # candidate regions with gene names
cut -f 5-7,9,12-14 ${folder}/am/cr.2.99.10kb_intersect.bed | sort | uniq > ${folder}/am/am.2.99.10kb_gr.bed  # gene regions with pop name column

cut -f 9 ${folder}/eu/cr.2.99.10kb_intersect.bed | sort | uniq > ${GO}/eu/eu.2.99.10kb_genes.bed # gene list for topGO
cut -f 1-4,9 ${folder}/eu/cr.2.99.10kb_intersect.bed > ${folder}/eu/eu.2.99.10kb_cr+genes.bed  # candidate regions with gene names
cut -f 5-7,9,12-14 ${folder}/eu/cr.2.99.10kb_intersect.bed | sort | uniq > ${folder}/eu/eu.2.99.10kb_gr.bed  # gene regions with pop name column

#-----------collapse identical candidate regions and condenses unique values (genes) (R)
setwd("/faststorage/project/Eels/eel_combined_Aja/selection")

library(tidyverse)

pops<- c("am","eu")

for (pop in pops) {
file<-paste(pop, "/", pop, ".2.99.10kb_cr+genes.bed", sep="")
out<-paste(pop, "/", pop, ".2.99.10kb_cr+genes.condensed.bed", sep="")

df <- read.table(file, header=F, quote="", sep="\t",colClasses="character")

new <- df %>%
  group_by(V1,V2,V3,V4) %>%
  summarise_all(~paste(unique(V5), collapse = ", "))

write.table(new, out,col.names=FALSE,row.names=FALSE,quote=F,sep = "\t")
}

#===========Plot XP-EHH, XP-nSL, XP-CLR, FST and dXY
folder="/faststorage/project/Eels/eel_combined_Aja"

xpehh=${folder}/XP-EHH/XP-EHH.norm.txt
xpnsl=${folder}/XP-nSL/XP-nSL.norm.txt
xpclr1=${folder}/XP-CLR/am/XP-CLR.25kb.am.txt
xpclr2=${folder}/XP-CLR/eu/XP-CLR.25kb.eu.txt
fst=${folder}/Fst.dxy/windowStats_w10000-s1000-m20/Fst.Dxy.pi.csv
plot=${folder}/selection/selection.tiff
    
sbatch -A Coregonus -t 12:00:00 --mem 32G --job-name plot_selection --wrap\
 "Rscript ${folder}/selection/plot_selection.compressed.r AM EU ${xpehh} ${xpnsl} ${xpclr1} ${xpclr2} ${fst} ${plot}"
 
#===========Plot iHS and nSL for European eel
folder="/faststorage/project/Eels/eel_combined_Aja"

ihs=${folder}/iHS/iHS.norm.txt
nsl=${folder}/nSL/nSL.norm.txt
plot=${folder}/selection/selection_iHS-nSL.tiff
    
sbatch -A Coregonus -t 12:00:00 --mem 32G --job-name plot_selection --wrap\
 "Rscript ${folder}/selection/plot_iHS-nSL.compressed.r ${ihs} ${nsl} ${plot}"


