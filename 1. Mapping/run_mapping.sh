## Mapping of reads

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#===========Index reference genome.
sbatch -A Eels-t 12:00:00 --wrap "bwa index /faststorage/project/Eels/eel_combined_Aja/ref/GCF_013347855.1.fa"

#===========Map (BWA mem) for loop.
ref="/faststorage/project/Eels/eel_combined_Aja/ref/GCF_013347855.1.fa"
from="/faststorage/project/Eels/eel_combined_Aja/raw_data/"
to="/faststorage/project/Eels/eel_combined_Aja/mem/"

mkdir $to

for NameRoot in `ls $from | grep "fq$" | cut -d _ -f 1,2 | sort | uniq`
do
  sbatch -A Eels -t 24:00:00 --mem 16G -c 16 --wrap\
  "bwa mem -M -t 16 $ref ${from}/${NameRoot}_R1.fq ${from}/${NameRoot}_R2.fq > ${to}/${NameRoot}.sam"
done

#===========Convert SAM to BAM, clean up and sort
folder="/faststorage/project/Eels/eel_combined_Aja/mem"

mkdir ${folder}/tmp1
mkdir ${folder}/tmp2

for NameRoot in `ls $folder | grep "sam$" | cut -d "." -f 1 | sort | uniq`
do
sbatch -A Eels -t 24:00:00 --job-name samtools --mem 16G --wrap\
 "samtools sort -n -T ${folder}/tmp1/${NameRoot} ${folder}/${NameRoot}.sam | \
 samtools fixmate - - -m | \
 samtools sort - -T ${folder}/tmp2/${NameRoot} \
 -O bam -o ${folder}/${NameRoot}.bam"
done

#===========Index BAM files
folder="/faststorage/project/Eels/eel_combined_Aja/mem"

for NameRoot in `ls $folder | grep "bam$" | cut -d "." -f 1 | sort | uniq`
do
sbatch -A Eels -t 24:00:00 --job-name samtools --mem 16G --wrap\
 "samtools index -@ 16 ${folder}/${NameRoot}.bam"
done

#===========Stats
from="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/mem"
to="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/STATS_mem"

mkdir $to

for bam in `ls $from | grep ".bam$"`
do
sbatch -A Eels -t 12:00:00 --wrap "samtools flagstat ${from}/${bam} > ${to}/${bam/.bam/.flagstat}"
sbatch -A Eels -t 12:00:00 --wrap "samtools stats -c 1,1000,1 ${from}/${bam} > ${to}/${bam/.bam/.stats}"
done

#-----------Summary of flagstat.
echo -e "Indi\tTotal\tSec\tSup\tDup\tMapped\tPaired\tR1\tR2\tPro\tPMapped\tSingle"\
 > STATS_mem/zzflagstat
for file in `ls STATS_mem | grep "\\.flagstat$"`
do
echo -e -n ${file/.flagstat/}"\t" >> STATS_mem/zzflagstat
awk 'NR==1{printf($1);next}{printf("\t%d",$1)}END{printf("\n")}' STATS_mem/${file}\
 >> STATS_mem/zzflagstat
done

#-----------Summary (R).
# Run "STATS.r".
