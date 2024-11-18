# Calling of variants

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#create a list of absolute paths to all bam files
ls /faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/mem/*.bam > "/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/mem/bamlist.txt"

#create a list of scaffold names
grep -o -E "^>\w+" "/faststorage/project/Eels/eel_combined_Aja/ref/GCF_013347855_1.fa" | tr -d ">" | tee scaffold_names.txt

#===========Call variants for each scaffold
mkdir "/faststorage/project/Eels/eel_combined_Aja/VCF"

ref="/faststorage/project/Eels/eel_combined_Aja/ref/GCF_013347855_1.fa"
from="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/mem/bamlist.txt"
to="/faststorage/project/Eels/eel_combined_Aja/VCF/raw_scaffolds/"

mkdir $to

cat scaffold_names.txt | while read line 
do
sbatch -A Eels -t 168:00:00 --job-name call_variants --wrap\
 "bcftools mpileup -Ou -r $line -q 20 --annotate FORMAT/DP,FORMAT/AD -f $ref -b ${from} | bcftools call -f GQ -v -m -O v -o ${to}/${line}.vcf"
done

#===========Merge scaffold VCF files
folder="/faststorage/project/Eels/eel_combined_Aja/VCF/raw_scaffolds"

sbatch -A Eels -t 24:00:00 -c 8 --mem 24G --job-name merge --wrap\
 "bcftools concat --threads 8 -f ${folder}/z.txt -O v -o ${folder}/Eels_raw.vcf"
 
#-----------Line count (how many variants remain?)
sbatch -A Eels -t 12:00:00 --job-name line_count --wrap\
 "bcftools view -H /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_raw.vcf | wc -l"

#-----------Rename individuals in header
file="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/individuals.txt"
input="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/VCF/Eels_raw.vcf"
output="/faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/VCF/Eels_raw.reheader.vcf"

sbatch -A Coregonus -t 12:00:00 --wrap\
 "bcftools reheader -s $file $input > $output"

#mv $output $input

#===========Inspect SNP depth distribution
awk '($0!~/^#/)&&($8~/^DP/){split($8,a,"[=;]");print a[2]}'\
 /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_raw.vcf > /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_raw.snp.depth

#-----------Plot (R)
a<-scan("/faststorage/project/Eels/eel_combined_Aja/VCF/Eels_raw.snp.depth")

tiff(file="/faststorage/project/Eels/eel_combined_Aja/VCF/Eels_raw.SNP_depth.tiff", units="cm", width = 12, height = 10, res=300)
par(mar=c(2,4,2,2)+.1,mgp=c(2,0.5,0), tck = -0.01, cex=0.6)
hist(a[a<2000],100,main="SNP depth distribution",xlab="SNP depth")
dev.off()

#===========Filter (indels, monomorphic, depth, MapQ)
sbatch -A Eels -t 12:00:00 --job-name vcf.sh --wrap\
 "vcfutils.pl varFilter -Q 20 -d 900 -D 1500 /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_raw.vcf | vcftools --vcf - --remove-indels --maf 0.0001 --recode --recode-INFO-all --out /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_filtered"

#-----------Line count (how many variants remain?)
sbatch -A Eels -t 12:00:00 --job-name line_count --wrap\
 "bcftools view -H /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_filtered.recode.vcf | wc -l"
 
#-------------Annotate SNP IDs
sbatch -A Eels -t 12:00:00 --job-name annotate --wrap\
 "bcftools annotate --set-id '%CHROM\_%POS' -O v -o /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_filtered.ann.vcf /faststorage/project/Eels/eel_combined_Aja/VCF/Eels_filtered.recode.vcf"

#------------Filter VCF file to remove singletons and doubletons, alleles that are not biallelic and SNPs with missing data
input="/faststorage/project/Eels/eel_combined_Aja/VCF/Eels_filtered.ann.vcf"
output="/faststorage/project/Eels/eel_combined_Aja/VCF/Eels_filtered.ann.mac3.max2.miss1"

sbatch -A Eels -t 12:00:00 --job-name filter --wrap\
 "vcftools --vcf ${input} --mac 3 --max-alleles 2 --max-missing 1 --recode --recode-INFO-all --out ${output}"
