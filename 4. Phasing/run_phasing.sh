# Statistical phasing of ancestral reference populations

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#==============Split VCF file into subsets (one with species and one with admixed) and scaffolds
input="/faststorage/project/Eels/eel_combined_Aja/VCF/chr01-19.filtered.ann.mac3.max2.miss1.recode.vcf"
ancestral="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am"
admixed="/faststorage/project/Eels/eel_combined_Aja/phasing/admixed"
filter="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_admixed.txt"

mkdir ${ancestral}
mkdir ${admixed}

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${input} --chr ${line} --remove ${filter} --recode --out ${ancestral}/${line}"
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${input} --chr ${line} --keep ${filter} --recode --out ${admixed}/${line}"
done

#=============Phase ancestral populations vcf file (per scaffold)
folder="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am"

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --mem 24G -c 6 --wrap\
 "/faststorage/project/Eels/eel_genome_Aja/phasing/bin/shapeit --input-vcf ${folder}/${line}.recode.vcf \
        --output-log ${folder}/${line}.phased.log \
        -O ${folder}/${line}.phased \
        --states 200 \
        --force \
        --thread 6"
done

#------------Convert HAPS/SAMPLE output to VCF
folder="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am"

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap \
 "/faststorage/project/Eels/eel_genome_Aja/phasing/bin/shapeit -convert --input-haps ${folder}/${line}.phased --output-vcf ${folder}/${line}.phased.vcf"
done

#------------Compress scaffold files
folder="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am/"

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap \
 "bgzip ${folder}/${line}.phased.vcf"
done

#------------Index scaffold files
folder="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am/"

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap \
 "tabix ${folder}/${line}.phased.vcf.gz"
done

#===========Merge scaffold VCF files
cd "/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am/"

sbatch -A Eels -t 24:00:00 -c 8 --mem 24G --job-name merge --wrap\
 "bcftools concat --threads 8 -f z.txt -O v -o /faststorage/project/Eels/eel_combined_Aja/eel_combined_Aja/VCF/chr01-19.filtered.ann.mac3.max2.miss1.phased.vcf"
