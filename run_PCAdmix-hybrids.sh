## Local ancestry inference using PCAdmix (non-phased data for admixed individuals, phased data for ancestral reference populations)

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#============Subset eu+am scaffold files into by species
mkdir "/faststorage/project/Eels/eel_combined_Aja/PCAdmix"
mkdir "/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids"

ancestral="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am"
european="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/european_phased"
american="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/american_phased"
amlist="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_american.txt"
eulist="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_european.txt"

mkdir ${european}
mkdir ${american}

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${ancestral}/${line}.phased.vcf --keep ${eulist} --recode --out ${european}/AA_${line}"
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${ancestral}/${line}.phased.vcf --keep ${amlist} --recode --out ${american}/AR_${line}"
done

#-------------Convert to beagle 3 format
vcf2beagle="/faststorage/project/Eels/eel_genome_Aja/PCAdmix_2021/vcf2beagle.jar"
admixed="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/admixed_nonphased"
european="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/european_phased"
american="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/american_phased"

mkdir ${admixed}

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap\
 "cat /faststorage/project/Eels/eel_combined_Aja/phasing/admixed/${line}.vcf | java -jar $vcf2beagle m ${admixed}/admixed_${line}.pcadmix"
sbatch -A Eels -t 12:00:00 --wrap\
 "cat ${american}/AR_${line}.recode.vcf | java -jar $vcf2beagle m ${american}/AR_${line}.pcadmix"
sbatch -A Eels -t 12:00:00 --wrap\
 "cat ${european}/AA_${line}.recode.vcf | java -jar $vcf2beagle m ${european}/AA_${line}.pcadmix"
done

#-------------unzip
sbatch -A Eels -t 12:00:00 --wrap\
 "gunzip /faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/*phased/*bgl.gz"
 
#=============Make map files
admixed="/faststorage/project/Eels/eel_combined_Aja/phasing/admixed/"
maps="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/maps"

mkdir ${maps}

cat ./chromosome_names.txt | while read line
do
awk '{print $1"\t"$1"_"$2"\t"0"\t"$2}' ${admixed}/$line.ordered.vcf | sed '/#/d' > ${maps}/$line.map
done

#=============Run PCAdmix
american="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/american_phased/"
european="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/european_phased/"
admixed="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/admixed_nonphased/"
maps="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/maps/"
output="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/wMb0.01-maf0.10-ld0"

mkdir $output

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 1:00:00 -p express --job-name pcadmix --mem 16G --wrap\
 "/faststorage/project/Eels/eel_genome_Aja/PCAdmix_2021/PCAdmix3_linux -anc ${european}/AA_$line.pcadmix.bgl ${american}/AR_$line.pcadmix.bgl -adm ${admixed}/admixed_$line.pcadmix.bgl -map ${maps}/$line.map -maf 0.10 -wMb 0.01 -ld 0 -o ${output}/$line"
done

#--------------Plot results
folder=/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/wMb0.01-maf0.05-ld0/
plot="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/visualize_PCAdmix.r"

mkdir ${folder}/plots

for i in {01..19}
do
sbatch -A Eels -t 1:00:00 --job-name plot --wrap "Rscript ${plot} ${i} ${folder}"
done
