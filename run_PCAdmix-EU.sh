## Local ancestry inference using PCAdmix (non-phased admixed files, phased ancestral files) for europeans analyzed as admixed

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#============Subset eu+am scaffold files by species
input="/faststorage/project/Eels/eel_combined_Aja/phasing/eu+am/"
eu_list="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_european.txt"
am_list="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_american.txt"
eu="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/eu_phased/"
am="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/am_phased/"

mkdir ${eu}
mkdir ${am}

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${input}/${line}.phased.vcf --keep ${am_list} --recode --out ${am}/AR_${line}"
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${input}/${line}.phased.vcf --keep ${eu_list} --recode --out ${eu}/AA_${line}"
done

#=============Filter phased scaffolds (eu and am) (maf 0.05) in preparation for PCAdmix with alternative admixed groups
eu_input="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/eu_phased"
eu_output="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/eu_phased_maf0.05"
am_input="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/am_phased"
am_output="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/am_phased_maf0.05"

mkdir ${eu_output}
mkdir ${am_output}

cat ./chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${am_input}/AR_${line}.recode.vcf --maf 0.05 --recode --out ${am_output}/AR_${line}"
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${eu_input}/AA_${line}.recode.vcf --maf 0.05 --recode --out ${eu_output}/AA_${line}" 
done

#----------------Extract new admixed population and new european population (exclude new admix-population) from phased scaffolds
for population in Bur Ell Gir-Sar Rin Seb Sto Val
do
mkdir /faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}

input="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/eu_phased/"
list="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_${population}.txt"
pop_output="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/${population}_phased"
eu_output="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/european-no-${population}_phased"

mkdir ${pop_output}
mkdir ${eu_output}

cat ../chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${input}/AA_${line}.recode.vcf --keep ${list} --recode --out ${pop_output}/${population}_${line}"
sbatch -A Eels -t 12:00:00 --wrap\
 "vcftools --vcf ${input}/AA_${line}.recode.vcf --remove ${list} --recode --out ${eu_output}/AA_${line}"
done
done

#-------------Convert to beagle 3 format (american)
folder="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/am_phased/"

cat ../chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --job-name vcf2beagle --wrap\
 "cat ${folder}/AR_${line}.recode.vcf | java -jar /faststorage/project/Eels/eel_genome_Aja/PCAdmix_2021/vcf2beagle.jar m ${folder}/AR_${line}.pcadmix"
done

#-------------Convert to beagle 3 format (eu + new "admixed" pops)
pop_folder="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/${population}_phased/"
eu_folder="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/european-no-${population}_phased/"

for population in Bur Ell Gir-Sar Rin Seb Sto Val
do
cat ../chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --job-name vcf2beagle --wrap\
 "cat ${eu_folder}/AA_$line.recode.vcf | java -jar /faststorage/project/Eels/eel_genome_Aja/PCAdmix/vcf2beagle.jar m ${eu_folder}/AA_$line.pcadmix"
sbatch -A Eels -t 12:00:00 --job-name vcf2beagle --wrap\
 "cat ${pop_folder}/${population}_${line}.recode.vcf | java -jar /faststorage/project/Eels/eel_genome_Aja/PCAdmix/vcf2beagle.jar m ${pop_folder}/${population}_$line.pcadmix"   
done
done

#-------------unzip (american)
sbatch -A Eels -t 12:00:00 --job-name unzip --wrap\
 "gunzip /faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/*phased/*bgl.gz"

#-------------unzip (eu + new "admixed" pops)
for population in Bur Ell Gir-Sar Rin Seb Sto Val
do
sbatch -A Eels -t 12:00:00 --job-name unzip --wrap\
 "gunzip /faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/*phased/*bgl.gz"
done

#==============Run PCAdmix
for population in Bur Ell Gir-Sar Rin Seb Sto Val
do

out=/faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/wMb0.01_prune-off
mkdir ${out}

cat ../chromosome_names.txt | while read line
do
sbatch -A Eels -t 12:00:00 --job-name pcadmix --mem 24G --wrap\
 "/faststorage/project/Eels/eel_genome_Aja/PCAdmix_2021/PCAdmix3_linux -anc /faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/european-no-${population}_phased/AA_$line.pcadmix.bgl /faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/am_phased/AR_$line.pcadmix.bgl -adm /faststorage/project/Eels/eel_combined_Aja/PCAdmix/${population}/${population}_phased/${population}_$line.pcadmix.bgl -map /faststorage/project/Eels/eel_combined_Aja/PCAdmix/hybrids/maps/$line.map -prune 0 -wMb 0.01 -o ${out}/$line"
done
done

#===============Collect vit file from each population into one file per chromosome
cd "/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed"

run=wMb0.01_maf0.05_prune-off
mkdir ${run}

cat ../../chromosome_names.txt | while read id
do 
cat ../Bur/${run}/${id}.vit.txt ../Ell/${run}/${id}.vit.txt ../Gir-Sar/${run}/${id}.vit.txt ../Rin/${run}/${id}.vit.txt ../Seb/${run}/${id}.vit.txt ../Sto/${run}/${id}.vit.txt ../Val/${run}/${id}.vit.txt > ${run}/vit_original/${id}-europeans.vit.txt
done

#=============Plot results (one file per figure)
folder="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/wMb0.01_maf0.05_prune-off/vit_coverage/"
plot="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/european-as-admixed/visualize_wMb0.01_maf0.05.r"

mkdir ${folder}/plots

for i in {01..19}
do
sbatch -A Eels -t 1:00:00 --job-name plot --wrap "Rscript ${plot} ${i} ${folder}"
done
