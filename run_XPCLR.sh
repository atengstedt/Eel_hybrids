## Perform selection scans for Eels: XP-CLR

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 3 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#============RUN XP-CLR per scaffold for both EU-AM and AM-EU
input="/faststorage/project/Eels/eel_combined_Aja/selection/chromosomes"
output="/faststorage/project/Eels/eel_combined_Aja/XP-CLR"
american="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_american.txt"
european="/faststorage/project/Eels/eel_combined_Aja/PCAdmix/Eel_list_european.txt"

mkdir ${output}
mkdir ${output}/am
mkdir ${output}/eu

cat chromosome_names.txt | while read line
do
sbatch -A Coregonus -t 72:00:00 --job-name XPCLR --mem 2G --wrap\
 "xpclr --input ${input}/${line}.phased.vcf --samplesA ${european} --samplesB ${american} --chr ${line} --phased --size 25000 --step 5000 --out ${output}/eu/${line}.25kb.XP-CLR.txt"
sbatch -A Coregonus -t 72:00:00 --job-name XPCLR --mem 2G --wrap\
 "xpclr --input ${input}/${line}.phased.vcf --samplesA ${american} --samplesB ${european} --chr ${line} --phased --size 25000 --step 5000 --out ${output}/am/${line}.25kb.XP-CLR.txt"
done
#------------replace random pos ID with chr name (I'm replacing with chromosome numbers for plotting)
folder="/faststorage/project/Eels/eel_combined_Aja/XP-CLR"

for chrom in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
  do
  awk -v dt="${chrom}" 'BEGIN{FS=OFS="\t"}NR>1{$2=dt}1' ${folder}/am/Chr_${chrom}.25kb.XP-CLR.txt > ${folder}/am/tmp.${chrom} && mv ${folder}/am/tmp.${chrom} ${folder}/am/Chr_${chrom}.25kb.XP-CLR.txt
  awk -v dt="${chrom}" 'BEGIN{FS=OFS="\t"}NR>1{$2=dt}1' ${folder}/eu/Chr_${chrom}.25kb.XP-CLR.txt > ${folder}/eu/tmp.${chrom} && mv ${folder}/eu/tmp.${chrom} ${folder}/eu/Chr_${chrom}.25kb.XP-CLR.txt
done

#-------------merge all scaffold result files
folder="/faststorage/project/Eels/eel_combined_Aja/XP-CLR"

cp ${folder}/am/Chr_01.25kb.XP-CLR.txt ${folder}/am/XP-CLR.25kb.am.txt
cp ${folder}/eu/Chr_01.25kb.XP-CLR.txt ${folder}/eu/XP-CLR.25kb.eu.txt

for i in 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19; do
  awk 'NR>1 {print}' ${folder}/am/Chr_${i}.25kb.XP-CLR.txt >> ${folder}/am/XP-CLR.25kb.am.txt
  awk 'NR>1 {print}' ${folder}/eu/Chr_${i}.25kb.XP-CLR.txt >> ${folder}/eu/XP-CLR.25kb.eu.txt
done

#-------------plot 
folder="/faststorage/project/Eels/eel_combined_Aja/XP-CLR"
xpclr1=${folder}/am/XP-CLR.25kb.am.txt
xpclr2=${folder}/eu/XP-CLR.25kb.eu.txt
plot1=${folder}/am/XP-CLR.25kb.am.tiff
plot2=${folder}/eu/XP-CLR.25kb.eu.tiff
    
sbatch -A Coregonus -t 12:00:00 --job-name plot_XPCLR --mem 8G --wrap\
 "Rscript ${folder}/plot_XPCLR.compressed.r am eu ${xpclr1} ${plot1}"
sbatch -A Coregonus -t 12:00:00 --job-name plot_XPCLR --mem 8G --wrap\
 "Rscript ${folder}/plot_XPCLR.compressed.r eu am ${xpclr2} ${plot2}"

#==============Extract candidate regions
folder="/faststorage/project/Eels/eel_combined_Aja/XP-CLR"
output="/faststorage/project/Eels/eel_combined_Aja/selection"

sbatch -A Coregonus -t 12:00:00 --job-name extract_top --mem 8G --wrap\
 "Rscript ${output}/calc_candidate_regions.XP-CLR.r am eu ${folder}/am/XP-CLR.25kb.am.txt 0.99 ${output}/am/cr.99.XP-CLR.25kb.AM.txt"
sbatch -A Coregonus -t 12:00:00 --job-name extract_top --mem 8G --wrap\
 "Rscript ${output}/calc_candidate_regions.XP-CLR.r eu am ${folder}/eu/XP-CLR.25kb.eu.txt 0.99 ${output}/eu/cr.99.XP-CLR.25kb.EU.txt"
