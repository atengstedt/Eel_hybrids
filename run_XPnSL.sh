## Perform selection scans for Eels: XP-nSL

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 3 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#------------Run XP-nSL
input="/faststorage/project/Eels/eel_combined_Aja/selection"
output="/faststorage/project/Eels/eel_combined_Aja/XP-nSL"

mkdir ${output}

cat chromosome_names.txt | while read line
do
am=${input}/am/VCF/${line}_am.vcf
eu=${input}/eu/VCF/${line}_eu.vcf
out=${output}/${line}.XP-nSL
            
sbatch -A Eels -t 12:00:00 --job-name XP-nSL --mem 24G --wrap\
 "./selscan-linux-2.0.0/selscan --xpnsl --pmap --vcf ${am} --vcf-ref ${eu} --out ${out}"
done

#============normalize XP-nSL and identify clusters of extreme scores
cd "/faststorage/project/Eels/eel_combined_Aja/XP-nSL/"
sbatch -A Eels -t 24:00:00 --job-name norm_XPEHH --wrap\
 "/faststorage/project/Eels/eel_combined_Aja/selscan-linux-2.0.0/norm --xpnsl --files *xpnsl.out"
# "/faststorage/project/Eels/eel_combined_Aja/selscan-linux-2.0.0/norm --xpnsl --bp-win --winsize 400000 --files *xpnsl.out"

#------------replace random pos ID with chr name (I'm replacing with chromosome numbers for plotting)
folder="/faststorage/project/Eels/eel_combined_Aja/XP-nSL"

for chrom in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
  do
  awk -v dt="${chrom}" 'BEGIN{FS=OFS="\t"}NR>1{$1=dt}1' ${folder}/Chr_${chrom}.XP-nSL.xpnsl.out.norm > ${folder}/tmp.${chrom} && mv ${folder}/tmp.${chrom} ${folder}/Chr_${chrom}.XP-nSL.xpnsl.out.norm
done

#------------merge all normalized result files
folder="/faststorage/project/Eels/eel_combined_Aja/XP-nSL/"

cp ${folder}/Chr_01.XP-nSL.xpnsl.out.norm ${folder}/XP-nSL.norm.txt

for i in {02..19}; do
  awk 'NR>1 {print}' ${folder}/Chr_${i}.XP-nSL.xpnsl.out.norm >> ${folder}/XP-nSL.norm.txt
done

# replace scientific notation of large numbers
#sed -i 's/1e+06/1000000/g' "/faststorage/project/Eels/eel_combined_Aja/XP-nSL/XP-nSL.norm.txt"

#------------plot 
folder="/faststorage/project/Eels/eel_combined_Aja/XP-nSL/"

cd ${folder}
xpnsl=${folder}/XP-nSL.norm.txt
plot=${folder}/XP-nSL.norm.tiff

sbatch -A Eels -t 12:00:00 --job-name plot_XPNSL --mem 44G --wrap\
 "Rscript plot_XPNSL.compressed.r ${xpnsl} ${plot}"

#-----------Divide the normalized result files in positive and negative normxpehh values
folder="/faststorage/project/Eels/eel_combined_Aja/XP-nSL"

awk '$9 > 0' ${folder}/XP-nSL.norm.txt | sed 1d | awk -F"\t" '{print $1 "\t" $2 "\t" $9}' > ${folder}/XP-nSL.norm.AM.txt
awk '$9 < 0' ${folder}/XP-nSL.norm.txt | awk -F"\t" '{print $1 "\t" $2 "\t" $9}' > ${folder}/XP-nSL.norm.EU.txt

#-----------Extract candidate regions
folder="/faststorage/project/Eels/eel_combined_Aja/"

sbatch -A Coregonus -t 24:00:00 --job-name calc_cand --mem 4G --wrap\
 "Rscript ${folder}/selection/calc_candidate_regions.r ${folder}/XP-nSL/XP-nSL.norm.AM.txt 10000 1000 0.99 ${folder}/selection/am/cr.99.XP-nSL.norm.10kb.AM.txt"
sbatch -A Coregonus -t 24:00:00 --job-name calc_cand --mem 4G --wrap\
 "Rscript ${folder}/selection/calc_candidate_regions.r ${folder}/XP-nSL/XP-nSL.norm.EU.txt 10000 1000 0.99 ${folder}/selection/eu/cr.99.XP-nSL.norm.10kb.EU.txt"
