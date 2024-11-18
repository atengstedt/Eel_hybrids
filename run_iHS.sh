## Perform selection scans for Eels: iHS

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 3 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#------------Run iHS
input="/faststorage/project/Eels/eel_combined_Aja/selection"
output="/faststorage/project/Eels/eel_combined_Aja/iHS"

mkdir ${output}

cat chromosome_names.txt | while read line
do
eu=${input}/eu/VCF/${line}_eu.vcf
out=${output}/${line}.iHS
            
sbatch -A Eels -t 12:00:00 --job-name iHS --mem 24G --wrap\
 "./selscan-linux-2.0.0/selscan --ihs --pmap --vcf ${eu} --out ${out}"
done

#============normalize iHS
cd "/faststorage/project/Eels/eel_combined_Aja/iHS" 

sbatch -A Eels -t 24:00:00 --job-name norm_iHS --wrap\
 "../selscan-linux-2.0.0/norm --ihs --files *.ihs.out"

#------------replace random pos ID with chr name (I'm replacing with chromosome numbers for plotting)
folder="/faststorage/project/Eels/eel_combined_Aja/iHS"

for chrom in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
  do
  awk -v dt="${chrom}" 'BEGIN{FS=OFS="\t"}NR>0{$1=dt}1' ${folder}/Chr_${chrom}.iHS.ihs.out.100bins.norm > ${folder}/tmp.${chrom} && mv ${folder}/tmp.${chrom} ${folder}/Chr_${chrom}.iHS.ihs.out.100bins.norm
done

#------------merge all normalized result files
folder="/faststorage/project/Eels/eel_combined_Aja/iHS"

cp ${folder}/Chr_01.iHS.ihs.out.100bins.norm ${folder}/iHS.norm.txt

for i in {02..19}; do
  awk 'NR>0 {print}' ${folder}/Chr_${i}.iHS.ihs.out.100bins.norm >> ${folder}/iHS.norm.txt
done

#------------plot 
folder="/faststorage/project/Eels/eel_combined_Aja/iHS"

cd ${folder}
xpehh=${folder}/iHS.norm.txt
plot=${folder}/iHS.norm.tiff

sbatch -A Eels -t 12:00:00 --job-name plot_XPEHH --mem 44G --wrap\
 "Rscript plot_iHS.compressed.r ${xpehh} ${plot}"
 
#-----------Extract candidate regions
folder="/faststorage/project/Eels/eel_combined_Aja/"

sbatch -A Coregonus -t 24:00:00 --job-name calc_cand --mem 16G --wrap\
 "Rscript ${folder}/selection/calc_candidate_regions.iHS-nSL.r ${folder}/iHS/iHS.norm.txt 10000 1000 0.99 ${folder}/iHS/iHS.norm.10kb.txt ${folder}/selection/eu/cr.99.iHS.norm.10kb.EU.txt"
