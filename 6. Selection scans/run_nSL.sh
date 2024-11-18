## Perform selection scans for Eels: nSL

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 3 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

#------------Run nSL
input="/faststorage/project/Eels/eel_combined_Aja/selection"
output="/faststorage/project/Eels/eel_combined_Aja/nSL"

mkdir ${output}

cat chromosome_names.txt | while read line
do
eu=${input}/eu/VCF/${line}_eu.vcf
out=${output}/${line}.nSL.maf0
            
sbatch -A Eels -t 12:00:00 --job-name nSL --mem 24G --wrap\
 "./selscan-linux-2.0.0/selscan --nsl --maf 0 --pmap --vcf ${eu} --out ${out}"
done

#============normalize nSL 
cd "/faststorage/project/Eels/eel_combined_Aja/nSL" 

sbatch -A Eels -t 24:00:00 --job-name norm_iHS --wrap\
 "../selscan-linux-2.0.0/norm --ihs --files *.nsl.out"

#------------replace random pos ID with chr name (I'm replacing with chromosome numbers for plotting)
folder="/faststorage/project/Eels/eel_combined_Aja/nSL"

for chrom in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
  do
  awk -v dt="${chrom}" 'BEGIN{FS=OFS="\t"}NR>0{$1=dt}1' ${folder}/Chr_${chrom}.nSL.nsl.out.100bins.norm > ${folder}/tmp.${chrom} && mv ${folder}/tmp.${chrom} ${folder}/Chr_${chrom}.nSL.nsl.out.100bins.norm
done

#------------merge all normalized result files
folder="/faststorage/project/Eels/eel_combined_Aja/nSL"

cp ${folder}/Chr_01.nSL.nsl.out.100bins.norm ${folder}/nSL.norm.txt

for i in {02..19}; do
  awk 'NR>0 {print}' ${folder}/Chr_${i}.nSL.nsl.out.100bins.norm >> ${folder}/nSL.norm.txt
done

#------------plot 
folder="/faststorage/project/Eels/eel_combined_Aja/nSL"

cd ${folder}
xpehh=${folder}/nSL.norm.txt
plot=${folder}/nSL.norm.tiff

sbatch -A Eels -t 12:00:00 --job-name plot_XPEHH --mem 44G --wrap\
 "Rscript plot_nSL.compressed.r ${xpehh} ${plot}"

#-----------Extract candidate regions
folder="/faststorage/project/Eels/eel_combined_Aja/"

sbatch -A Coregonus -t 24:00:00 --job-name calc_cand --mem 16G --wrap\
 "Rscript ${folder}/selection/calc_candidate_regions.iHS-nSL.r ${folder}/nSL/nSL.norm.txt 10000 1000 0.99 ${folder}/nSL/nSL.norm.10kb.txt ${folder}/selection/eu/cr.99.nSL.norm.10kb.EU.txt"
