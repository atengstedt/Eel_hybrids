## FST scan

##########################
# Check job info in batch.
for id in `ls *.out | cut -d "-" -f 2 | cut -d "." -f 1`
do
jobinfo $id | grep State
done | sort | uniq -c
##########################

### Run in environment thesis

mkdir Fst.dxy

#----------Calc dxy in sliding windows (speciationgenomics) from files with monomorphic sites
# convert VCF to geno.gz
input="/faststorage/project/Eels/eel_combined_Aja/VCF/for_HO_2/"
output="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy/geno-files"

mkdir ${output}

cat ./chromosome_names.txt | while read line 
do
sbatch -A Eels -t 12:00:00 --job-name vcf2geno --wrap\
 "vcftools --vcf ${input}/${line}.vcf --recode --stdout | ~/miniconda3/envs/thesis/bin/python ./Fst.dxy/parseVCF.py | ~/miniconda3/envs/thesis/bin/bgzip > ${output}/${line}.geno.gz"
done

# calculate pi, dxy and Fst
input="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy/geno-files"
output="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy/windowStats_w10000-s1000-m20/"

mkdir $output

cat ./chromosome_names.txt | while read line 
do
sbatch -A Eels -t 12:00:00 -c 6 --wrap\
 "python ./Fst.dxy/popgenWindows.py -g ${input}/${line}.geno.gz -o ${output}/${line}.Fst.Dxy.pi.csv -f phased -w 10000 -m 20 -s 1000 -p AR -p AA --popsFile ${input}/pop.info -T 6"
done

#------------replace random pos ID with chr name (I'm replacing with chromosome numbers for plotting)
folder="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy/windowStats_w10000-s1000-m20"

for chrom in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
  do
  awk -v dt="${chrom}" 'BEGIN{FS=OFS=","}NR>1{$1=dt}1' ${folder}/Chr_${chrom}.Fst.Dxy.pi.csv > ${folder}/tmp.${chrom} && mv ${folder}/tmp.${chrom} ${folder}/Chr_${chrom}.Fst.Dxy.pi.csv
done

# merge all chromosome files
folder="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy"

cp ${folder}/windowStats_w10000-s1000-m20/Chr_01.Fst.Dxy.pi.csv ${folder}/windowStats_w10000-s1000-m20/Fst.Dxy.pi.csv
    
for i in 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19; do
awk 'NR>1 {print}' ${folder}/windowStats_w10000-s1000-m20/Chr_${i}.Fst.Dxy.pi.csv >> ${folder}/windowStats_w10000-s1000-m20/Fst.Dxy.pi.csv
done

#Plot 
folder="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy"

stats=${folder}/windowStats_w10000-s1000-m20/Fst.Dxy.pi.csv
plot1=${folder}/windowStats_w10000-s1000-m20/Fst.tiff
plot2=${folder}/windowStats_w10000-s1000-m20/dxy.tiff
    
sbatch -A Eels -t 12:00:00 --job-name plot_Fst --wrap\
 "Rscript ${folder}/plot_FST-dxy.compressed.r am eu ${stats} ${plot1} ${plot2}"


#-----------Extract candidate regions
folder="/faststorage/project/Eels/eel_combined_Aja/Fst.dxy"
output="/faststorage/project/Eels/eel_combined_Aja/selection"

stats=${folder}/windowStats_w10000-s1000-m20/Fst.Dxy.pi.csv
cr1=${folder}/windowStats_w10000-s1000-m20/cr.95.Fst.txt
cr2=${folder}/windowStats_w10000-s1000-m20/cr.95.dxy.txt
    
sbatch -A Coregonus -t 12:00:00 --job-name extract_top --mem 8G --wrap \
 "Rscript ${output}/calc_candidate_regions.FST.r am eu ${stats} 0.95 ${cr1} ${cr2}"
 
