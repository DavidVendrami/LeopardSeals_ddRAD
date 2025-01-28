###############################################################
######## Leopard seal RADseq bioinformatic analyses ###########
###############################################################

##### 1. QCs - part 1
# Run fastqc and multiqc on all samples to check quality
fastqc *.gz
multiqc .

# All good. Let's also plot the number of reads per sample.
cut -d$'\t' -f1,5 QCs/multiqc_data/multiqc_fastqc.txt | sed 's/ /_/' > Reads_perSample.txt
Rscript Plot_reads_per_sample.r

# It appears there are 5 samples with less than 1M reads:
# 128: but luckily we have a good replicate (128b)!
# 93139: no rpelicate
# 105: no replicate
# 78140: no replicate (and basically no sequencing)
# empty: expected to be empty (negative control)
# Let's map them anyway but let's not use them for SNP calling.

##### 2. Mapping (submitted via sbatch)
#!/bin/bash
#SBATCH --job-name=Mapping.job
#SBATCH --output=Mapping.out
#SBATCH --error=Mapping.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

# index reference genome 
bwa index Hydrurga_leptonyx.fasta

# Then align. Pipe the output to samtools to convert to bam and sort. 
for i in /prj/furseal-genome/David/Leopard_Seal/Raw_reads/*.fastq.gz
do
gzip -d $i
/vol/biotools/bin/bwa mem -t 8 /prj/furseal-genome/David/Leopard_Seal/Reference/Hydrurga_leptonyx.fasta ${i%.gz} | \
/vol/biotools/bin/samtools view -b -h | /vol/biotools/bin/samtools sort -o Mappings/${i%.fq.gz}_sorted.bam
gzip ${i%.gz}
done

##### 3. QCs - part 2
cd Mappings

for i in *.bam 
do
samtools index $i 
samtools stats $i > QCs_part2/$i.samtools.stats
samtools idxstats $i > QCs_part2/$i.samtools.idxstats
samtools flagstat $i > QCs_part2/$i.samtools.flagstats
done

multiqc .

# Everything went as expected with really high mapping rate (~99%)
# The only exception is sample 96, which yielded much lower mapping rate (84.9%). Why? Hard to say. 
# Let's genotype it anyway to see if something funny emerges.

##### 4. Genotyping
# Move bad samples into a separate folder
mkdir bad_quality_samples
mv empty_1_sorted.bam bad_quality_samples/
mv 128_1_sorted.bam 
mv 128_1_sorted.bam bad_quality_samples/
mv 93139_1_sorted.bam bad_quality_samples/
mv 105_1_sorted.bam bad_quality_samples/
mv 78140_1_sorted.bam bad_quality_samples/

# Produce population map
ls -1 *.bam | sed 's/.bam/    1/g' > pop_map.txt

# Genotype using Stacks ref_map.pl (submitted as a slurm job)
#!/bin/bash
#SBATCH --job-name=Genotyping.job
#SBATCH --output=/prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ref_map_output/genotyping.out
#SBATCH --error=/prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ref_map_output/genotyping.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

# Genotype with Stacks ref_map.pl
/prj/mar-in-gen/bin/bin/ref_map.pl --samples /prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ --popmap /prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/pop_map.txt -o /prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ref_map_output -T 8 -X "populations: --vcf"

##### 5. Repeat genotyping with min MQ of 20
#!/bin/bash
#SBATCH --job-name=Genotyping.job
#SBATCH --output=/prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ref_map_output/genotyping.out
#SBATCH --error=/prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ref_map_output/genotyping.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8

# Genotype with Stacks ref_map.pl
/prj/mar-in-gen/bin/bin/ref_map.pl --samples /prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ --popmap /prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/pop_map.txt -o /prj/furseal-genome/David/Leopard_Seal/Raw_reads/Mappings/ref_map_output_MQ20 -T 8 -X "gstacks: --min-mapq 20" -X "populations: --vcf"
# Very small difference in number of genotyped SNPs, so no need to proceed with this?

##### 6. Filtering (interactively with slxterm)
# Let's explore a bit retained numbers of SNPs for GQ and DP = 5 and 10; and max-missing = 0.8 and 0.9
mkdir Expl_filtering

# Retain only bi-allelic SNPs and genotypes with coverage of at least 5 or 10
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 5 --minGQ 5 --recode --out Expl_filtering/Leop_Seal_Biall_GQDP5
vcftools --vcf populations.snps.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 10 --minGQ 10 --recode --out Expl_filtering/Leop_Seal_Biall_GQDP10

# Let's explore the number of retained SNPs on the two DP datasets and remove SNPs with missingness of 20% and 10%.
cd Expl_filtering
vcftools --vcf Leop_Seal_Biall_GQDP5.recode.vcf --max-missing 0.8 --out Leop_Seal_Biall_GQDP5_MD80 # Retained = 33,157; out of 465,523
vcftools --vcf Leop_Seal_Biall_GQDP5.recode.vcf --max-missing 0.9 --out Leop_Seal_Biall_GQDP5_MD90 # Retained = 26,905; out of 465,523
vcftools --vcf Leop_Seal_Biall_GQDP10.recode.vcf --max-missing 0.8 --out Leop_Seal_Biall_GQDP10_MD80 # Retained = 26,547; out of 465,523
vcftools --vcf Leop_Seal_Biall_GQDP10.recode.vcf --max-missing 0.9 --out Leop_Seal_Biall_GQDP10_MD90 # Retained = 20,525; out of 465,523

# Let's keep the DP10 files
vcftools --vcf Leop_Seal_Biall_GQDP10.recode.vcf --max-missing 0.8 --recode --out Leop_Seal_Biall_GQDP10_MD80 # Retained = 26,547; out of 465,523
vcftools --vcf Leop_Seal_Biall_GQDP10.recode.vcf --max-missing 0.9 --recode --out Leop_Seal_Biall_GQDP10_MD90 # Retained = 20,525; out of 465,523

# Let's filter out potentially paralogous loci (i.e.: too high coverage, > twice the mean)
vcftools --vcf Leop_Seal_Biall_GQDP10_MD80.recode.vcf --site-mean-depth --out Leop_Seal_Biall_GQDP10_MD80 # Threshold = 129.996
vcftools --vcf Leop_Seal_Biall_GQDP10_MD90.recode.vcf --site-mean-depth --out Leop_Seal_Biall_GQDP10_MD90 # Threshold = 149.7219
Rscript Plot_depth.r

vcftools --vcf Leop_Seal_Biall_GQDP10_MD80.recode.vcf --max-meanDP 129.996 --recode --out Leop_Seal_Biall_GQDP10_MD80_maxDP129 # Retained = 25,263; out of 26,547
vcftools --vcf Leop_Seal_Biall_GQDP10_MD90.recode.vcf --max-meanDP 149.7219 --recode --out Leop_Seal_Biall_GQDP10_MD90_maxDP149 # Retained = 19,997; out of 20,525

# MAF
vcftools --vcf Leop_Seal_Biall_GQDP10_MD80_maxDP129.recode.vcf --maf 0.01 --recode --out Leop_GQDP10_MD80_maxDP129_MAF01 # Retained = 13,448; out of 26,547
vcftools --vcf Leop_Seal_Biall_GQDP10_MD90_maxDP149.recode.vcf --maf 0.01 --recode --out Leop_GQDP10_MD90_maxDP149_MAF01 # Retained = 11,112; out of 20,525

# Convert to plink
/prj/mar-in-gen/bin/plink/plink --vcf Leop_GQDP10_MD80_maxDP129_MAF01.recode.vcf --recode --double-id --aec --out Leop_GQDP10_MD80_maxDP129_MAF01
/prj/mar-in-gen/bin/plink/plink --vcf Leop_GQDP10_MD90_maxDP149_MAF01.recode.vcf --recode --double-id --aec --out Leop_GQDP10_MD90_maxDP149_MAF01

# HWE (mild, aimed at removing potentially-still-present PCR duplicate-based artifacts)
/prj/mar-in-gen/bin/plink/plink --file Leop_GQDP10_MD80_maxDP129_MAF01 --hwe 0.001 midp --aec --make-bed --out Leop_GQDP10_MD80_maxDP129_MAF01_HWE001 # Retained = 13,306; out of 13,448
/prj/mar-in-gen/bin/plink/plink --file Leop_GQDP10_MD90_maxDP149_MAF01 --hwe 0.001 midp --aec --make-bed --out Leop_GQDP10_MD90_maxDP149_MAF01_HWE001 # Retained = 11,018; out of 11,112

# Let's quickly check if any obvious difference exist among the 4 datasets by running a PCA and calcualting rxys
cp Leop_GQDP10_MD80_maxDP129_MAF01_HWE001* Final_4/
cp Leop_GQDP10_MD90_maxDP149_MAF01_HWE001* Final_4/
/prj/mar-in-gen/bin/plink/plink --bfile Leop_GQDP10_MD80_maxDP129_MAF01_HWE001 --aec --recodeA --out Final_4/Leop_GQDP10_MD80_maxDP129_MAF01_HWE001_LD2
/prj/mar-in-gen/bin/plink/plink --bfile Leop_GQDP10_MD90_maxDP149_MAF01_HWE001 --aec --recodeA --out Final_4/Leop_GQDP10_MD90_maxDP149_MAF01_HWE001_LD2

# Relatedness
# See 'temp_script_to_check_relatedness.sh' script
# Uuuuh, it seems there are a lot of duplicates (18)! 
# Only few of them are technical replicates, most of them are unknown recaptures! Really cool.

# PCA
Rscript Preliminary_PCAs.r
# The three sets of outliers are the animals that were sequenced 3 times.
# It seems that missingness 90% do a better job in minimizing error rate.

####### 8. Remove duplicates (keep individual with higher genotyping rate) and re-run filters (check that there are no further duplicates)
# Let's use the DP10 - MD 90% dataset to choose the repicates with more SNPs as it has lowest eror rate for known replicate.
# error rates (based on only known replicate):
R
data<-read.table('Leop_GQDP5_MD90_maxDP_MAF01_HWE001_LD2.raw', h=T)
info<-data[,c(1:6)]
geno<-data[,-c(1:6)]
miss<-apply(geno,1,function(x) length(which(is.na(x))))
info$Miss<-miss
write.table(info,"Missingness.txt",quote=F,col.names=T,row.names=F,sep='\t')
# process it to retain one individual per replicate.

# Create two vcf files: one for analyses (without duplicates) and one only with the duplicates to later calculate genotyping error rate
cut -d '       ' -f1 ../Expl_filtering/Final_4/Missingness.txt > Inds_to_keep.txt
vcftools --vcf ../populations.snps.vcf --keep Inds_to_keep.txt --recode --out Leop_NoDups
vcftools --vcf ../populations.snps.vcf --keep Dups.txt --recode --out Leop_AllDups

# Now, re-filter everything:
# DP and GQ
vcftools --vcf Leop_NoDups.recode.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 5 --minGQ 5 --recode --out Leop_NoDups_Biall_GQDP5
vcftools --vcf Leop_NoDups.recode.vcf --min-alleles 2 --max-alleles 2 --remove-indels --minDP 10 --minGQ 10 --recode --out Leop_NoDups_Biall_GQDP10

# max-missing
vcftools --vcf Leop_NoDups_Biall_GQDP5.recode.vcf --max-missing 0.8 --recode --out Leop_NoDups_Biall_GQDP5_MD80 # Retained = 33,205; out of 465,523
vcftools --vcf Leop_NoDups_Biall_GQDP5.recode.vcf --max-missing 0.9 --recode --out Leop_NoDups_Biall_GQDP5_MD90 # Retained = 26,977; out of 465,523
vcftools --vcf Leop_NoDups_Biall_GQDP10.recode.vcf --max-missing 0.8 --recode --out Leop_NoDups_Biall_GQDP10_MD80 # Retained = 26,633; out of 465,523
vcftools --vcf Leop_NoDups_Biall_GQDP10.recode.vcf --max-missing 0.9 --recode --out Leop_NoDups_Biall_GQDP10_MD90 # Retained = 20,519; out of 465,523

# MaxDP
vcftools --vcf Leop_NoDups_Biall_GQDP5_MD80.recode.vcf --site-mean-depth --out Leop_NoDups_Biall_GQDP5_MD80 # Threshold = 110.937
vcftools --vcf Leop_NoDups_Biall_GQDP5_MD90.recode.vcf --site-mean-depth --out Leop_NoDups_Biall_GQDP5_MD90 # Threshold = 127.0952
vcftools --vcf Leop_NoDups_Biall_GQDP10_MD80.recode.vcf --site-mean-depth --out Leop_NoDups_Biall_GQDP10_MD80 # Threshold = 129.7845
vcftools --vcf Leop_NoDups_Biall_GQDP10_MD90.recode.vcf --site-mean-depth --out Leop_NoDups_Biall_GQDP10_MD90 # Threshold = 149.7095
vcftools --vcf Leop_NoDups_Biall_GQDP5_MD80.recode.vcf --max-meanDP 110.937 --recode --out Leop_NoDups_Biall_GQDP5_MD80_maxDP110 # Retained = 30,415
vcftools --vcf Leop_NoDups_Biall_GQDP5_MD90.recode.vcf --max-meanDP 127.0952 --recode --out Leop_NoDups_Biall_GQDP5_MD90_maxDP127 # Retained = 25,570
vcftools --vcf Leop_NoDups_Biall_GQDP10_MD80.recode.vcf --max-meanDP 129.7845 --recode --out Leop_NoDups_Biall_GQDP10_MD80_maxDP129 # Retained = 25,364
vcftools --vcf Leop_NoDups_Biall_GQDP10_MD90.recode.vcf --max-meanDP 149.7095 --recode --out Leop_NoDups_Biall_GQDP10_MD90_maxDP149 # Retained = 20,008

# MAF
vcftools --vcf Leop_NoDups_Biall_GQDP5_MD80_maxDP110.recode.vcf --maf 0.01 --recode --out Leop_NoDups_Biall_GQDP5_MD80_maxDP110_MAF01 # Retained = 15,006
vcftools --vcf Leop_NoDups_Biall_GQDP5_MD90_maxDP127.recode.vcf --maf 0.01 --recode --out Leop_NoDups_Biall_GQDP5_MD90_maxDP127_MAF01 # Retained = 13,082
vcftools --vcf Leop_NoDups_Biall_GQDP10_MD80_maxDP129.recode.vcf --maf 0.01 --recode --out Leop_NoDups_Biall_GQDP10_MD80_maxDP129_MAF01 # Retained = 12,854
vcftools --vcf Leop_NoDups_Biall_GQDP10_MD90_maxDP149.recode.vcf --maf 0.01 --recode --out Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01 # Retained = 10,596

# Convert to plink
/prj/mar-in-gen/bin/plink/plink --vcf Leop_NoDups_Biall_GQDP5_MD80_maxDP110_MAF01.recode.vcf --recode --double-id --aec --out Leop_NoDups_Biall_GQDP5_MD80_maxDP110_MAF01
/prj/mar-in-gen/bin/plink/plink --vcf Leop_NoDups_Biall_GQDP5_MD90_maxDP127_MAF01.recode.vcf --recode --double-id --aec --out Leop_NoDups_Biall_GQDP5_MD90_maxDP127_MAF01
/prj/mar-in-gen/bin/plink/plink --vcf Leop_NoDups_Biall_GQDP10_MD80_maxDP129_MAF01.recode.vcf --recode --double-id --aec --out Leop_NoDups_Biall_GQDP10_MD80_maxDP129_MAF01
/prj/mar-in-gen/bin/plink/plink --vcf Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01.recode.vcf --recode --double-id --aec --out Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01

# HWE
/prj/mar-in-gen/bin/plink/plink --file Leop_NoDups_Biall_GQDP5_MD80_maxDP110_MAF01 --hwe 0.001 midp --aec --make-bed --out Leop_NoDups_Biall_GQDP5_MD80_maxDP110_MAF01_HWE001 # Retained = 14,892
/prj/mar-in-gen/bin/plink/plink --file Leop_NoDups_Biall_GQDP5_MD90_maxDP127_MAF01 --hwe 0.001 midp --aec --make-bed --out Leop_NoDups_Biall_GQDP5_MD90_maxDP127_MAF01_HWE001 # Retained = 13,021
/prj/mar-in-gen/bin/plink/plink --file Leop_NoDups_Biall_GQDP10_MD80_maxDP129_MAF01 --hwe 0.001 midp --aec --make-bed --out Leop_NoDups_Biall_GQDP10_MD80_maxDP129_MAF01_HWE001 # Retained = 12,784
/prj/mar-in-gen/bin/plink/plink --file Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01 --hwe 0.001 midp --aec --make-bed --out Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01_HWE001 # Retained = 10,555

######### 9. Error rate -> choose dataset that minimizes error rate
# Make list of SNPs
/prj/mar-in-gen/bin/plink/plink --bfile Leop_NoDups_Biall_GQDP5_MD80_maxDP110_MAF01_HWE001 --recode --aec --out DP5_MD80_loci 
/prj/mar-in-gen/bin/plink/plink --bfile Leop_NoDups_Biall_GQDP5_MD90_maxDP127_MAF01_HWE001 --recode --aec --out DP5_MD90_loci
/prj/mar-in-gen/bin/plink/plink --bfile Leop_NoDups_Biall_GQDP10_MD80_maxDP129_MAF01_HWE001 --recode --aec --out DP10_MD80_loci 
/prj/mar-in-gen/bin/plink/plink --bfile Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01_HWE001 --recode --aec --out DP10_MD90_loci
cut -d '      ' -f1,4 DP5_MD80_loci.map > DP5_MD80_loci.txt
cut -d '      ' -f1,4 DP5_MD90_loci.map > DP5_MD90_loci.txt
cut -d '      ' -f1,4 DP10_MD80_loci.map > DP10_MD80_loci.txt
cut -d '      ' -f1,4 DP10_MD90_loci.map > DP10_MD90_loci.txt
vcftools --vcf Leop_AllDups.recode.vcf --positions DP5_MD80_loci.txt --recode --out DP5_MD80_dups
vcftools --vcf Leop_AllDups.recode.vcf --positions DP5_MD90_loci.txt --recode --out DP5_MD90_dups
vcftools --vcf Leop_AllDups.recode.vcf --positions DP10_MD80_loci.txt --recode --out DP10_MD80_dups
vcftools --vcf Leop_AllDups.recode.vcf --positions DP10_MD90_loci.txt --recode --out DP10_MD90_dups
/prj/mar-in-gen/bin/plink/plink --vcf DP5_MD80_dups.recode.vcf --recodeA --double-id --aec --out DP5_MD80_dups
/prj/mar-in-gen/bin/plink/plink --vcf DP5_MD90_dups.recode.vcf --recodeA --double-id --aec --out DP5_MD90_dups
/prj/mar-in-gen/bin/plink/plink --vcf DP10_MD80_dups.recode.vcf --recodeA --double-id --aec --out DP10_MD80_dups
/prj/mar-in-gen/bin/plink/plink --vcf DP10_MD90_dups.recode.vcf --recodeA --double-id --aec --out DP10_MD90_dups
# DP5 - 80: 0.004972216 +/- 0.004
# DP5 - 90: 0.003637767 +/- 0.004
# DP10 - 80: 0.003635332 +/- 0.004
# DP10 - 90: 0.002858608 +/- 0.003 -> We'll work with this!

# As a final step let's check missing data per sample
/grp/animalbehaviour/davidlee/bin/plink/plink --bfile Leop_NoDups_Biall_GQDP10_MD90_maxDP149_MAF01_HWE001 --missing --aec --out missD
# Remove 7 samples with > 20% missing data
129_1_sorted
137_1_sorted
78128_1_sorted
93140_1_sorted
93143_1_sorted
93145_1_sorted
96_1_sorted
