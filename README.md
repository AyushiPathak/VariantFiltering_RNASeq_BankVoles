# VariantFiltering_RNASeq_BankVoles
Variant calling and filtering analysis for RNA-seq data of Bank Voles. This pipeline is a part of the Allele Specific Expression pipeline. | LU Master’s Thesis in Bioinformatics 2023 | BINP51

# Project Description

- what pipeline/ project does
- Technologies Used
- Challenges faces and featured you want to implement in future

Table of contents

Installation and Running the project

How to use the Project

Credits

Badges

Include Tests

# Uppmax Set-up

## Sequencing stats

```bash
for file in /proj/uppstore2017199/b2016119_nobackup/ayushi/1_merged_files/*fastq.gz ; do name=$( echo $file | cut -d '/' -f 7 | cut -d '.' -f 1); echo $name; gzip $file | wc -l ; done
for file in /proj/uppstore2017199/b2016119_nobackup/ayushi/1_merged_files/*fastq ; do cat $file | grep '@'| wc -l ; done >  /proj/uppstore2017199/b2016119_nobackup/ayushi/1_merged_files/seqStat.txt
```

# Quality Check Trimming

## Pre-trimming fastqc

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 1-00:00:00
#SBATCH -J fastqc.job
#SBATCH -o fastqc.out
#SBATCH -e fastqc.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools FastQC/0.11.9

raw_read_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/1_merged_files"
out_dir_fastqc="/proj/uppstore2017199/b2016119_nobackup/ayushi/3_fastqc"

#fastqc of untrimmed reads and copying the reads
fastqc -t 4 $raw_read_dir/*_R1.fastq.gz -o $out_dir_fastqc/R1.untrimmedfastqc
fastqc -t 4 $raw_read_dir/*_R2.fastq.gz -o $out_dir_fastqc/R2.untrimmedfastqc
```

## Trimming

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J trimmimg.job
#SBATCH -o trimmimg.out
#SBATCH -e trimmimg.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools trimmomatic/0.39

raw_read_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/1_merged_files"
out_dir_trimmed="/proj/uppstore2017199/b2016119_nobackup/ayushi/2_trimmed_reads"

cd $raw_read_dir

files=$('*fastq.gz')

for file in $files
do
cp $file $SNIC_TMP
done

cd $SNIC_TMP

prefixes=$(ls *.fastq.gz  | cut -d "_" -f 1-2) #P6207_112

for prefix in $prefixes
do
trimmomatic PE -threads 20 -basein ${prefix}_R1.fastq.gz -baseout ${prefix}_trimmed.fastq.gz ILLUMINACLIP:$TRIMMOMATIC_HOME/adapters/TruSeq2-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25
cp ${prefix}_trimmed_[12]P.fastq.gz $out_dir_trimmed
done
```

## Post-trimming quality check

```bash

```

# Alignment

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 10-00:00:00
#SBATCH -J mapping.job
#SBATCH -o mapping.out
#SBATCH -e mapping.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#load the modules
module load bioinfo-tools star/2.7.9a

#Provide trimmed read dir
trimmed_read_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/2_trimmed_reads"

#Provide the directory of the genome
genome_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/resources"

#Provide the name of the genome fasta file
genome="GCF_902806735.1_Bank_vole1_10x_genomic.fna"

#Provide GTF 3_fasta_files
GTF_file='genomic.gtf'

#Specify an output directory - need to create it beforehand
output_dir_indexing="/proj/uppstore2017199/b2016119_nobackup/ayushi/5_indexing"
output_dir_mapping="/proj/uppstore2017199/b2016119_nobackup/ayushi/6_mapping"

#Copy trimmed reads to the working directory
cd $trimmed_read_dir
cp *P*.fastq $SNIC_TMP

#Copy the genome assembly to the working directory
cd $genome_dir
cp $genome $SNIC_TMP
cp $GTF_file $SNIC_TMP

#Go to the working directory
cd $SNIC_TMP

#Loop through the samples and map them separately
samples=$(ls *.fastq | cut -d '.' -f 1 | sort | uniq)

for sample in $samples
do
STAR --runThreadN 20 --genomeDir $output_dir_indexing  --readFilesIn ${sample}.P1.fastq ${sample}.P2.fastq  --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix $output_dir_mapping/${sample}_mapped --outSAMtype BAM SortedByCoordinate
done
```

## Alignment Metrics

```bash
cd /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling
mkdir sizeMetrics
mkdir alignmentMetrics
```

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-06:00:00
#SBATCH -J 9.2.1_gatk_CollectAlignmentSummaryMetrics.job
#SBATCH -o 9.2.1_gatk_CollectAlignmentSummaryMetrics.out
#SBATCH -e 9.2.1_gatk_CollectAlignmentSummaryMetrics.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#load the modules
module load bioinfo-tools
module load  GATK/4.3.0.0

#Provide the directory of the input files
input_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/8.2.1_sorted_bam_bamaddrg/sorting_chr_pos"

#Provide the directory of the input files
output_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling"

#reference genome
ref_genome="/proj/uppstore2017199/b2016119_nobackup/ayushi/resources"

#Copy input files to the working directory
cd $input_dir
cp *.bam $SNIC_TMP
cp *.bai $SNIC_TMP

#Go to wirking directory
cd $SNIC_TMP

#Loop through the samples to feed them in a loop
samples=$(ls *.sorted.bam | cut -d '.' -f 1)

for sample in $samples
do
gatk CollectAlignmentSummaryMetrics -R $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna -I ${sample}.sorted.bam  -O $output_dir/alignmentMetrics/${sample}.alignemntMetrics.txt
gatk CollectInsertSizeMetrics  -R $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna -I ${sample}.sorted.bam  -O $output_dir/alignmentMetrics/${sample}.sizeMetrics.txt -H $output_dir/alignmentMetrics/${sample}.sizeMetrics.hist.pdf
done
```

```bash
module load bioinfo-tools
module load bamtools/2.5.2
module load samtools/1.9
module load picard/2.23.4
files=$(ls /proj/uppstore2017199/b2016119_nobackup/ayushi/6_mapping/*mappedAligned.sortedByCoord.out.bam); for file in $files; do name=$(echo $file | cut -d '/' -f 7| cut -d '_' -f 1-2); echo $name; samtools flagstat $file ; done > /proj/uppstore2017199/b2016119_nobackup/ayushi/6_mapping/Falgstat.txt
```

# Removing Duplicates

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 00-06:00:00
#SBATCH -J removing_dupes.job
#SBATCH -o removing_dupes.out
#SBATCH -e removing_dupes.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#load the modules
module load bioinfo-tools
module load picard/2.27.5

#Provide the directory of the input files
input_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/6_mapping"

#Provide the directory of the input files
output_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/7_removing_duplicates"

#Copy input files to the working directory
cd $input_dir
cp *.sortedByCoord.out.bam $SNIC_TMP

#Go to wirking directory
cd $SNIC_TMP

#Loop through the samples to feed them in a loop
samples=$(ls *.sortedByCoord.out.bam | cut -d 'm' -f 1 | sort | uniq)

for sample in $samples
do
java -jar $PICARD_ROOT/picard.jar MarkDuplicates I=${sample}mappedAligned.sortedByCoord.out.bam O=${sample}.bam M=${sample}.txt
cp ${sample}.bam $output_dir
cp *.txt $output_dir
done
```

# Adding Read Groups

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 10:00:00
#SBATCH -J 8.2.2_bamaddrg.job
#SBATCH -o 8.2.2_bamaddrg.out
#SBATCH -e 8.2.2_bamaddrg.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# script for bamaddrg for left over samples
# sorting done on the basis of file names

#loading the modules
module load bioinfo-tools
module load samtools/1.9
echo "Modules Loaded!"

#activate the environment
conda activate bamaddrg-env
echo "Activated Conda environment!"
conda list

#Provide the directory of the input files
input_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/7_removing_duplicates"

#Provide the directory of the added read groups files
output_dir_rg="/proj/uppstore2017199/b2016119_nobackup/ayushi/8.2_adding_readgroup_bamaddrg"

#Provide the directory of the sorted files
output_dir_sorted='/proj/uppstore2017199/b2016119_nobackup/ayushi/8.2.1_sorted_bam_bamaddrg'

#Copy input files to the working directory
cd $input_dir
cp *.bam $SNIC_TMP

#Go to working directory
cd $SNIC_TMP

#Loop through the samples to feed them in a loop
samples=$(ls *.bam | cut -d '_' -f 1-2)

for sample in $samples
do
bamaddrg -b ${sample}_.bam > $output_dir_rg/${sample}.addrg.bam
done

conda deactivate
```

# Sorting and Indexing

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 10:00:00
#SBATCH -J 8.2.4_indexing.job
#SBATCH -o 8.2.4_indexing.out
#SBATCH -e 8.2.4_indexing.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# first indexing script
# files are sorted --> chr position

#load the modules
module load bioinfo-tools
module load samtools/1.9

#Provide the directory for the input files
input_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/8.2.1_sorted_bam_bamaddrg/sorting_chr_pos"

#Copy input files to the working directory
cd $input_dir
cp *.bam $SNIC_TMP

#Go to working directory
cd $SNIC_TMP

#Loop through the samples to feed them in a loop
samples=$(ls *.bam | cut -d '.' -f 1-2)

for sample in $samples
do
samtools index ${sample}.bam
cp *bai $input_dir
done
```

# Variant Calling and Filtering

## Freebayes

### Variant Calling

```python
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-15:00:00
#SBATCH -J 9.1_freebayes_VC.job
#SBATCH -o 9.1_freebayes_VC.out
#SBATCH -e 9.1_freebayes_VC.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# First script for freebayes

# Load the modules
module load bioinfo-tools
module load freebayes/1.3.2
module load htslib/1.17
module load gnuparallel

# Provide the directory of the input files
input_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/8.2.1_sorted_bam_bamaddrg/sorting_chr_pos"

# Provide the directory of the input files
output_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling"

# Reference genome directory
ref_genome="/proj/uppstore2017199/b2016119_nobackup/ayushi/resources"

# Copy input files to the working directory
cd $input_dir
cp *.sorted.bam $SNIC_TMP
cp *bai $SNIC_TMP

# Go to wirking directory
cd $SNIC_TMP

freebayes-parallel <(fasta_generate_regions.py $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna.fai 100000) 15 -f $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna --genotype-qualities -L $output_dir/listOfBamFiles.txt > $output_dir/freebayes_output_all.vcf
bgzip -f $output_dir/freebayes_output_all.vcf.gz

# Indexing vcf file
tabix -p vcf $output_dir/freebayes_output_all.vcf.gz
```

### Statistical Analysis

```bash
cd /proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling
mkdir StatAnalysis
```

```bash
bcftools stats /proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/freebayes_output_all.vcf.gz > /proj/uppstore2017199/b2016119_noback
up/ayushi//9.1_freebays_variant_calling/queryingStat.txt
cat /proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/queryingStat.txt | grep -E '^SN'
```

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-1:00:00
#SBATCH -J 10_freebayesVcf_StatAnalysis.job
#SBATCH -o 10_freebayesVcf_StatAnalysis.out
#SBATCH -e 10_freebayesVcf_StatAnalysis.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#load the modules
module load bioinfo-tools
module load vcftools/0.1.16

input_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/freebayes_output_all.vcf.gz'
output='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/StatAnalysis/freebayesVcf_StatAnalysis'

vcftools --gzvcf ${input_dir} --freq2 --out ${output} --max-alleles 2
vcftools --gzvcf ${input_dir} --depth --out ${output}
vcftools --gzvcf ${input_dir} --site-mean-depth --out ${output}
vcftools --gzvcf ${input_dir} --site-quality --out ${output}
vcftools --gzvcf ${input_dir} --missing-indv --out ${output}
vcftools --gzvcf ${input_dir} --missing-site --out ${output}
vcftools --gzvcf ${input_dir} --het --out ${output}
```

```bash
awk '$6 <1e-5' hwe_output.hwe > Variants_to_filter_hwe.txt
cat Variants_to_filter_hwe.txt  | wc -l
#444
```

### Visualizing with R

Copying the files from server to local

```bash
cd Downloads/BINP51/
scp ayuship@rackham.uppmax.uu.se://proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/VCF_Stat* .
```

Preparing the R Studio to work in right folder

```r
getwd()
setwd("C:/Users/Ayushi Pathak/Downloads")
```

```r
# Install and load ggplot2 if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

# Install and load gridExtra if not already installed
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

library(gridExtra)

# import data sets
VCF_StatAnalysis_lmiss <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_lmiss", row.names=NULL)
VCF_StatAnalysis_lqual <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_lqual", row.names=NULL)
VCF_StatAnalysis_het <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_het", row.names=NULL)
VCF_StatAnalysis_frq <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_frq", row.names=NULL)
VCF_StatAnalysis_idepth <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_idepth")
VCF_StatAnalysis_imiss <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_imiss")
VCF_StatAnalysis_ldepth_mean <- read.delim("C:/Users/Ayushi Pathak/Downloads/VCF_StatAnalysis_ldepth_mean")
VCF_StatAnalysis_hwe <- read.delim("C:/Users/Ayushi Pathak/Downloads/hwe_output.hwe")

## Variant based statistics

# Variant quality
new_header_lqual<-c('CHROM','POS','QUAL')
colnames(VCF_StatAnalysis_lqual)<-new_header_lqual
head(VCF_StatAnalysis_lqual,n=10)

var_qual<-ggplot(VCF_StatAnalysis_lqual,aes(QUAL))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "QUAL", y = "DENSITY")

summary(VCF_StatAnalysis_lqual$QUAL)

var_qual_1000<-var_qual+xlim(0,1000) + theme_light()
var_qual_100<-var_qual+xlim(0,100) + theme_light()

#Variant mean depth

new_header_meandepth<-c("chr", "pos", "mean_depth", "var_depth")
colnames(VCF_StatAnalysis_ldepth_mean)<-new_header_meandepth
head(VCF_StatAnalysis_ldepth_mean,n=10)

mean_depth<-ggplot(VCF_StatAnalysis_ldepth_mean,aes(mean_depth))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "MEAN_DEPTH", y = "DENSITY")

summary(VCF_StatAnalysis_ldepth_mean$mean_depth)

mean_depth_100<-mean_depth+ theme_light() + xlim(0, 100)
mean_depth_50<-mean_depth+ theme_light() + xlim(0, 50)

# Variant missingness

new_header_lmiss<- c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss")
colnames(VCF_StatAnalysis_lmiss)<-new_header_lmiss
head(VCF_StatAnalysis_lmiss,n=10)

var_miss<-ggplot(VCF_StatAnalysis_lmiss,aes(fmiss))+geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "MISSING FREQUENCY", y = "DENSITY")

summary(VCF_StatAnalysis_lmiss$fmiss)

# Minor allele frequency

new_header<-c('CHROM','POS','ALLELE','CHR','FREQ1','FREQ2')
colnames(VCF_StatAnalysis_frq)<-new_header
head(VCF_StatAnalysis_frq,n=10)
VCF_StatAnalysis_frq$maf<-apply(VCF_StatAnalysis_frq[,c('FREQ1','FREQ2')],1,function(row) min(row))

var_frq<- ggplot(VCF_StatAnalysis_frq,aes(maf)) +geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "MAF", y = "DENSITY")
head(VCF_StatAnalysis_frq)
summary(VCF_StatAnalysis_frq$maf)

#HWE

head(VCF_StatAnalysis_hwe,n=10)
var_hwe<- ggplot(VCF_StatAnalysis_hwe,aes(P_HWE)) +geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "HWE", y = "DENSITY")
head(VCF_StatAnalysis_hwe)
summary(VCF_StatAnalysis_hwe$P_HWE)

#hwe_100<-var_hwe +  xlim(0,0.00001) # to look in the amount of variants we will filter
#hwe_100

## Individual based statistics

# Mean depth per individual

new_header_idepth<-c("ind", "nsites", "depth")
colnames(VCF_StatAnalysis_idepth)<-new_header_idepth
head(VCF_StatAnalysis_idepth,n=19)

ind_depth<-ggplot(VCF_StatAnalysis_idepth,aes(depth))+geom_histogram(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "DEPTH", y = "COUNT")

# Proportion of missing data per individual

new_header_imiss<-c("ind", "ndata", "nfiltered", "nmiss", "imiss")
colnames(VCF_StatAnalysis_imiss)<-new_header_imiss

head(VCF_StatAnalysis_imiss,n=19)

ind_miss<-ggplot(VCF_StatAnalysis_imiss,aes(imiss))+geom_histogram(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+labs(x = "MISSING DATA / INDIVIDUAL", y = "COUNT")

# Heterozygosity and inbreeding coefficient per individual

new_header_het<-c("ind","ho", "he", "nsites", "f")
colnames(VCF_StatAnalysis_het)<-new_header_het
head(VCF_StatAnalysis_het,n=10)

ind_het<-ggplot(VCF_StatAnalysis_het,aes(f)) +geom_histogram(fill = "#1F4E79", colour = "black", alpha = 0.3) + theme_light()+ theme_light()+labs(x = "INBREEDING COEFFICIENT", y = "COUNT")

##arranging the graphs
combine_plot_f<- var_qual_100 / mean_depth_50 / var_miss/ var_frq/ ind_miss / var_hwe  +plot_layout(ncol = 2, nrow = 3)
combine_plot_f

#var_qual_100 + mean_depth_50 + var_miss + var_frq  + ind_miss + var_hwe
```

### Filtering SNPs and biallelic genes

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-1:00:00
#SBATCH -J 11_snp_only_biallele_freebayes.job
#SBATCH -o 11_snp_only_biallele_freebayes.out
#SBATCH -e 11_snp_only_biallele_freebayes.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL
# to remove multi allelic variants
# to remove indels

# loading the modules
module load bioinfo-tools
module load vcftools/0.1.16

# working directory
dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling'

vcftools --gzvcf ${dir}/freebayes_output_all.vcf.gz --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out ${dir}/snps_only_biallele_freebayes
```

### VCF filtering

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 5
#SBATCH -t 00-1:00:00
#SBATCH -J 12.1_vcfFilteringFreebayes.job
#SBATCH -o 12.1_vcfFilteringFreebayes.out
#SBATCH -e 12.1_vcfFilteringFreebayes.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# to filter the variants using decided parameters

# loading the modules
module load bioinfo-tools
module load vcftools/0.1.16
module load vcflib
module load bcftools/1.9
module load htslib/1.8
module load python/3.8.7

# input directory
in_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling'

# output directory
O='/proj/uppstore2017199/b2016119_nobackup/ayushi/10_variant_filtering_freebayes'

# setting filters
maf=0.05
max_miss=0.25
min_depth=8
hwe=0.00001

# copy files in tmp folder
cd $in_dir
cp snps_only_biallele_freebayes.recode.vcf $SNIC_TMP

#go to working directory
cd $SNIC_TMP

#get the file and prefix
file='snps_only_biallele_freebayes.recode.vcf'
sample=$(ls *freebayes.recode.vcf | cut -d '_' -f 4 | cut -d '.' -f 1)

vcffilter -f 'SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & QUAL > 30' ${file} > $O/${sample}_vcffiltered.vcf
vcftools --vcf $O/${sample}_vcffiltered.vcf --maf $maf --max-missing $max_miss --minDP $min_depth --hwe $hwe --recode --recode-INFO-all --stdout | gzip -c > $O/${sample}_vcftoolsfiltered.vcf.gz
vcftools --gzvcf $O/${sample}_vcftoolsfiltered.vcf.gz --exclude-positions  /proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/Variants_to_filter_hwe.txt --recode --out  $O/${sample}_hwe_filtered
vcfallelicprimitives --keep-geno --keep-info $O/${sample}_hwe_filtered.recode.vcf > $O/${sample}_vcfAllelicPrimitivesfiltered.vcf
bgzip -c $O/${sample}_vcfAllelicPrimitivesfiltered.vcf > $O/${sample}_filtered.vcf.gz
bcftools sort -Oz $O/${sample}_filtered.vcf.gz -o $O/${sample}_filtered_sorted.vcf.gz
tabix -f -p vcf $O/${sample}_filtered_sorted.vcf.gz
```

Run 2

```bash
# loading the modules
module load bioinfo-tools
module load vcftools/0.1.16
module load vcflib
module load bcftools/1.9
module load htslib/1.8
module load python/3.8.7

# input directory
in_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling'

# output directory
O='/proj/uppstore2017199/b2016119_nobackup/ayushi/10_v2_variant_filtering_freebayes'

# setting filters
maf=0.05
max_miss=0.5
min_depth=8
hwe=0.0001

# copy files in tmp folder
cd $in_dir
cp snps_only_biallele_freebayes.recode.vcf $SNIC_TMP

#go to working directory
cd $SNIC_TMP

#get the file and prefix
file='snps_only_biallele_freebayes.recode.vcf'
sample='v2freebayes'

vcffilter -f 'SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & QUAL > 40' ${file} > $O/${sample}_vcffiltered.vcf
vcftools --vcf $O/${sample}_vcffiltered.vcf --maf $maf --max-missing $max_miss --minDP $min_depth --hwe $hwe --recode --recode-INFO-all --stdout | gzip -c > $O/${sample}_vcftoolsfiltered.vcf.gz
vcftools --gzvcf $O/${sample}_vcftoolsfiltered.vcf.gz --exclude-positions  /proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/Variants_to_filter_hwe.txt --recode --out  $O/${sample}_hwe_filtered
vcfallelicprimitives --keep-geno --keep-info $O/${sample}_hwe_filtered.recode.vcf > $O/${sample}_vcfAllelicPrimitivesfiltered.vcf
bgzip -c $O/${sample}_vcfAllelicPrimitivesfiltered.vcf > $O/${sample}_filtered.vcf.gz
bcftools sort -Oz $O/${sample}_filtered.vcf.gz -o $O/${sample}_filtered_sorted.vcf.gz
tabix -f -p vcf $O/${sample}_filtered_sorted.vcf.gz
```

Run3

```bash
# loading the modules
module load bioinfo-tools
module load vcftools/0.1.16
module load vcflib
module load bcftools/1.9
module load htslib/1.8
module load python/3.8.7

# input directory
in_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling'

# output directory
O='/proj/uppstore2017199/b2016119_nobackup/ayushi/10_v3_variant_filtering_freebayes'

# setting filters
maf=0.05
max_miss=0.75
min_depth=15
hwe=0.0001

# copy files in tmp folder
cd $in_dir
cp snps_only_biallele_freebayes.recode.vcf $SNIC_TMP

#go to working directory
cd $SNIC_TMP

#get the file and prefix
file='snps_only_biallele_freebayes.recode.vcf'
sample='v3freebayes'

vcffilter -f 'SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1 & QUAL > 60' ${file} > $O/${sample}_vcffiltered.vcf
vcftools --vcf $O/${sample}_vcffiltered.vcf --maf $maf --max-missing $max_miss --minDP $min_depth --hwe $hwe --recode --recode-INFO-all --stdout | gzip -c > $O/${sample}_vcftoolsfiltered.vcf.gz
vcftools --gzvcf $O/${sample}_vcftoolsfiltered.vcf.gz --exclude-positions  /proj/uppstore2017199/b2016119_nobackup/ayushi/9.1_freebays_variant_calling/Variants_to_filter_hwe.txt --recode --out  $O/${sample}_hwe_filtered
vcfallelicprimitives --keep-geno --keep-info $O/${sample}_hwe_filtered.recode.vcf > $O/${sample}_vcfAllelicPrimitivesfiltered.vcf
bgzip -c $O/${sample}_vcfAllelicPrimitivesfiltered.vcf > $O/${sample}_filtered.vcf.gz
bcftools sort -Oz $O/${sample}_filtered.vcf.gz -o $O/${sample}_filtered_sorted.vcf.gz
tabix -f -p vcf $O/${sample}_filtered_sorted.vcf.gz
```

## GATK

### Split N Cigar reads

```bash
mkdir calledVariants SplitNCigarReads
```

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-20:00:00
#SBATCH -J 9.2.2_GATK_splitNCigars_haplotypeCaller.job
#SBATCH -o 9.2.2_GATK_splitNCigars_haplotypeCaller.out
#SBATCH -e 9.2.2_GATK_splitNCigars_haplotypeCaller.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# model script with 101

#load the modules
module load bioinfo-tools
module load GATK/4.3.0.0

#Provide the directory of the input files
input_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/8.2.1_sorted_bam_bamaddrg/sorting_chr_pos"

#Provide the directory of the input files
output_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling//SplitNCigarReads"
output_dir2="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants"

#reference genome
ref_genome="/proj/uppstore2017199/b2016119_nobackup/ayushi/resources"

#Copy input files to the working directory
cd $input_dir
cp P6207_101.sorted.bam $SNIC_TMP
cp P6207_101.sorted.bam.bai $SNIC_TMP

#Go to working directory
cd $SNIC_TMP

#Loop through the samples to feed them in a loop
samples=$(ls *.sorted.bam | cut -d '.' -f 1)

for sample in $samples
do
gatk SplitNCigarReads -R $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna -I ${sample}.sorted.bam -O $output_dir/${sample}.bam
gatk HaplotypeCaller -R $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna -I $output_dir/${sample}.bam  -O $output_dir2/${sample}.rawVariants.vcf -ERC GVCF
done
```

### Haplotype caller

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-23:40:00
#SBATCH -J 9.2.2_gatk_haplotypeCaller_rerun.job
#SBATCH -o 9.2.2_gatk_haplotypeCaller_rerun.out
#SBATCH -e 9.2.2_gatk_haplotypeCaller_rerun.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# gatk haplotype caller after Split N Cigar Reads

#load the modules
module load bioinfo-tools
module load GATK/4.3.0.0

#Provide the directory of the input files
output_dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling//SplitNCigarReads"
output_dir2="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants"

#reference genome
ref_genome="/proj/uppstore2017199/b2016119_nobackup/ayushi/resources"

#Copy input files to the working directory
cd $output_dir
cp *bam $SNIC_TMP
cp *bam.bai $SNIC_TMP

#Go to working directory
cd $SNIC_TMP

#Loop through the samples to feed them in a loop
samples=$(ls *bam | cut -d '.' -f 1)

for sample in $samples
do
gatk HaplotypeCaller -R $ref_genome/GCF_902806735.1_Bank_vole1_10x_genomic.fna -I $output_dir/${sample}.bam  -O $output_dir2/${sample}.rawVariants.vcf -ERC GVCF
done
```

### Merging the VCF files

```bash
#run in interactive
cd /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants
module load htslib/1.9
for file in *.vcf ; do bgzip $file ; done
for file in *.vcf.gz ; do tabix -p vcf ${file} -f ; done

```

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-05:00:00
#SBATCH -J 9.2.3_mergingFiles_bcftools.job
#SBATCH -o 9.2.3_mergingFiles_bcftools.out
#SBATCH -e 9.2.3_mergingFiles_bcftools.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#load the modules
module load bioinfo-tools
module load bcftools/1.8
module load python/3.8.7

#Provide the directory of the input and output files
dir="/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants"

#Copy input files to the working directory
cd $dir
cp *.rawVariants.vcf.gz $SNIC_TMP
cp *.rawVariants.vcf.gz.tbi $SNIC_TMP

cd $SNIC_TMP

bcftools merge \
	P6207_101.rawVariants.vcf.gz \
	P6207_104.rawVariants.vcf.gz \
	P6207_105.rawVariants.vcf.gz \
	P6207_109.rawVariants.vcf.gz \
	P6207_112.rawVariants.vcf.gz \
	P6207_113.rawVariants.vcf.gz \
	P6207_116.rawVariants.vcf.gz \
	P6207_117.rawVariants.vcf.gz \
	P6207_120.rawVariants.vcf.gz \
	P6207_121.rawVariants.vcf.gz \
	P6207_123.rawVariants.vcf.gz \
	P6207_125.rawVariants.vcf.gz \
	P6207_128.rawVariants.vcf.gz \
	P6207_129.rawVariants.vcf.gz \
	P6207_132.rawVariants.vcf.gz \
	P6207_134.rawVariants.vcf.gz \
	P6207_135.rawVariants.vcf.gz \
	P6207_136.rawVariants.vcf.gz \
	-o $dir/mergedVCF_bcftools.vcf
```

### Statistical analysis

```bash
cd /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling
mkdir StatAnalysis
```

```bash
grep -v '#' biallele_snp_gatk.vcf | wc -l
#2199466
```

```bash
bcftools query -l  /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/DBImport/genotypeOutput.vcf.gz
```

```bash
bcftools stats /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/DBImport/genotypeOutput.vcf.gz > /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/DBImport/queryingStat.txt
bcftools stats genotypeOutput.vcf.gz | grep -E '^SN'
```

```bash
bcftools query biallele_snp_gatk.vcf -f '%FS\t%SOR\t%ReadPosRankSum\t%QUAL\t%AF\t%DP\t%InbreedingCoeff\n' > snps.metrics.txt
echo -e "FS\tSOR\tReadPosRankSum\tQUAL\tAF\tDP\tInbreedingCoeff" > header
cat header snps.metrics.txt > snps.metrics.tsv
```

```bash
scp ayuship@rackham.uppmax.uu.se:/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/snps.metrics.tsv .
```

```bash
# Install and load ggplot2 if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

# Install and load gridExtra if not already installed
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra")
}

library(gridExtra)

# import data sets
#snps.metrics<- read.delim("C:/Users/Ayushi Pathak/Downloads/snps.metrics.tsv", row.names=NULL)
snps.metrics <- read.delim("C:/Users/Ayushi Pathak/Downloads/snps.metrics.tsv", na.strings=".")

# Fisher Strand
FS<-ggplot(snps.metrics,aes(FS))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+ labs(y='DENSITY')
#FS100<-FS+xlim(0,100) + theme_light()
#FS50<-FS+xlim(0,50) + theme_light()
FS25<-FS+xlim(0,25) + theme_light()
FS10<-FS+xlim(0,10) + theme_light()
summary(snps.metrics$FS)

#FS100
#FS50
#FS10

# Strand Odds Ratio
SOR<-ggplot(snps.metrics,aes(SOR))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+ labs(y='DENSITY')
SOR10<-SOR+xlim(0,10) + theme_light()
summary(snps.metrics$SOR)

# Read Position Rank Sum Test
ReadPosRankSum<-ggplot(snps.metrics,aes(ReadPosRankSum))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+ labs(y='DENSITY')
ReadPosRankSum100<-ReadPosRankSum+xlim(-10,10) + theme_light()
summary(snps.metrics$ReadPosRankSum)

# Quality
QUAL<-ggplot(snps.metrics,aes(QUAL))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+ labs(y='DENSITY')
QUAL100<- QUAL + xlim(0,150)
summary(snps.metrics$QUAL)

# Allele Frequency
snps.metrics$AF<-as.numeric(snps.metrics$AF)
AF<- ggplot(snps.metrics,aes(x=AF))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light() + labs(y='DENSITY')
AF1<-AF +xlim(0,1) + theme_light()
summary(snps.metrics$AF)

# Depth
DP <- ggplot(snps.metrics,aes(DP))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+ labs(y='DENSITY')
DP100<- DP+xlim(0,100)
summary(snps.metrics$DP)

# HWE
HWE<- ggplot(snps.metrics,aes(InbreedingCoeff))+ geom_density(fill = "#29ADB2", colour = "black", alpha = 0.3)+ theme_light()+ labs(y='DENSITY')

library(patchwork)

combine_plot <- FS10 / SOR10 /ReadPosRankSum100 /QUAL100 /AF1/ DP100 / HWE +plot_layout(ncol = 2, nrow = 4)
combine_plot
```

### Joint Genotyping

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-7:00:00
#SBATCH -J 11.2_jointGenotyping.job
#SBATCH -o 11.2_jointGenotyping.out
#SBATCH -e 11.2_jointGenotyping.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# loading the modules
module load bioinfo-tools
module load GATK/4.3.0.0

# files variables
ref='/proj/uppstore2017199/b2016119_nobackup/ayushi/resources/GCF_902806735.1_Bank_vole1_10x_genomic.fna '
in_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants/mergedVCF_bcftools.vcf'
out_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling'

gatk --java-options "-Xmx4g" GenotypeGVCFs -R $ref -V $in_dir -O $out_dir/genotypedPhased.vcf.gz

```

### Filtering SNPs and biallelic genes

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 00-1:00:00
#SBATCH -J 11.2_snp_only_biallele_gatk.job
#SBATCH -o 11.2_snp_only_biallele_gatk.out
#SBATCH -e 11.2_snp_only_biallele_gatk.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# to remove multi allelic variants
# to remove indels

# loading the modules
module load bioinfo-tools
module load GATK/4.3.0.0

# files variables
ref='/proj/uppstore2017199/b2016119_nobackup/ayushi/resources'
in_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/DBImport'
out_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling'

cd in_dir
cp genotypeOutput* $SNIC_TMP

cd $ref
cp GCF_902806735.1_Bank_vole1_10x_genomic.fna  $SNIC_TMP
cp GCF_902806735.1_Bank_vole1_10x_genomic.dict $SNIC_TMP

cd $SNIC_TMP

gatk SelectVariants -R ${ref}/GCF_902806735.1_Bank_vole1_10x_genomic.fna -V ${in_dir}/genotypeOutput.vcf.gz --restrict-alleles-to BIALLELIC -O ${out_dir}/biallele_gatk.vcf
gatk SelectVariants -R ${ref}/GCF_902806735.1_Bank_vole1_10x_genomic.fna -V ${out_dir}/biallele_gatk.vcf --select-type SNP -O ${out_dir}/biallele_snp_gatk.vcf
```

### VCF filtering

```bash
# add the inbreeding coefficient in GATK parameters
gatk CalculateGenotypePosteriors -R /proj/uppstore2017199/b2016119_nobackup/ayushi/resources/GCF_902806735.1_Bank_vole1_10x_genomic.fna -V biallele_snp_gatk.vcf -O posteriors.vcf
# creating a table to compare Inbreeding Coefficient values
#example
gatk VariantsToTable -V posteriors.vcf -F CHROM -F POS -F ID -GF GT -GF GQ -GF AD -GF PL -GF PS -O posteriors.txt
gatk VariantsToTable -V posteriors.vcf -F CHROM -F POS -F ID -F InbreedingCoeff -O InbreedingCoeff.txt
```

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 5
#SBATCH -t 00-1:00:00
#SBATCH -J 12.1_vcfFilteringGATK.job
#SBATCH -o 12.1_vcfFilteringGATK.out
#SBATCH -e 12.1_vcfFilteringGATK.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

# to filter the variants using decided parameters

# loading the modules
module load bioinfo-tools
module load GATK/4.3.0.0

# input directory
in='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/posteriors.vcf'

# output directory
O='/proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK'

gatk VariantFiltration -V $in\
	-filter "FS > 5.0 " --filter-name 'FS5' \
	-filter "SOR > 3.0" --filter-name 'SOR3' \
	-filter "ReadPosRankSum < -3.0" --filter-name 'RPRSum-3' \
	-filter "ReadPosRankSum > 3.0" --filter-name 'RPRSum3' \
	-filter "QUAL < 50.0" --filter-name 'QUAL50' \
  -filter "AF < 0.05" --filter-name 'AF0.05' \
  -filter "DP < 8" --filter-name 'DP8' \
  -filter "InbreedingCoeff < 0.00001" --filter-name 'InbreedCo0.00001' \
  -O $O/GATK_filteredRun1.vcf
```

Run 2

```bash
gatk VariantFiltration -V $in\
	-filter "FS > 1.5 " --filter-name 'FS1.5' \
	-filter "SOR > 1.5" --filter-name 'SOR1.5' \
	-filter "ReadPosRankSum < -1.5" --filter-name 'RPRSum-1.5' \
	-filter "ReadPosRankSum > 1.5" --filter-name 'RPRSum1.5' \
	-filter "QUAL < 70.0" --filter-name 'QUAL70' \
  -filter "AF < 0.05" --filter-name 'AF0.05' \
  -filter "DP < 8" --filter-name 'DP8' \
  -filter "InbreedingCoeff < 0.0001" --filter-name 'InbreedCo0.0001' \
  -O $O/GATK_filteredRun2.vcf
```

Run 3

```bash

gatk VariantFiltration -V $in\
	-filter "FS > 0.75 " --filter-name 'FS0.75' \
	-filter "SOR > 1.0" --filter-name 'SOR1.0' \
	-filter "ReadPosRankSum < -1.0" --filter-name 'RPRSum-1' \
	-filter "ReadPosRankSum > 1.0" --filter-name 'RPRSum1' \
	-filter "QUAL < 85.0" --filter-name 'QUAL85' \
  -filter "AF < 0.05" --filter-name 'AF0.05' \
  -filter "DP < 15" --filter-name 'DP15' \
  -filter "InbreedingCoeff < 0.0001" --filter-name 'InbreedCo0.0001' \
  -O $O/GATK_filteredRun2.vcf
```

### Variants filtered by each filter

```bash
#to retreive the header tag from vcf file
grep '##FILTER' GATK_filteredRun1.vcf | cut -d '=' -f 3 | cut -d ',' -f 1 > filter_tab_Run1.txt
# to count the filtered variants
filters=$(grep '##FILTER' GATK_filteredRun1.vcf | cut -d '=' -f 3 | cut -d ',' -f 1); for filter in $filters ; do grep -v '^#'  GATK_filteredRun1.vcf | grep -c $filter; done > file01tab.txt
paste filter_tab_Run1.txt  file01tab.txt |  awk '{print $1"\t"$2}' > merged_file.txt
```

```bash
# Data
data<-data.frame(
		Filters=c('AF','DP','FS','InbreedCo','PASS','QUAL','RPRSum+','RPRSum-','SOR'),
    Run1=c(1040621,674393,145785,286649,718284,211516,2834,4195,70641),
    Run2=c(1040621,674393,387615,286649,471858,468683,45767,59069,675833),
    Run3=c(1040621,973680,453599,286649,150843,625252,232434,312028,1018710))

# Convert data to long format using tidyr
library(tidyr)
data_long <- gather(data, key = "Variable", value = "Run", -Filters)

library(ggplot2)
ggplot(data_long, aes(x = Filters, y = Run, fill = Variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Variants Filtered in Different Runs",
       x = "Filters",
       y = "Variants Filtered") +
  scale_fill_manual(values = c("Run1" = "#1F4E79", "Run2" = "#29ADB2", "Run3" = "#C5E898"))
```

# SNP Detection by RNA-Seq: Concordance with those Detected by DNA-Seq

### GFF File

```bash
cd /proj/uppstore2017199/b2016119_nobackup/ayushi
# create the CDS file
awk '$3 ~ /CDS/ ' genomic.gtf > CDS.gtf
# create the CDS+exon file
awk '$3 ~ /CDS/ ' genomic.gtf > CDS.exon.gtf
awk '$3 ~ /exon/' genomic.gtf >> CDS.exon.gtf
awk '$3 ~ /transcript/' genomic.gtf >> CDS.exon.gtf
```

### GTF File to BED File

```python
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $1, $4 - 1, $5, $9, ".", $7}' CDS.gtf > CDS.bed
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $1, $4 - 1, $5, $9, ".", $7}' CDS.exon.gtf > CDS.exon.bed
```

### Genomic VCF file

```python
bcftools view -R CDS.bed freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz | bcftools sort | bgzip -c >  genomicVariants.vcf.gz
bcftools view -R CDS.exon.bed freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz | bcftools sort | bgzip -c >CDS.exon.vcf.gz
```

### RNA VCF File - common variant file

```python
bcftools view -R  genomicVariants.vcf.gz 10_variant_filtering/freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c >  transcritomicVariants.vcf.gz
```

### File size comparison

```python
# **11G** Oct 21 08:42 freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz

# **130M** Oct 24 15:30 CDS.bed

# **108M** Oct 24 16:55 genomicVariants.vcf.gz
bcftools view -H genomicVariants.vcf.gz | wc
# **222035** 8881400 420201799

# **86M** Sep 20 15:28 freebayes_filtered_sorted.vcf.gz
bcftools view -H 10_variant_filtering/freebayes_filtered_sorted.vcf.gz | wc
# **223059** 6022593 288090764

#**23M** Oct 24 17:03 transcritomicVariants.vcf.gz
bcftools view -H transcritomicVariants.vcf.gz | wc
# **53846** 1453842 70759113
```

gft file—>CDS gtf file —> bed file —>   —> vcf file with variants in CDS region—>   —> common variants from RNA (transcriptomic) vcf file and DNA (genomic) vcf file in vcf format

genomic vcf file —————————>^

transcriptomic called variants vcf file ——————————————————-—> ^

GATK versus DNA variants

run 01

```bash
cat GATK_filteredRun1.vcf | grep -v '^#' | grep 'PASS' > GATK_filteredPass1.txt
grep '^#' GATK_filteredRun1.vcf > GATK_filteredPass1.header
cat GATK_filteredPass1.header GATK_filteredPass1.txt > GATK_filteredPass1.vcf
bgzip GATK_filteredPass1.vcf
tabix -f -p vcf GATK_filteredPass1.vcf.gz
bcftools view -H GATK_filteredPass1.vcf.gz | wc -l
# 718,284
```

```bash
cd ..
bcftools view -R  genomicVariants.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass1.vcf.gz | bcftools sort | bgzip -c > gatkVersusDnaVariants1.vcf.gz
bcftools view gatkVersusDnaVariants1.vcf.gz | wc -l
# 31,497
bcftools view -R  CDS.exon.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass1.vcf.gz | bcftools sort | bgzip -c > gatk.versus.cdsExon1.vcf.gz
bcftools view gatk.versus.cdsExon1.vcf.gz | wc -l
# 83934
bcftools view -R freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass1.vcf.gz | bcftools sort | bgzip -c > gatk.versus.all1.vcf.gz
bcftools view -H gatk.versus.all1.vcf.gz | wc -l
#
```

run 02

```bash
cat GATK_filteredRun2.vcf | grep -v '^#' | grep 'PASS' > GATK_filteredPass2.txt
grep '^#' GATK_filteredRun2.vcf > GATK_filteredPass2.header
cat GATK_filteredPass2.header GATK_filteredPass2.txt > GATK_filteredPass2.vcf
bgzip GATK_filteredPass2.vcf
tabix -f -p vcf GATK_filteredPass2.vcf.gz
bcftools view -H GATK_filteredPass2.vcf.gz | wc -l
# 471,858

```

```bash
cd ..
bcftools view -R  genomicVariants.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass2.vcf.gz | bcftools sort | bgzip -c > gatkVersusDnaVariants2.vcf.gz
bcftools view gatkVersusDnaVariants2.vcf.gz | wc -l
# 19,294
bcftools view -R  CDS.exon.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass2.vcf.gz | bcftools sort | bgzip -c > gatk.versus.cdsExon2.vcf.gz
bcftools view gatk.versus.cdsExon2.vcf.gz | wc -l
#50925
bcftools view -R freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass2.vcf.gz | bcftools sort | bgzip -c > gatk.versus.all2.vcf.gz
bcftools view gatk.versus.all2.vcf.gz | wc -l
# 249183
```

run 03

```bash
# file prep
cd  /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK
bgzip GATK_filteredRun3.vcf
tabix -f -p vcf GATK_filteredRun3.vcf.gz > GATK_filteredRun3.vcf.gz.tbi
bcftools view -H GATK_filteredRun3.vcf.gz | grep 'PASS' > GATK_filteredPass3.txt
bcftools view -h GATK_filteredRun3.vcf.gz  > GATK_filteredPass3.header
cat GATK_filteredPass3.header GATK_filteredPass3.txt > GATK_filteredPass3.vcf
bgzip GATK_filteredPass3.vcf
tabix -f -p vcf GATK_filteredPass3.vcf.gz
bcftools view -H GATK_filteredPass3.vcf.gz | wc -l
# 150,843
```

```bash
cd ..
bcftools view -R  genomicVariants.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass3.vcf.gz | bcftools sort | bgzip -c > gatkVersusDnaVariants3.vcf.gz
bcftools view gatkVersusDnaVariants3.vcf.gz | wc -l
# 6,003
bcftools view -R  CDS.exon.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass3.vcf.gz | bcftools sort | bgzip -c > gatk.versus.cdsExon3.vcf.gz
bcftools view gatk.versus.cdsExon3.vcf.gz | wc -l
# 16759
bcftools view -R freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK/GATK_filteredPass3.vcf.gz | bcftools sort | bgzip -c > gatk.versus.all3.vcf.gz
bcftools view gatk.versus.all3.vcf.gz | wc -l
# 77807
```

freebayes versus DNA run

```bash
# only CDS
bcftools view -R  genomicVariants.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_variant_filtering_freebayes_original/freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayesVaesusDnaVaraints1.vcf.gz
view -H  freebayesVaesusDnaVaraints1.vcf.gz | wc -l
# 53,801
bcftools view -R  genomicVariants.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v2*/v2freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayesVaesusDnaVaraints2.vcf.gz
view -H  freebayesVaesusDnaVaraints2.vcf.gz | wc -l
# 50,223
bcftools view -R  genomicVariants.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v3_variant_filtering_freebayes/v3freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayesVaesusDnaVaraints3.vcf.gz
bcftools view -H  freebayesVaesusDnaVaraints3.vcf.gz | wc -l
# 39,530

# CDS and Exons
bcftools view -R  CDS.exon.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_variant_filtering_freebayes_original/freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayes1.vs.cdsexon.vcf.gz
view -H  freebayes1.vs.cdsexon.vcf.gz | wc -l
#
bcftools view -R  CDS.exon.vcf.gz  /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v2*/v2freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayes2.vs.cdsexon.vcf.gz
view -H  freebayes2.vs.cdsexon.vcf.gz| wc -l
#
bcftools view -R  CDS.exon.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v3_variant_filtering_freebayes/v3freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayes3.vs.cdsexon.vcf.gz
bcftools view -H  freebayes3.vs.cdsexon.vcf.gz | wc -l
# 75048

# Whole Genome
bcftools view -R  freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_variant_filtering_freebayes_original/freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayes1.vs.allgenome.vcf.gz
view -H  freebayes1.vs.cdsexon.vcf.gz | wc -l
#
bcftools view -R  freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v2*/v2freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayes2.vs.allgenome.vcf.gz
view -H  freebayes2.vs.cdsexon.vcf.gz| wc -l
#
bcftools view -R  freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v3_variant_filtering_freebayes/v3freebayes_filtered_sorted.vcf.gz | bcftools sort | bgzip -c > freebayes3.vs.allgenome.vcf.gz
bcftools view -H  freebayes3.vs.cdsexon.vcf.gz | wc -l
#
```

VISUALISATION

```bash
#Prepping the genomic data base
bcftools view -H freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz | cut -f 1-2 |tr '\t' ':' > allGenomeChrPos.txt
tabix -p vcf CDS.exon.vcf.gz
bcftools view -H CDS.exon.vcf.gz |cut -f 1-2 |tr '\t' ':' > cdsExonChrPos.txt
bcftools view -H genomicVariants.vcf.gz | cut -f 1-2 |tr '\t' ':' > DNA_filtered1chrpos.txt

#Prepping the GATK files
cd /proj/uppstore2017199/b2016119_nobackup/ayushi/10.1_variant_filtering_GATK
bcftools view -H GATK_filteredPass1.vcf.gz| cut -f 1-2 |tr '\t' ':' > GATK_filteredPass1chrpos.txt
bcftools view -H GATK_filteredPass2.vcf.gz| cut -f 1-2 |tr '\t' ':' > GATK_filteredPa
ss2chrpos.txt
bcftools view -H GATK_filteredPass3.vcf.gz| cut -f 1-2 |tr '\t' ':' > GATK_filteredPa
ss3chrpos.txt

#Prepping the Freebayes file
bcftools view -H /proj/uppstore2017199/b2016119_nobackup/ayushi/10_variant_filtering_freebayes_original/freebayes_filtered_sorted.vcf.gz | cut -f 1-2 |tr '\t' ':' > F_filtered1chrpos.txt
bcftools view -H /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v2*/v2freebayes_filtered_sorted.vcf.gz | cut -f 1-2 |tr '\t' ':' > F_filtered2chrpos.txt
bcftools view -H /proj/uppstore2017199/b2016119_nobackup/ayushi/10_v3_variant_filtering_freebayes/v3freebayes_filtered_sorted.vcf.gz | cut -f 1-2 |tr '\t' ':' > F_filtered3chrpos.txt
# output files saved in GATK folder
```

```bash
# Load the VennDiagram package and other dependent packages
library(grid)
library(futile.logger)
library(VennDiagram)

#setwd("C:/Users/Ayushi Pathak/Downloads")

Gfile1<-'C:/Users/Ayushi Pathak/Downloads/GATK_filteredPass1chrpos.txt'
Gfile2<-'C:/Users/Ayushi Pathak/Downloads/GATK_filteredPass2chrpos.txt'
Gfile3<-'C:/Users/Ayushi Pathak/Downloads/GATK_filteredPass3chrpos.txt'
Ffile1<-'C:/Users/Ayushi Pathak/Downloads/F_filtered1chrpos.txt'
Ffile2<-'C:/Users/Ayushi Pathak/Downloads/F_filtered2chrpos.txt'
Ffile3<-'C:/Users/Ayushi Pathak/Downloads/F_filtered3chrpos.txt'
CDS<-'C:/Users/Ayushi Pathak/Downloads/DNA_filtered1chrpos.txt'
CDSExon<-'C:/Users/Ayushi Pathak/Downloads/cdsExonChrPos.txt'
GENOME<-'C:/Users/Ayushi Pathak/Downloads/allGenomeChrPos.txt'

G1<-readLines(Gfile1, warn = FALSE, encoding = "UTF-8")
G2<-readLines(Gfile2, warn = FALSE, encoding = "UTF-8")
G3<-readLines(Gfile3, warn = FALSE, encoding = "UTF-8")
F1<-readLines(Ffile1, warn = FALSE, encoding = "UTF-8")
F2<-readLines(Ffile2, warn = FALSE, encoding = "UTF-8")
F3<-readLines(Ffile3, warn = FALSE, encoding = "UTF-8")
DNA1<-readLines(CDS, warn = FALSE, encoding = "UTF-8")
DNA2<-readLines(CDSExon, warn = FALSE, encoding = "UTF-8")
DNA3<-readLines(GENOME, warn = FALSE, encoding = "UTF-8")

# Create a list of sets for GCS region
chromosome_sets1<-list(G1,F1,DNA1)
chromosome_sets2<-list(G2,F2,DNA1)
chromosome_sets3<-list(G3,F3,DNA1)

# Create Venn diagram
diagram1 <- venn.diagram(chromosome_sets1, category.names = c("GatkRun1","FreebayesRun1","DNA"), filename = "1cds.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))
diagram2<-venn.diagram(chromosome_sets2, category.names = c("GATKRun2","FreebayesRun2","DNA"), filename = "2cds.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram3<-venn.diagram(chromosome_sets3, category.names = c("GATKRun3","FreebayesRun3","DNA"), filename = "3cds.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        

# Create a list of sets for GCSExon region
chromosome_sets1.1<-list(G1,F1,DNA2)
chromosome_sets2.1<-list(G2,F2,DNA2)
chromosome_sets3.1<-list(G3,F3,DNA2)

# Create Venn diagram
diagram1.1<-venn.diagram(chromosome_sets1.1, category.names = c("GatkRun1","FreebayesRun1","DNA"), filename = "1cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))
diagram2.1<-venn.diagram(chromosome_sets2.1, category.names = c("GATKRun2","FreebayesRun2","DNA"), filename = "2cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram3.1<-venn.diagram(chromosome_sets3.1, category.names = c("GATKRun3","FreebayesRun3","DNA"), filename = "3cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram1.12<-venn.diagram(list(F1,DNA2), category.names = c("F1","DNA+Exon"), filename = "F1cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))
diagram1.22<-venn.diagram(list(F2,DNA2), category.names = c("F1","DNA+Exon"), filename = "F2cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))
diagram1.32<-venn.diagram(list(F3,DNA2), category.names = c("F3","DNA+Exon"), filename = "F3cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))

# Create a list of sets for Genome
chromosome_sets1.2<-list(G1,F1,DNA3)
chromosome_sets2.2<-list(G2,F2,DNA3)
chromosome_sets3.2<-list(G3,F3,DNA3)

# Create Venn diagram
diagram1.2<-venn.diagram(chromosome_sets1.2, category.names = c("GatkRun1","FreebayesRun1","DNA"), filename = "1genome.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))
diagram2.2<-venn.diagram(chromosome_sets2.2, category.names = c("GATKRun2","FreebayesRun2","DNA"), filename = "2genome.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram3.2<-venn.diagram(chromosome_sets3.2, category.names = c("GATKRun3","FreebayesRun3","DNA"), filename = "3genome.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram1.21<-venn.diagram(list(F1,DNA3), category.names = c("F1","WGS"), filename = "F1genome.png",  col="Black",fill=c("#29ADB2", "#C5E898"))
diagram2.22<-venn.diagram(list(F2,DNA3), category.names = c("F2","WGS"), filename = "F2genome.png",  col="Black",fill=c("#29ADB2", "#C5E898"))                        
diagram3.23<-venn.diagram(list(F3,DNA3), category.names = c("F3","WGS"), filename = "F3genome.png",  col="Black",fill=c("#29ADB2", "#C5E898"))     
diagram1.24<-venn.diagram(list(G1,DNA3), category.names = c("G1","WGS"), filename = "G1genome.png",  col="Black",fill=c("#29ADB2", "#C5E898"))
diagram2.25<-venn.diagram(list(G2,DNA3), category.names = c("G2","WGS"), filename = "G2genome.png",  col="Black",fill=c("#29ADB2", "#C5E898"))                        
diagram3.26<-venn.diagram(list(G3,DNA3), category.names = c("G3","WGS"), filename = "G3genome.png",  col="Black",fill=c("#29ADB2", "#C5E898"))
```

## Specificity, Sensitivity and Precision Analysis

```bash
# Install and load ggplot2 if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Example counts of true positives, false positives, true negatives, and false negatives for six runs
runs <- c("GATK1", "GATK2", "GATK3", "Freebayes1", "Freebayes2", "Freebayes3")

TP <- c(83245,
        50236,
        16069,
        111281,
        101398,
        75048)
FP <- c(635039,
        421622,
        134774,
        111476,
        74865,
        29475)
TN <- c(150, 160, 155, 145, 162, 158)
FN <- c(431091,
        464100,
        498267,
        403055,
        412938,
        439288)

# Calculate sensitivity and precision for each run
sensitivity <- TP / (TP + FN)
precision <- TP / (TP + FP)

# Define custom colors for Sensitivity and Precision
custom_colors <- c("#1F4E79", "#29ADB2")

# Create a data frame for plotting
data <- data.frame(
  Run = rep(runs, 2),
  Metric = rep(c("Sensitivity", "Precision"), each = 6),
  Value = c(sensitivity, precision)
  #Count = rep(c(TP, FP), times = 6)
)

# Create a bar plot using ggplot2 with custom colors
ggplot(data, aes(x = Metric, y = Value, fill = Metric, group = Run)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  ylim(0, 1) +
  labs(title = "Performance Metrics for Six Runs", y = "Metric Value") +
  facet_wrap(~ Run, scales = "free_x", ncol = 6) +
  theme_minimal() +
  theme(strip.text.x = element_text(size = 10, face = "bold"))
```
