# READ ME: Allele-Specific Expression in Bank Voles


### Variant calling and filtering analysis for RNA-seq data of Bank Voles. This pipeline is a part of the Allele Specific Expression pipeline. | LU Master’s Thesis in Bioinformatics 2023 | BINP51

---

# Uppmax Set-up

What all Simple Linux Utility for Resource Management (SLURM) does?
- job submission using sbatch command
- resource allocation
- job scheduling
- job monitoring

### Basic SLURM Commands

sbatch : run a bash script

```bash
sbatch my_job_script_file.sh
```

interactive : run an interactive session

```bash
interactive -A <project_name> -t DD-hh:mm:ss
```

jobinfo : to get status and other information about job

```bash
jobinfo -u <user_name>
#OR
jobinfo -A account
```

scancel : to cancel a job

```bash
scancel -n <job_name>
```

### Script Components

#!/bin/bash -l : shebang to make file executable
#SBATCH -A snic2022-5-71 : project name
#SBATCH -p node : nodes
#SBATCH -n 20 : number of nodes
#SBATCH -t 2-00:00:00 : time for the run
#SBATCH -J hifiasm_wikstroemia_hmc20s55.job : job name; it will help to search the project status
#SBATCH -o hifiasm_wikstroemia_hmc20s55.out : output file name
#SBATCH -e hifiasm_wikstroemia_hmc20s55.err : : error file name
#SBATCH [--mail-user=ayushipathakofficial@gmail.com](mailto:--mail-user=ayushipathakofficial@gmail.com) email for the job updates like queued, running, ending and cancelled
#SBATCH --mail-type=ALL:  what kind of mailing list you prefer

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
cd 2_trimmed_reads/
fastqc *.fastq
multiqc *
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

GATK does not work well for tools outside its environment. So use DBImport.

Prepare the interval list

```bash
cd /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants
module load bioinfo-tools
module load bcftools/1.9
# extracting chr list from all the vcf files
files=$(ls *.vcf.gz); for file in $files; do ls $file | bcftools query -H -f "%CHROM\n" $file ; done > gatk.chr.txt
# files=$(ls *.vcf.gz); for file in $files; do ls $file; bcftools view -H $file | bcftools query -f '%CHROM\n' | sort | uniq > $file.txt ; done
# combining the list, sorting and uniq
#cat *.uniq > gatk.chr.txt

#sorting and uniq
cat gatk.chr.txt | sort | uniq -c > gatk.chr.txt.stats
cat gatk.chr.txt | sort | uniq > gatk.chr.txt.sorted.uniq
```

Combined the **Joint Genotyping** step in this script as it needs same temporary directory.

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 00-10:00:00
#SBATCH -J 9.2.4_mergingFiles_DBImport.job
#SBATCH -o 9.2.4_mergingFiles_DBImport.out
#SBATCH -e 9.2.4_mergingFiles_DBImport.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#load the modules
module load bioinfo-tools
module load  GATK/4.3.0.0

#Provide the directory of the input and output files
in_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/calledVariants'
out_dir='/proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling'
res='/proj/uppstore2017199/b2016119_nobackup/ayushi/resources'

#Copy input files to the working directory
cd $in_dir
cp *.rawVariants.vcf.gz $SNIC_TMP
cp *.rawVariants.vcf.idx $SNIC_TMP
cp *.rawVariants.vcf.gz.tbi $SNIC_TMP
cd $res
cp GCF_902806735.1_Bank_vole1_10x_genomic.fna $SNIC_TMP
cp GCF_902806735.1_Bank_vole1_10x_genomic.fna.fai $SNIC_TMP
cp GCF_902806735.1_Bank_vole1_10x_genomic.dict $SNIC_TMP

cd $SNIC_TMP

gatk GenomicsDBImport \
                        -V P6207_101.rawVariants.vcf.gz \
                        -V P6207_104.rawVariants.vcf.gz \
                        -V P6207_105.rawVariants.vcf.gz \
                        -V P6207_109.rawVariants.vcf.gz \
                        -V P6207_112.rawVariants.vcf.gz \
                        -V P6207_113.rawVariants.vcf.gz \
                        -V P6207_116.rawVariants.vcf.gz \
                        -V P6207_117.rawVariants.vcf.gz \
                        -V P6207_120.rawVariants.vcf.gz \
                        -V P6207_121.rawVariants.vcf.gz \
                        -V P6207_123.rawVariants.vcf.gz \
                        -V P6207_125.rawVariants.vcf.gz \
                        -V P6207_128.rawVariants.vcf.gz \
                        -V P6207_129.rawVariants.vcf.gz \
                        -V P6207_132.rawVariants.vcf.gz \
                        -V P6207_134.rawVariants.vcf.gz \
                        -V P6207_135.rawVariants.vcf.gz \
                        -V P6207_136.rawVariants.vcf.gz \
                        --genomicsdb-workspace-path /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/DBImport/my_database \
      --tmp-dir $SNIC_TMP \
                        --merge-input-intervals \
      -L NC_024538.1 -L NW_025965004.1 -L NW_025965005.1 -L NW_025965006.1 -L NW_025965007.1 -L NW_025965008.1 -L NW_025965009.1 -L NW_025965010.1 -L NW_025965011.1 -L NW_025965012.1 -L NW_025965013.1 -L NW_025965014.1 -L NW_025965015.1 -L NW_025965016.1 -L NW_025965017.1 -L NW_025965018.1 -L NW_025965019.1 -L NW_025965020.1 -L NW_025965021.1 -L NW_025965022.1 -L NW_025965023.1 -L NW_025965024.1 -L NW_025965025.1 -L NW_025965026.1 -L NW_025965027.1 -L NW_025965028.1 -L NW_025965029.1 -L NW_025965030.1 -L NW_025965031.1 -L NW_025965032.1 -L NW_025965033.1 -L NW_025965034.1 -L NW_025965035.1 -L NW_025965036.1 -L NW_025965037.1 -L NW_025965038.1 -L NW_025965039.1 -L NW_025965040.1 -L NW_025965041.1 -L NW_025965042.1 -L NW_025965043.1 -L NW_025965044.1 -L NW_025965045.1 -L NW_025965046.1 -L NW_025965047.1 -L NW_025965048.1 -L NW_025965049.1 -L NW_025965050.1 -L NW_025965051.1 -L NW_025965052.1 -L NW_025965053.1 -L NW_025965054.1 -L NW_025965055.1 -L NW_025965056.1 -L NW_025965057.1 -L NW_025965058.1 -L NW_025965059.1 -L NW_025965060.1 -L NW_025965061.1 -L NW_025965062.1 -L NW_025965063.1 -L NW_025965064.1 -L NW_025965065.1 -L NW_025965066.1 -L NW_025965067.1 -L NW_025965068.1 -L NW_025965069.1 -L NW_025965070.1 -L NW_025965071.1 -L NW_025965072.1 -L NW_025965073.1 -L NW_025965074.1 -L NW_025965075.1 -L NW_025965076.1 -L NW_025965077.1 -L NW_025965078.1 -L NW_025965079.1 -L NW_025965080.1 -L NW_025965081.1 -L NW_025965082.1 -L NW_025965083.1 -L NW_025965084.1 -L NW_025965085.1 -L NW_025965086.1 -L NW_025965087.1 -L NW_025965088.1 -L NW_025965089.1 -L NW_025965090.1 -L NW_025965091.1 -L NW_025965092.1 -L NW_025965093.1 -L NW_025965094.1 -L NW_025965095.1 -L NW_025965096.1 -L NW_025965097.1 -L NW_025965098.1 -L NW_025965099.1 -L NW_025965100.1 -L NW_025965101.1 -L NW_025965102.1 -L NW_025965103.1 -L NW_025965104.1 -L NW_025965105.1 -L NW_025965106.1 -L NW_025965107.1 -L NW_025965108.1 -L NW_025965109.1 -L NW_025965110.1 -L NW_025965111.1 -L NW_025965112.1 -L NW_025965113.1 -L NW_025965114.1 -L NW_025965115.1 -L NW_025965116.1 -L NW_025965117.1 -L NW_025965118.1 -L NW_025965119.1 -L NW_025965120.1 -L NW_025965121.1 -L NW_025965122.1 -L NW_025965123.1 -L NW_025965124.1 -L NW_025965125.1 -L NW_025965126.1 -L NW_025965127.1 -L NW_025965128.1 -L NW_025965129.1 -L NW_025965130.1 -L NW_025965131.1 -L NW_025965132.1 -L NW_025965133.1 -L NW_025965134.1 -L NW_025965135.1 -L NW_025965136.1 -L NW_025965137.1 -L NW_025965138.1 -L NW_025965139.1 -L NW_025965140.1 -L NW_025965141.1 -L NW_025965142.1 -L NW_025965143.1 -L NW_025965144.1 -L NW_025965145.1 -L NW_025965146.1 -L NW_025965147.1 -L NW_025965148.1 -L NW_025965149.1 -L NW_025965150.1 -L NW_025965151.1 -L NW_025965152.1 -L NW_025965153.1 -L NW_025965154.1 -L NW_025965155.1 -L NW_025965156.1 -L NW_025965157.1 -L NW_025965158.1 -L NW_025965159.1 -L NW_025965160.1 -L NW_025965161.1 -L NW_025965162.1 -L NW_025965163.1 -L NW_025965164.1 -L NW_025965165.1 -L NW_025965166.1 -L NW_025965167.1 -L NW_025965168.1 -L NW_025965169.1 -L NW_025965170.1 -L NW_025965171.1 -L NW_025965172.1 -L NW_025965173.1 -L NW_025965174.1 -L NW_025965175.1 -L NW_025965176.1 -L NW_025965177.1 -L NW_025965178.1 -L NW_025965179.1 -L NW_025965180.1 -L NW_025965181.1 -L NW_025965182.1 -L NW_025965183.1 -L NW_025965184.1 -L NW_025965185.1 -L NW_025965186.1 -L NW_025965187.1 -L NW_025965188.1 -L NW_025965189.1 -L NW_025965190.1 -L NW_025965191.1 -L NW_025965192.1 -L NW_025965193.1 -L NW_025965194.1 -L NW_025965195.1 -L NW_025965196.1 -L NW_025965197.1 -L NW_025965198.1 -L NW_025965199.1 -L NW_025965200.1 -L NW_025965201.1 -L NW_025965202.1 -L NW_025965203.1 -L NW_025965204.1 -L NW_025965205.1 -L NW_025965206.1 -L NW_025965207.1 -L NW_025965208.1 -L NW_025965209.1 -L NW_025965210.1 -L NW_025965211.1 -L NW_025965212.1 -L NW_025965213.1 -L NW_025965214.1 -L NW_025965215.1 -L NW_025965216.1 -L NW_025965217.1 -L NW_025965218.1 -L NW_025965219.1 -L NW_025965220.1 -L NW_025965221.1 -L NW_025965222.1 -L NW_025965223.1 -L NW_025965224.1 -L NW_025965225.1 -L NW_025965226.1 -L NW_025965227.1 -L NW_025965228.1 -L NW_025965229.1 -L NW_025965230.1 -L NW_025965231.1 -L NW_025965232.1 -L NW_025965233.1 -L NW_025965234.1 -L NW_025965235.1 -L NW_025965236.1 -L NW_025965237.1 -L NW_025965238.1 -L NW_025965239.1 -L NW_025965240.1 -L NW_025965241.1 -L NW_025965242.1 -L NW_025965243.1 -L NW_025965244.1 -L NW_025965245.1 -L NW_025965246.1 -L NW_025965247.1 -L NW_025965248.1 -L NW_025965249.1 -L NW_025965250.1 -L NW_025965251.1 -L NW_025965252.1 -L NW_025965253.1 -L NW_025965254.1 -L NW_025965255.1 -L NW_025965256.1 -L NW_025965257.1 -L NW_025965258.1 -L NW_025965259.1 -L NW_025965260.1 -L NW_025965261.1 -L NW_025965262.1 -L NW_025965263.1 -L NW_025965264.1 -L NW_025965265.1 -L NW_025965266.1 -L NW_025965267.1 -L NW_025965268.1 -L NW_025965269.1 -L NW_025965270.1 -L NW_025965271.1 -L NW_025965272.1 -L NW_025965273.1 -L NW_025965274.1 -L NW_025965275.1 -L NW_025965276.1 -L NW_025965277.1 -L NW_025965278.1 -L NW_025965279.1 -L NW_025965280.1 -L NW_025965281.1 -L NW_025965282.1 -L NW_025965283.1 -L NW_025965284.1 -L NW_025965285.1 -L NW_025965286.1 -L NW_025965287.1 -L NW_025965288.1 -L NW_025965289.1 -L NW_025965290.1 -L NW_025965291.1 -L NW_025965292.1 -L NW_025965293.1 -L NW_025965294.1 -L NW_025965295.1 -L NW_025965296.1 -L NW_025965297.1 -L NW_025965298.1 -L NW_025965299.1 -L NW_025965300.1 -L NW_025965301.1 -L NW_025965302.1 -L NW_025965303.1 -L NW_025965304.1 -L NW_025965305.1 -L NW_025965306.1 -L NW_025965307.1 -L NW_025965308.1 -L NW_025965309.1 -L NW_025965310.1 -L NW_025965311.1 -L NW_025965312.1 -L NW_025965313.1 -L NW_025965314.1 -L NW_025965315.1 -L NW_025965316.1 -L NW_025965317.1 -L NW_025965318.1 -L NW_025965319.1 -L NW_025965320.1 -L NW_025965321.1 -L NW_025965322.1 -L NW_025965323.1 -L NW_025965324.1 -L NW_025965325.1 -L NW_025965326.1 -L NW_025965327.1 -L NW_025965328.1 -L NW_025965329.1 -L NW_025965330.1 -L NW_025965331.1 -L NW_025965332.1 -L NW_025965333.1 -L NW_025965334.1 -L NW_025965335.1 -L NW_025965336.1 -L NW_025965337.1 -L NW_025965338.1 -L NW_025965339.1 -L NW_025965340.1 -L NW_025965341.1 -L NW_025965342.1 -L NW_025965343.1 -L NW_025965344.1 -L NW_025965345.1 -L NW_025965346.1 -L NW_025965347.1 -L NW_025965348.1 -L NW_025965349.1 -L NW_025965350.1 -L NW_025965351.1 -L NW_025965352.1 -L NW_025965353.1 -L NW_025965354.1 -L NW_025965355.1 -L NW_025965356.1 -L NW_025965357.1 -L NW_025965358.1 -L NW_025965359.1 -L NW_025965360.1 -L NW_025965361.1 -L NW_025965362.1 -L NW_025965363.1 -L NW_025965364.1 -L NW_025965365.1 -L NW_025965366.1 -L NW_025965367.1 -L NW_025965368.1 -L NW_025965369.1 -L NW_025965370.1 -L NW_025965371.1 -L NW_025965372.1 -L NW_025965373.1 -L NW_025965374.1 -L NW_025965375.1 -L NW_025965376.1 -L NW_025965377.1 -L NW_025965378.1 -L NW_025965379.1 -L NW_025965380.1 -L NW_025965381.1 -L NW_025965382.1 -L NW_025965383.1 -L NW_025965384.1 -L NW_025965385.1 -L NW_025965386.1 -L NW_025965387.1 -L NW_025965388.1 -L NW_025965389.1 -L NW_025965390.1 -L NW_025965391.1 -L NW_025965392.1 -L NW_025965393.1 -L NW_025965394.1 -L NW_025965395.1 -L NW_025965396.1 -L NW_025965397.1 -L NW_025965398.1 -L NW_025965399.1 -L NW_025965400.1 -L NW_025965401.1 -L NW_025965402.1 -L NW_025965403.1 -L NW_025965404.1 -L NW_025965405.1 -L NW_025965406.1 -L NW_025965407.1 -L NW_025965408.1 -L NW_025965409.1 -L NW_025965410.1 -L NW_025965411.1 -L NW_025965412.1 -L NW_025965413.1 -L NW_025965414.1 -L NW_025965415.1 -L NW_025965416.1 -L NW_025965417.1 -L NW_025965418.1 -L NW_025965419.1 -L NW_025965420.1 -L NW_025965421.1 -L NW_025965422.1 -L NW_025965423.1 -L NW_025965424.1 -L NW_025965425.1 -L NW_025965426.1 -L NW_025965427.1 -L NW_025965428.1 -L NW_025965429.1 -L NW_025965430.1 -L NW_025965431.1 -L NW_025965432.1 -L NW_025965433.1 -L NW_025965434.1 -L NW_025965435.1 -L NW_025965436.1 -L NW_025965437.1 -L NW_025965438.1 -L NW_025965439.1 -L NW_025965440.1 -L NW_025965441.1 -L NW_025965442.1 -L NW_025965443.1 -L NW_025965444.1 -L NW_025965445.1 -L NW_025965446.1 -L NW_025965447.1 -L NW_025965448.1 -L NW_025965449.1 -L NW_025965450.1 -L NW_025965451.1 -L NW_025965452.1 -L NW_025965453.1 -L NW_025965454.1 -L NW_025965455.1 -L NW_025965456.1 -L NW_025965457.1 -L NW_025965458.1 -L NW_025965459.1 -L NW_025965460.1 -L NW_025965461.1 -L NW_025965462.1 -L NW_025965463.1 -L NW_025965464.1 -L NW_025965465.1 -L NW_025965466.1 -L NW_025965467.1 -L NW_025965468.1 -L NW_025965469.1 -L NW_025965470.1 -L NW_025965471.1 -L NW_025965472.1 -L NW_025965473.1 -L NW_025965474.1 -L NW_025965475.1 -L NW_025965476.1 -L NW_025965477.1 -L NW_025965478.1 -L NW_025965479.1 -L NW_025965480.1 -L NW_025965481.1 -L NW_025965482.1 -L NW_025965483.1 -L NW_025965484.1 -L NW_025965485.1 -L NW_025965486.1 -L NW_025965487.1 -L NW_025965488.1 -L NW_025965489.1 -L NW_025965490.1 -L NW_025965491.1 -L NW_025965492.1 -L NW_025965493.1 -L NW_025965494.1 -L NW_025965495.1 -L NW_025965496.1 -L NW_025965497.1 -L NW_025965498.1 -L NW_025965499.1 -L NW_025965500.1 -L NW_025965501.1 -L NW_025965502.1 -L NW_025965503.1 -L NW_025965504.1 -L NW_025965505.1 -L NW_025965506.1 -L NW_025965507.1 -L NW_025965508.1 -L NW_025965509.1 -L NW_025965510.1 -L NW_025965511.1 -L NW_025965512.1 -L NW_025965513.1 -L NW_025965514.1 -L NW_025965515.1 -L NW_025965516.1 -L NW_025965517.1 -L NW_025965518.1 -L NW_025965519.1 -L NW_025965520.1 -L NW_025965521.1 -L NW_025965522.1 -L NW_025965523.1 -L NW_025965524.1 -L NW_025965525.1 -L NW_025965526.1 -L NW_025965527.1 -L NW_025965528.1 -L NW_025965529.1 -L NW_025965530.1 -L NW_025965531.1 -L NW_025965532.1 -L NW_025965533.1 -L NW_025965534.1 -L NW_025965535.1 -L NW_025965536.1 -L NW_025965537.1 -L NW_025965538.1 -L NW_025965539.1 -L NW_025965540.1 -L NW_025965541.1 -L NW_025965542.1 -L NW_025965543.1 -L NW_025965544.1 -L NW_025965545.1 -L NW_025965546.1 -L NW_025965547.1 -L NW_025965548.1 -L NW_025965549.1 -L NW_025965550.1 -L NW_025965551.1 -L NW_025965552.1 -L NW_025965553.1 -L NW_025965554.1 -L NW_025965555.1 -L NW_025965556.1 -L NW_025965557.1 -L NW_025965558.1 -L NW_025965559.1 -L NW_025965560.1 -L NW_025965561.1 -L NW_025965562.1 -L NW_025965563.1 -L NW_025965564.1 -L NW_025965565.1 -L NW_025965566.1 -L NW_025965567.1 -L NW_025965568.1 -L NW_025965569.1 -L NW_025965570.1 -L NW_025965571.1 -L NW_025965572.1 -L NW_025965573.1 -L NW_025965574.1 -L NW_025965575.1 -L NW_025965576.1 -L NW_025965577.1 -L NW_025965578.1 -L NW_025965579.1 -L NW_025965580.1 -L NW_025965581.1 -L NW_025965582.1 -L NW_025965583.1 -L NW_025965584.1 -L NW_025965585.1 -L NW_025965586.1 -L NW_025965587.1 -L NW_025965588.1 -L NW_025965589.1 -L NW_025965590.1 -L NW_025965591.1 -L NW_025965592.1 -L NW_025965593.1 -L NW_025965594.1 -L NW_025965595.1 -L NW_025965596.1 -L NW_025965597.1 -L NW_025965598.1 -L NW_025965599.1 -L NW_025965600.1 -L NW_025965601.1 -L NW_025965602.1 -L NW_025965603.1 -L NW_025965604.1 -L NW_025965605.1 -L NW_025965606.1 -L NW_025965607.1 -L NW_025965608.1 -L NW_025965609.1 -L NW_025965610.1 -L NW_025965611.1 -L NW_025965612.1 -L NW_025965613.1 -L NW_025965614.1 -L NW_025965615.1 -L NW_025965616.1 -L NW_025965617.1 -L NW_025965618.1 -L NW_025965619.1 -L NW_025965620.1 -L NW_025965621.1 -L NW_025965622.1 -L NW_025965623.1 -L NW_025965624.1 -L NW_025965625.1 -L NW_025965626.1 -L NW_025965627.1 -L NW_025965628.1 -L NW_025965629.1 -L NW_025965630.1 -L NW_025965631.1 -L NW_025965632.1 -L NW_025965633.1 -L NW_025965634.1 -L NW_025965635.1 -L NW_025965636.1 -L NW_025965637.1 -L NW_025965638.1

## for GenotypeGVCFs command look for the scratch directory on you uppmax account. Usually it gives you an error. Follow the code below of the error persist.
cp -r /proj/uppstore2017199/b2016119_nobackup/ayushi/9.2_GATK_variant_calling/DBImport/my_database  $SNIC_TMP

gatk GenotypeGVCFs \
                        -R GCF_902806735.1_Bank_vole1_10x_genomic.fna \
                        -V gendb://my_database \
                        --output $out_dir/DBImport/genotypeOutput.vcf.gz
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

FILTERING THE EXPRESSED GENES

using featureCounts

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 5
#SBATCH -t 10-00:00:00
#SBATCH -J expressedGenesForSensitivityAnlaysis.job
#SBATCH -o expressedGenesForSensitivityAnlaysis.out
#SBATCH -e expressedGenesForSensitivityAnlaysis.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load subread
module load star/2.7.9a

featureCounts -p -T 4 -a /proj/uppstore2017199/b2016119_nobackup/ayushi/resources/genomic.gtf -t exon -F GTF -g gene_id -o /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/featureCountsAnalysis
```

Filters

```jsx
grep -v '#' featureCountsAnalysis | grep -v '^Geneid' | awk '{for(i=6;i<=24;i++) if($i < 4) $i=0}1' featureCountsAnalysis > featureCountsAnalysis.filter1
```

```jsx
awk '{count=0; for(i=6; i<=24; i++) if($i == 0) count++; if(count < 13) print}' featureCountsAnalysis.filter1 > featureCountsAnalysis.filter2
```

```bash
/proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/featureCountsAnalysis.final
/proj/uppstore2017199/b2016119_nobackup/ayushi/genomic.gff

```

```bash
# Filter the gtf file to contain only the **exons** from these genes. This will give you an **annotation file for the exons from expressed genes.**
awk '$3 ~ /exon/' genomic.gtf > exon.gtf
# Convert the gtf file to bed file
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $1, $4 - 1, $5, $9, ".", $7}' exon.gtf > exon.bed
# Filtering the bed file to contain only expressed genes

#bcftools view -R exon.bed /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/featureCountsAnalysis.final | bcftools sort >  /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/expressed_annotation.txt
```

```bash

grep -v '^#' featureCountsAnalysis.filter2 | grep -v '^geneid' | cut -f 1 > listExpressedGenes.txt
cat listExpressedGenes.txt | tr '\t' ';' > listExpressedGenes2.txt
cat listExpressedGenes2.txt | tr ' ' ';' > listExpressedGenes3.txt
cut -d ';' -f 1 listExpressedGenes3.txt> listExpressedGenes4.txt
###
LOC125387327
LOC125398152
Elmod3
Capg
Sh2d6
LOC125405528
LOC125406016
Mat2a
Ggcx
LOC125408558
Vamp8
Vamp5
Rnf181
Tmem150a
CUNH2orf68
Usp39###
```

```bash
#!/bin/bash -l
#SBATCH -A naiss2023-22-97
#SBATCH -p node
#SBATCH -n 5
#SBATCH -t 12:00:00
#SBATCH -J expressed.job
#SBATCH -o expressed.out
#SBATCH -e expressed.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#for value in $(cat /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/expressed_genes.txt); do grep -E "gene_id \"$value\"" /proj/uppstore2017199/b2016119_nobackup/ayushi/exon.gtf; done > /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/expressed.gtf

while read -r value; do sed -n "/gene_id \"$value\"/p" /proj/uppstore2017199/b2016119_nobackup/ayushi/exon.gtf; done < /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/expressed_genes.txt > /proj/uppstore2017199/b2016119_nobackup/ayushi/20_expressedGenesForSensitivityAnlaysis/expressed.gtf
```

```bash
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $1, $4 - 1, $5, $9, ".", $7}' expressed.gtf > expressed.bed
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ {print $1, $4 - 1, $5, $9, ".", $7}' expressed6.gtf > expressed6.bed

module load bioinfo-tools bcftools
cd 20_expressedGenesForSensitivityAnlaysis
bcftools view -R expressed.bed ../freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz | bcftools sort | bgzip -c >  expressed.vcf.gz
bcftools view -R expressed6.bed ../freebayes.bv_chromium.filt1.decomposed.filt.snps.rep.all.vcf.gz | bcftools sort | bgzip -c >  expressed6.vcf.gz

bcftools view -H expressed.vcf.gz | cut -f 1-2 |tr '\t' ':' > expressedChrPos.txt
bcftools view -H expressed6.vcf.gz | cut -f 1-2 | tr '\t' ':' > expressedChrPos6.txt

cd Downloads
scp ayuship@rackham.uppmax.uu.se:/proj/uppstore2017199/b2016119_nobackup/ayushi/20*/expressedChrPos.txt .
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

setwd("C:/Users/Ayushi Pathak/Downloads")

Gfile1<-'C:/Users/Ayushi Pathak/Downloads/GATK_filteredPass1chrpos.txt'
Gfile2<-'C:/Users/Ayushi Pathak/Downloads/GATK_filteredPass2chrpos.txt'
Gfile3<-'C:/Users/Ayushi Pathak/Downloads/GATK_filteredPass3chrpos.txt'
Ffile1<-'C:/Users/Ayushi Pathak/Downloads/F_filtered1chrpos.txt'
Ffile2<-'C:/Users/Ayushi Pathak/Downloads/F_filtered2chrpos.txt'
Ffile3<-'C:/Users/Ayushi Pathak/Downloads/F_filtered3chrpos.txt'
CDS<-'C:/Users/Ayushi Pathak/Downloads/DNA_filtered1chrpos.txt'
CDSExon<-'C:/Users/Ayushi Pathak/Downloads/cdsExonChrPos.txt'
GENOME<-'C:/Users/Ayushi Pathak/Downloads/allGenomeChrPos.txt'
expressed<-'C:/Users/Ayushi Pathak/Downloads/expressedChrPos.txt'

G1<-readLines(Gfile1, warn = FALSE, encoding = "UTF-8")
G2<-readLines(Gfile2, warn = FALSE, encoding = "UTF-8")
G3<-readLines(Gfile3, warn = FALSE, encoding = "UTF-8")
F1<-readLines(Ffile1, warn = FALSE, encoding = "UTF-8")
F2<-readLines(Ffile2, warn = FALSE, encoding = "UTF-8")
F3<-readLines(Ffile3, warn = FALSE, encoding = "UTF-8")
DNA1<-readLines(CDS, warn = FALSE, encoding = "UTF-8")
DNA2<-readLines(CDSExon, warn = FALSE, encoding = "UTF-8")
DNA3<-readLines(GENOME, warn = FALSE, encoding = "UTF-8")
expressed<-readLines(expressed, warn = FALSE, encoding = "UTF-8")

# Create a list of sets for CDS region
chromosome_sets1<-list(G1,F1,DNA1)
chromosome_sets2<-list(G2,F2,DNA1)
chromosome_sets3<-list(G3,F3,DNA1)

# Create Venn diagram
diagram1 <- venn.diagram(chromosome_sets1, category.names = c("GatkRun1","FreebayesRun1","DNA"), filename = "1cds.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))
diagram2<-venn.diagram(chromosome_sets2, category.names = c("GATKRun2","FreebayesRun2","DNA"), filename = "2cds.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram3<-venn.diagram(chromosome_sets3, category.names = c("GATKRun3","FreebayesRun3","DNA"), filename = "3cds.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        

# Create a list of sets for CDSExon region
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
diagram1.14<-venn.diagram(list(G1,DNA2), category.names = c("G1","DNA+Exon"), filename = "G1cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))
diagram1.25<-venn.diagram(list(G2,DNA2), category.names = c("G1","DNA+Exon"), filename = "G2cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))
diagram1.36<-venn.diagram(list(G3,DNA2), category.names = c("G3","DNA+Exon"), filename = "G3cdsexon.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))

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

#####################################

# Create a list of sets for expressed region
chromosome_sets1<-list(G1,F1,expressed)
chromosome_sets2<-list(G2,F2,expressed)
chromosome_sets3<-list(G3,F3,expressed)

# Create Venn diagram
diagram1 <- venn.diagram(chromosome_sets1, category.names = c("GatkRun1","FreebayesRun1","ExpressedGenes"), filename = "1ex6.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))
diagram2<-venn.diagram(chromosome_sets2, category.names = c("GATKRun2","FreebayesRun2","ExpressedGenes"), filename = "2ex6.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        
diagram3<-venn.diagram(chromosome_sets3, category.names = c("GATKRun3","FreebayesRun3","ExpressedGenes"), filename = "3ex6.png",  col="Black",fill=c("#1F4E79", "#29ADB2", "#C5E898"))                        

chromosome_sets1.1<-list(F1,expressed)
chromosome_sets2.1<-list(F2,expressed)
chromosome_sets3.1<-list(F3,expressed)

diagram1<-venn.diagram(chromosome_sets1.1, category.names = c("Freebayes1","ExpressedGenes"), filename = "1ex6.1.png",  col="Black",fill=c("#1F4E79", "#29ADB2" ))
diagram2<-venn.diagram(chromosome_sets2.1, category.names = c("Freebayes2","ExpressedGenes"), filename = "2ex6.1.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))                        
diagram3<-venn.diagram(chromosome_sets3.1, category.names = c("Freebayes3","ExpressedGenes"), filename = "3ex6.1.png",  col="Black",fill=c("#1F4E79", "#29ADB2"))
```

## Specificity and Precision Analysis

```bash
# Install and load ggplot2 if not already installed
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

# Example counts of true positives, false positives, true negatives, and false negatives for six runs
runs <- c("GATK1", "GATK2", "GATK3", "Freebayes1", "Freebayes2", "Freebayes3")

# For Expressed Genes
TP <- c(82965,50032,16035,111262,101398,75048)

FP <- c(635319,421806,134808,111495,74865,29475)

FN <- c(339592,372525,406522,311295,321159,347509)

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
