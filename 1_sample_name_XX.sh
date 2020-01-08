#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name='sample_name_XX'
#SBATCH --output='sample_name_XX.log'
#SBATCH --mem-per-cpu=100G
#SBATCH --partition=shortq
#SBATCH -m block
/bin/hostname


module load star/2.5.2b
module load htslib/1.3.1

genome=hg38

# path to a directory with raw (unmapped) bam files
local_raw_data_location=MGs_RNAseq/raw_bsf_data

# specify output directory path
mapping_results_dir=MGs_RNAseq/mapping_results

# specify a sample name
sample_name=sample_name_XX

# create output directories 
mkdir -p ${mapping_results_dir}
mkdir -p ${mapping_results_dir}/fastq
mkdir -p ${mapping_results_dir}/fastqc
mkdir -p ${mapping_results_dir}/star_mapping_htseq



# Generate fastq file
fastq_out=$(echo ${mapping_results_dir}/fastq/${sample_name}_R1.fastq)

`/cm/shared/apps/samtools/1.3.1/bin/samtools view ${local_raw_data_location}/${sample_name}.bam | awk -v fastq_out=${fastq_out} '{ print "@"$1"\n"$10"\n+\n"$11 > fastq_out; }'`




# Trimming
# Generate *trimmed.fastq

`/cm/shared/apps/java/jdk/1.7.0_80/bin/java -Xmx60000m -jar /cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar SE -phred33 -threads 2 ${mapping_results_dir}/fastq/${sample_name}_R1.fastq ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq HEADCROP:13 ILLUMINACLIP:/data/groups/lab_bock/shared/resources/adapters/epignome_adapters_2_add.fa:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:18`


# Generate Fastqc report
`fastqc --noextract --outdir ${mapping_results_dir}/fastqc/ ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq`



# STAR mapping to hg38/GRCh38 assembly of the human reference genome and to SIRV-Set 3 sequences

`STAR --runThreadN 4 --genomeDir /data/groups/lab_winter/reference_files/indices/STAR/ --readFilesIn ${mapping_results_dir}/fastq/${sample_name}_R1_trimmed.fastq  --outFileNamePrefix ${mapping_results_dir}/star_mapping_htseq/${sample_name}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate`



# Generate `sample_name_XX/read_counts.txt`

`htseq-count -m intersection-nonempty -s yes -f bam -r pos ${mapping_results_dir}/star_mapping_htseq/${sample_name}_Aligned.sortedByCoord.out.bam /data/groups/lab_winter/reference_files/genomes/hg38_SIRV.gtf > ${mapping_results_dir}/star_mapping_htseq/${sample_name}_read_counts.txt`




 