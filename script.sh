# ###########################################
# #  - NGS-2 - Assignment                   #
# #  - April 29,2019                        #
# #  - St.Name : Ahmed Omar                 #
# #  - Nile University                      #
# ###########################################
# #!/bin/bash

# ----------------------------------------------------
# # Install missing packages 
# conda install -c bioconda star 
# sudo apt install picard-tools

# # ----------------------------------------------------
# # define folders path 
Assignment=~/Desktop/assignment
STAR=/home/ahmed/miniconda3/envs/ngs2/bin/STAR
picard=/home/ahmed/miniconda3/envs/ngs2/bin/picard
gatk=/home/ahmed/miniconda3/envs/ngs2/bin/gatk
samtools=/home/ahmed/miniconda3/envs/ngs2/bin/gatk

# # -----------------------------------------------------
# # rename files 
# mv $Assignment/data/S2_SRR8797509_1.part_001.part_001.fastq $Assignment/data/S1_L001_R1_001.fastq
# mv $Assignment/data/S2_SRR8797509_2.part_001.part_001.fastq $Assignment/data/S1_L001_R2_001.fastq

# mv $Assignment/data/S1_SRR8797509_1.part_001.part_001.fastq $Assignment/data/S2_L001_R1_001.fastq
# mv $Assignment/data/S1_SRR8797509_2.part_001.part_001.fastq $Assignment/data/S2_L001_R2_001.fastq

# # ----------------------------------------------------
# # Downloading reference
# mkdir -p $Assignment/ref/
# wget -P $Assignment/ref/ http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa

# # ----------------------------------------------------
# ## Alignment using STAR
# # Indexing step 1$Assignment 
# GENOME_DIR="$Assignment/ref/STAR_index/"
# mkdir -p $GENOME_DIR && cd $GENOME_DIR
# ln -s $Assignment/ref/chr22_with_ERCC92.fa .
# cd -
# $STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $Assignment/ref/chr22_with_ERCC92.fa  --runThreadN 4

# # ----------------------------------------------------
# # execut Alignment steps 1: 
# PASS_1_DIR=$Assignment/results/STAR_results/1pass
# mkdir -p $PASS_1_DIR
# for R1 in $Assignment/data/*_R1_001.fastq;do
#     mkdir -p $PASS_1_DIR/$(basename $R1 _R1_001.fastq)
#     cd $PASS_1_DIR/$(basename $R1 _R1_001.fastq)
#     R2=$(echo $R1 | sed 's/_R1_/_R2_/')
#     echo $R1 $R2
#     STAR --genomeDir $Assignment/ref/star-index/ --readFilesIn $Assignment/data/$R1 $Assignment/data/$R2 --runThreadN 4
#     cd -
# done

# # ----------------------------------------------------
# ## Indexing step 1$Assignment 
# for S in {1$Assignment2}; do
#     GENOME_DIR="$Assignment/ref/STAR-index-2/S"$S
#     mkdir -p $GENOME_DIR && cd $GENOME_DIR
#     cd -
#     $STAR --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $Assignment/ref/chr22_with_ERCC92.fa --sjdbFileChrStartEnd $Assignment/STAR_results/1pass/S${S}_L001/SJ.out.tab --sjdbOverhang 75 --runThreadN 4
# done

# # ----------------------------------------------------
# # execut Alignment steps 2: 
# PASS_2_DIR=$Assignment/results/STAR_results/2pass
# mkdir -p $PASS_2_DIR
# for R1 in $Assignment/data/*_R1_001.fastq;do
#     mkdir -p $PASS_2_DIR/$(basename $R1 _R1_001.fastq)
#     cd $PASS_2_DIR/$(basename $R1 _R1_001.fastq)
#     R2=$(echo $R1 | sed 's/_R1_/_R2_/')
#     echo $R1 $R2
#     echo $Assignment/ref/star-index-2/$(basename $R1 | cut -d"_" -f1)
#     STAR --genomeDir $Assignment/ref/star-index-2/$(basename $R1 | cut -d"_" -f1) --readFilesIn $Assignment/$R1 $Assignment/$R2 --runThreadN 4
#     cd -
# done

# # ----------------------------------------------------
# # Add [Read group information] 
# mkdir -p $Assignment/results/picard
# for i in S1_L001 S2_L001;do
#     SM=$(basename $i | cut -d"_" -f1)                                          ##sample ID
#     LB=$i                                        ##library ID
#     PL="Illumina"                                                           ##platform (e.g. illumina, solid)
#     RGID=$(cat $Assignment/data/${i}_R1_001.fastq | head -n1 | sed 's/ /_/g' | cut -d "_" -f1)       ##read group identifier 
#     PU=$RGID.$LB                                                            ##Platform Unit
#     echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"

#     $picard AddOrReplaceReadGroups I=$Assignment/STAR_results/2pass/$i/Aligned.out.sam O=$Assignment/results/picard/${SM}_rg_added_sorted.bam SO=coordinate RGID=$RGID RGLB=$LB RGPL=$PL RGPU=$PU RGSM=$SM 

#     $picard MarkDuplicates I=$Assignment/results/picard/${SM}_rg_added_sorted.bam O=$Assignment/results/picard/${SM}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$Assignment/results/picard/${SM}_output.metrics 
# done


# # ----------------------------------------------------
# #Creating the fasta sequence dictionary file
# picard-tools CreateSequenceDictionary \
#     R=$Assignment/ref/chr22_with_ERCC92.fa\
#     O=$Assignment/ref/chr22_with_ERCC92.fa.dict
# $samtools faidx $Assignment/ref/chr22_with_ERCC92.fa

# # # ----------------------------------------------------
# #Split'N'Trim
# mkdir -p $Assignment/results/gatk_res
# for i in S1 S2;do
#     $gatk SplitNCigarReads -R $Assignment/ref/chr22_with_ERCC92.fa -I $Assignment/results/picard/${i}_dedupped.bam -O $Assignment/results/gatk_res/${i}_split.bam
# done

# download vcf file fro reference 
# wget -P $Assignment/ref ftp://ftp.ensembl.org/pub/release-96/variation/vcf/homo_sapiens/homo_sapiens-chr22.vcf.gz

# cd $Assignment/ref
# gunzip -k homo_sapiens-chr22.vcf.gz 
# grep "^#" homo_sapiens-chr22.vcf > chr22.vcf
# grep "^22" homo_sapiens-chr22.vcf | sed 's/^22/chr22/' >> chr22.vcf
# gatk IndexFeatureFile -F chr22.vcf
# cd -

for sample in S1 S2;do
    $gatk --java-options "-Xmx2G" BaseRecalibrator -R $Assignment/ref/chr22_with_ERCC92.fa -I $Assignment/results/picard/${sample}_dedupped.bam --known-sites $Assignment/ref/chr22.vcf -O $Assignment/results/gatk_res/${sample}.report
    $gatk --java-options "-Xmx2G" ApplyBQSR -R $Assignment/ref/chr22_with_ERCC92.fa -I $Assignment/results/picard/${sample}_dedupped.bam -bqsr $Assignment/results/gatk_res/${sample}.report -O $Assignment/results/gatk_res/${sample}.bqsr.bam --add-output-sam-program-record --emit-original-quals
done 

$gatk HaplotypeCaller -R $Assignment/ref/chr22_with_ERCC92.fa -I $Assignment/results/picard/S1_dedupped.bam -O $Assignment/results/gatk_res/S1_dedupped.vcf
$gatk HaplotypeCaller -R $Assignment/ref/chr22_with_ERCC92.fa -I $Assignment/results/picard/S2_dedupped.bam -O $Assignment/results/gatk_res/S2_dedupped.vcf

$gatk VariantFiltration -R $Assignment/ref/chr22_with_ERCC92.fa -V $Assignment/results/gatk_res/S1_dedupped.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $Assignment/results/gatk_res/S1.vcf
$gatk VariantFiltration -R $Assignment/ref/chr22_with_ERCC92.fa -V $Assignment/results/gatk_res/S1_dedupped.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $Assignment/results/gatk_res/S1.vcf
