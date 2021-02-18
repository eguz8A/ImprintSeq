#!/bin/bash
############################################################################################################################
## IMPRINTSEQ_STEP 1_Fastq_to_MedianMethPerDMR
############################################################################################################################

#0. Define variables
hg19 = "location of hg19 bisulfite reference genome" #http://felixkrueger.github.io/Bismark/Docs/

#1. Library data
#Create list of samples and analysis folders
ls *_R1_001.fastq.gz > samples.txt
vim -c "%s/_R1_001.fastq.gz//g|wq" samples.txt
mkdir bam_files
mkdir methylation_files
mkdir CpGmeth_files
for k in `cat samples.txt`; do
    ## Creating Temp Files
	if [ -e .BS_${k}_TempFiles ]; then
        rm -r .BS_${k}_TempFiles
	fi
        mkdir .BS_${k}_TempFiles
	cd .BS_${k}_TempFiles
	cp ../${k}_R1_001.fastq.gz ${k}_R1_001.fastq.gz
	cp ../rpmt_5prime_10bp.py rpmt_5prime_10bp.py
	echo -e `date +[%D-%R]` "\tStarting Bisulfite sequencing analysis of ${k}"

#2. Quality and filtering
    #Quality        
	fastqc ${k}_R1_001.fastq.gz
    # Trimming using cores (-j) quality (-q) 30, adapter min-overlap (-O) and min-length (-m) 30. 
	cutadapt --adapter AGATCGGAAGAGC -m 30 -O 1 -q 30 -o ${k}_R1_trimmed.fq  ${k}_R1_001.fastq.gz  > /dev/null 2>&1
    # hard trim 10bp on 5' end (in-house script from NuGen-TECAN)
	python rpmt_5prime_10bp.py ${k}_R1_trimmed.fq
    #Quality
	fastqc ${k}_R1_trimmed_rpmt.fq

#3. Alignment, remove PCR duplicates and methylation extraction per CpG
    #Alignment
	bismark --pbat --multicore 8 -p 4 --bowtie2 --genome_folder $hg19 ${k}_R1_trimmed_rpmt.fq  > /dev/null 2>&1
    #Deduplication
	deduplicate_bismark --barcode --bam ${k}_R1_trimmed_rpmt_bismark_bt2.bam  > /dev/null 2>&1
    #Methylation extraction with Bismark
	bismark_methylation_extractor --multicore 8 --cytosine_report --bedGraph --gzip --genome_folder $hg19 ${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bam  > /dev/null 2>&1
	cp ${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bismark.cov.gz ../methylation_files/${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bismark.cov.gz


#4. Create a sorted index bam file for IGV or other softwares. 
    #Add read groups to bam file with picard to make bam file readable by samtools, igv and/or gatk
	java -Xmx100g -jar /data/Resources/Software/Javas/picard.jar AddOrReplaceReadGroups I=${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bam O=${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS.bam RGPL=illumina RGLB=LaneX RGPU=NONE RGSM=${k}  > /dev/null 2>&1
	samtools sort ${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS.bam > ${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS_sorted.bam
	samtools index ${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS_sorted.bam
	cp ${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS_sorted.bam ../bam_files/${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS_sorted.bam
	cp ${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS_sorted.bam.bai ../bam_files/${k}_trimmed_rpmt_bismark_bt2.dedup_RRBS_sorted.bam.bai
	echo -e `date +[%D-%R]` "\tBisulfite sequencing analysis of ${k} finished"

#5.Methylation analysis: Filtering by coverage 100 and by 63 imprinting DMRs
	echo -e `date +[%D-%R]` "\tCpG methylation data of ${k} filtered by coverage" | tee -a ../${k}_CpGfiltering.log
	echo -e `date +[%D-%R]` "\tNumber CpG before filtering in ${k}" | tee -a ../${k}_CpGfiltering.log
	gunzip ${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bismark.cov.gz
	wc -l ${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bismark.cov
	sed -i "s/^/${k}\t/" ${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bismark.cov
    #calculate and filter by coverage 100
	awk 'BEGIN {FS="\t"} {OFS="\t"} {print $0, $6+$7}' ${k}_R1_trimmed_rpmt_bismark_bt2.dedup_RRBS.bismark.cov | awk 'BEGIN {FS="\t"} {OFS="\t"} { if ($8 > 30) print }' > ${k}_Meth+DP100.txt
	echo -e `date +[%D-%R]` "\tNumber CpG DP > 100 in ${k}" | tee -a ../${k}_CpGfiltering.log
	wc -l ${k}_Meth+DP100.txt
    #Format file by sorting columns
	awk 'BEGIN {FS="\t"} {OFS="\t"} {print $2, $3, $4, $5, $6, $7, $8, $1}' ${k}_Meth+DP100.txt | awk 'BEGIN {FS="\t"} {gsub(/_/,"\t",$8); print}' | sort -k1,1 -k2,2n | awk 'BEGIN {FS=" "} {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8}' | sort -k1,1 -k2,2n > ${k}_Meth+DP100.sorted.txt
    #Filter by target regions, 63 imprinting DMRs
	bedtools intersect -wa -wb -a Targets_hg19.bed -b ${k}_Meth+DP100.sorted.txt | sort | uniq -u | awk 'BEGIN {FS=" "} {OFS="\t"} {print $5, $6, $7, $4, $8, $9, $10, $11, $12}' > ${k}_DP100_Targets.txt
	echo -e `date +[%D-%R]` "\tNumber CpG on Target with DP>100 in ${k}" | tee -a ../${k}_CpGfiltering.log
	wc -l ${k}_DP100_Targets.txt
	cp ${k}_DP100_Targets.txt ../CpGmeth_files/${k}_DP100_Targets.txt
	cd ..
done

#6. Combining all individual files in 1
cd .BS_${k}_TempFiles
cat *_DP100_Targets.txt > Library_DP100_Targets.txt
awk 'BEGIN {FS="\t"} {OFS="\t"} {print "chr" $1 ":" $2 "_" $4, $5, $9}' Library_DP100_Targets.txt > Library_CpGmeth.txt
cp Library_CpGmeth.txt ../CpGmeth_files/Library_CpGmeth.txt
echo -e `date +[%D-%R]` "\tCreated file with CpGmeth per all samples in Library" | tee -a ../Library_CpGmeth.log

#7. Methylation correction by MethylCal
#https://github.com/lb664/MethylCal