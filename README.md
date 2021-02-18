# ImprintSeq
ImprintSeq is a hybridization-based custom-sequencing panel to interrogate 63 imprinted DMRs (iDMRs). This panel uses 3,989 probes to cover 203,279 base pairs across 63 genomic regions involved in genomic imprinting. Probes were designed for both strands of bisulfite converted DNA. The design interrogates over 9,257 CpG sites. ImprintSeq is a hybridization-based custom-sequencing panel to interrogate 63 imprinted DMRs (iDMRs). This panel uses 3,989 probes to cover 203,279 base pairs across 63 genomic regions involved in genomic imprinting. Probes were designed for both strands of bisulfite converted DNA. The design interrogates over 9,257 CpG sites. 
Bisulfite sequencing libraries were pooled, normalized to 4nM and sequenced in a lane of Illumina HiSeq4000 sequencing machine (Illumina, San Diego, CA) using 10% PhiX. The BCL files produced by the Hiseq4000 instrument were demultiplexed into fastq files using Illumina's bcl2fastq v2.19 (Illumina, San Diego, CA) with all the default settings (including adaptor trimming). After demultiplexing, the quality of the sequencing data was assessed by FastQC software. Then adapter removal and quality filtering were performed with Cutadapt1 (--m 30 --q 30) and 10 bp on 5’ end were removed in an extra trimming step using NuGen in-house Python script. After re-evaluating quality by FastQC, data was aligned with human build hg19 using Bismark2 software (--pbat option), PCR duplicates were removed with deduplicate_bismark option and finally, methylation information was extracted with bismark_methylation_extractor tool. After the extraction of methylation with Bismark, we filtered out from the analysis any CpG with less than 100 reads of coverage to achieve accurate DNA methylation measurements.
