The tumor-normal paired data of pancreatic cancer whole-exome sequencing (WES) was obtained from National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) data portal (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA386668). To prepare these raw datasets (FASTQ files) into the analysis-ready training and testing dataset, which can be found in the training_data folder and tumor-only_data folder, we utilized the following alignment and variant calling tools:

# 1.	WES Alignment with Bowtie2
We performed end-to-end read alignment using BOWTIE2 version 2.2.6 (http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) for the tumor reads and the paired normal reads. The reference genome file we used was GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.

1)	Map reads against reference genome:
```
bowtie2 --end-to-end --very-fast --rg-id [ID FOR THE PAIRED-END READS] -x GRCH38 -q -1 [INPUT FASTQ FILE PAIR 1] -2 [INPUT FASTQ FILE PAIR 2] | samtools view - -Sb -h -t GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai -o [OUTPUT BAM FILE]
```

2)	Sort the output BAM file with SAMTOOLS:
```
samtools sort [INPUT BAM FILE] -o [OUTPUT SORTED BAM FILE] -m 8000000000
```

3)	Remove PCR duplicates with SAMTOOLS:
```
samtools rmdup [INPUT SORTED BAM FILE] [OUTPUT SORTED DEDUPLICATED BAM FILE]
```
The sorted and deduplicated BAM files were then used as the input for the variant calling process using three different variant calling tools, which were Mutect2, HaplotypeCaller, and BCFtools.

# 2.	Variant Calling
## 2.1.	Mutect2 
Mutect2 version 4.1.9.0 (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) was run in the tumor-normal matched mode to retrieve the somatic mutation variants calls from the tumor-normal paired data. This tumor-normal matched mode required a Panel of Normals (PON). PON was generated from the six normal BAM files by first running Mutect2 in the tumor-only mode for each normal sample and then combine the normal calls using the function CreateSomaticPanelOfNormals from Mutect2. Finally, we ran Mutect2 for each tumor sample with its matching normal sample using the generated PON.
1)	Run Mutect2 w/matched normal for the benchmark set:
```
gatk --java-options "-Xmx8g" Mutect2 -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -I [INPUT SORTED DEDUPLICATED TUMOR BAM FILE] -I [INPUT SORTED DEDUPLICATED NORMAL BAM FILE] -tumor [ID TUMOR BAM FILE] -normal [ID NORMAL BAM FILE] -pon [PON VCF.gz FILE] --germline-resource somatic-hg38_af-only-gnomad.hg38.vcf --af-of-alleles-not-in-resource 0.0000025 -L Homo_sapiens_assembly38_exome.targets.interval_list -O [OUTPUT VCF FILE]
```

2)	Run FilterMutectCalls to filter somatic variants, germline variants, and artifacts in the Mutect2 VCF callset:
```
gatk FilterMutectCalls -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -V [INPUT MUTECT2 UNFILTERED VCF] -O [OUTPUT MUTECT2 FILTERED VCF]
```

3)	Filter out the indels from the Mutect2 filtered VCFs callset.
4)	Annotate each of the SNP-only VCFs with ANNOVAR to acquire the functional prediction features:
```
perl table_annovar.pl [INPUT MUTECT2 FILTERED VCF] humandb/ -buildver hg38 -out [OUTPUT ANNOTATED VCF] --remove --protocol refGene,exac03,avsnp150,dbnsfp33a,gnomad_exome,cosmic92_coding,clinvar_20210123 --operation gx,f,f,f,f,f,f -nastring . -vcfinput -polish -xref example/gene_fullxref.txt
```

5)	Variant collapsing

The multiple VCF files, which were annotated with functional prediction features, will be combined using variant collapsing.  According to the study by Kalatskaya et al. [1], the variants identified as somatic or germline by Mutect2 can be assumed to have the same position of genome and substitution pattern. With variant collapsing, we created a dataset of unique variants.

## 2.2.	HaplotypeCaller 
Utilizing HaplotypeCaller version 4.1.9.0, we want to acquire the variants statistical features of FS, MQ, MQRankSum, QD, and ReadRankSum which we will integrate to the dataset of variants called by Mutect2. To acquire these features, we ran HaplotypeCaller on each tumor and normal sample’s BAM file using the GVCF (Genomic VCF) workflow for scalable variant calling in exome sequence data.
1)	Run the HaplotypeCaller on each tumor and normal samples BAM files to create single-sample gVCFs, with the option --emitRefConfidence GVCF, and using the .g.vcf extension for the output file:
```
gatk --java-options "-Xmx4g" HaplotypeCaller -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -I [INPUT SORTED DEDUPLICATED TUMOR BAM FILE] -O [OUTPUT .g.vcf] -A StrandBiasBySample -ERC GVCF
```
2)	Aggregate the multiple GVCF files:
```
gatk --java-options "-Xmx96g -Xms96g" CombineGVCFs -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -V [INPUT .g.vcf] -V [INPUT .g.vcf] -V [INPUT .g.vcf] -V [INPUT .g.vcf] -O [OUTPUT FILE COHORT .g.vcf]
```

3)	Joint genotyping:
```
gatk --java-options "-Xmx4g" GenotypeGVCFs -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -V [INPUT FILE COHORT .g.vcf] -O [OUTPUT FINAL COHORT VCF]
```

4)	Subset to SNPs-only callset with SelectVariants:
```
gatk SelectVariants -V [INPUT FINAL COHORT VCF] -select-type SNP -O [OUTPUT SNP-ONLY VCF]
```

5)	Hard-filtering variant:
```
gatk VariantFiltration -V [INPUT SNP-ONLY VCF] -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -O [OUTPUT FILTERED SNP-ONLY VCF]
```

## 2.3.	BCFtools
Variant calling with BCFtools version 1.11
1)	Create a list of bams to use:
```
ls *.bam > [OUTPUT BAMLIST .txt]
```

2)	Pile the multiple samples, call variants according to the targeted regions, and pipe it to bcftools to create a VCF file:
```
bcftools mpileup -d 250 -R [INPUT TARGETED REGIONS BED FILE] -B -Ou -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -b [INPUT BAMLIST .txt] | bcftools call -mv -O v -o [OUTPUT VCF]
```

3)	Filter query for the variants calling results:
```
bcftools filter -sLowQual -g3 -G10 -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15)' [INPUT VCF] > [OUTPUT FILTERED VCF]
```

# 3.	Dataset Preprocessing
## 3.1.	Variants Dataset Integration
The three variants dataset from Mutect2, HaplotypeCaller, and BCFtools were integrated based on the chromosome number (CHROM), the reference position (POS), the start position of the variants (START), the end position of the variants (END), the reference allele (REF), and alternative allele (ALT). Then, we only included the variants labeled as “PASS”, or “germline” by Mutect2, because these labels were taken into consideration for the variants pre-labeling step.

## 3.2.	Imputation of Missing Prediction Features
The missing prediction features were imputed using the k-Nearest Neighbors approach, with k = 41.

## 3.3.	Variant Pre-Labeling
The variants labeled as “PASS” by the Mutect2 filter and labeled as “confirmed somatic variants” by COSMIC, were assigned the primary label “somatic mutation”. The variants “unknown” to COSMIC but labeled as “germline” by Mutect2 and known to dbSNP common (MAF>=1%), were labeled as “germline”.  
