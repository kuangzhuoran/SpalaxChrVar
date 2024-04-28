#1.Processing of reference genomes
bwa index ref.fa
gatk CreateSequenceDictionary -R ref.fa -O ref.dict

#2.mapping
#Take the sample 'SBN2' as an example.
bwa mem -t 10 -R '@RG\tID:SBN2\tPL:illumina\tLB:SBN2\tSM:SBN2' ref.fa 1.fq 2.fq | samtools view --threads 20 -Sb - | samtools sort -@ 36 -m 2G -O bam - -o SBN2.bam
gatk MarkDuplicates -I SBN2.bam  -O SBN2.markdup.bam -M SBN2.markdup_metrics.txt
samtools index SBN2.markdup.bam

#3.variants detection
gatk HaplotypeCaller --java-options -Xmx30G -R ref.fa --emit-ref-confidence GVCF -I SBN2.markdup.bam -O SBN2.gvcf
gatk CombineGVCFs -R ref.fa --variant SBN2.gvcf --variant SBN3.gvcf --variant SBN4.gvcf ...... -O hebing.gvcf
gatk GenotypeGVCFs -R ref.fa -V hebing.gvcf  -O hebing.vcf
bgzip -f hebing.vcf
tabix -p vcf hebing.vcf.gz

#4.filter
gatk SelectVariants -select-type SNP -V hebing.vcf.gz -O hebing.snp.vcf.gz
gatk VariantFiltration -V hebing.snp.vcf.gz --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O hebing.snp.hardfilter.vcf.gz
vcftools --remove-indels --recode --recode-INFO-all --gzvcf hebing.snp.filter.vcf.gz --maf 0.05 --minDP 5 --maxDP 50 --minGQ 20 --hwe 0.001 --max-missing 0.8 --stdout  > raw1.snp.vcf
bcftools filter -g 5 -O v raw1.snp.vcf -o raw2.snp.vcf

grep -v 'Filter' raw2.snp.vcf > final.SNP.vcf


#5.Imputation by beagle
#Optional, part of the analysis requires the use of beagle to process vcf files
vcftools --vcf final.SNP.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out 2allele
java -jar -Xmn12G -Xms24G -Xmx48G beagle.r1399.jar gt=2allele.recode.vcf out=final.SNP.phased.vcf