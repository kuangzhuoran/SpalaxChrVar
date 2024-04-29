#calculation of Fst
#final.SNP.vcf
bgzip final.SNP.vcf
tabix -p vcf final.SNP.vcf.gz
pixy --stats fst --bypass_invariant_check 'yes' --vcf final.SNP.vcf.gz --populations gal.gol.pop --window_size 50000 --n_cores 10 --output_prefix gal.gol


#########################################


#estimation of gene flow (RFmix)
#java -jar -Xmn12G -Xms24G -Xmx48G beagle.r1399.jar gt=2allele.recode.vcf out=final.SNP.phased
bgzip final.SNP.phased.vcf
tabix -p vcf final.SNP.phased.vcf.gz
#split vcf
for i in {1..29}
do
vcftools --gzvcf final.SNP.phased.vcf.gz --chr Chr${i} --recode --recode-INFO-all --stdout > Chr${i}.vcf
vcftools --vcf Chr${i}.vcf --keep ref.pop --recode --recode-INFO-all --out Chr${i}.ref.vcf
vcftools --vcf Chr${i}.vcf --keep target.pop --recode --recode-INFO-all --out Chr${i}.target.vcf
#
done
#get genetic map
for i in {1..29}
do
plink --vcf Chr${i}.vcf --recode --out Chr${i} --allow-extra-chr
awk '{print $1"\t"$4"\t"$4*0.000554779412}' Chr${i}.map  | sort -Vk 1 | awk '{print "Chr"$1"\t"$2"\t"$3}' > Chr${i}.genetic.map
done
#calculate RFmix
#https://github.com/slowkoni/rfmix
mkdir output
for i in {1..29}
do
rfmix -f Chr${i}.target.vcf -r Chr${i}.ref.vcf -g Chr${i}.genetic.map -m ref.pop -o output/Chr${i} --chromosome=Chr${i}
done


#########################################


#calculate XP-EHH
vcftools --vcf final.SNP.phased.vcf --keep carmeli.pop --recode --recode-INFO-all --stdout > carmeli.vcf
vcftools --vcf final.SNP.phased.vcf --keep judaei.pop  --recode --recode-INFO-all --stdout >  judaei.vcf

for i in {1..29}
do
vcftools --vcf carmeli.vcf --chr Chr${i} --recode --recode-INFO-all --stdout > carmeli.split/carmeli.Chr${i}.vcf
vcftools --vcf judaei.vcf --chr Chr${i} --recode --recode-INFO-all --stdout > judaei.split/judaei.Chr${i}.vcf
done

#get genetic map
for i in {1..29}
do
plink --vcf Chr${i}.vcf --recode --out Chr${i} --allow-extra-chr
awk '{print $1"\t"$4"\t"$4*0.000554779412}' Chr${i}.map  | sort -Vk 1 | awk '{print "Chr"$1"\t"$2"\t"$3}' > Chr${i}.genetic.map
done

for i in {1..29}
do
./selscan --xpehh --vcf carmeli.split/carmeli.Chr${i}.vcf --vcf-ref judaei.split/judaei.Chr${i}.vcf --map Chr${i}.genetic.map --out xpehh/xpehh.Chr${i}.out
#carmeli is the species which chromosome fusion has occurred
done
cd xpehh
./norm --xpehh --files Chr*.out.xpehh.out --bp-win --winsize 50000

for i in {1..29}
do
awk '{print "Chr"'${i}'"\t"$1"\t"$2"\t"$5}' Chr${i}.norm.50kb.windows | awk '{if($4==-1)print $1"\t"$2"\t"$3"\t""0"; else print $1"\t"$2"\t"$3"\t"$4}' > Chr${i}.norm.50kb.judaei
done
cat Chr*.norm.50kb.judaei > all.norm.50kb.judaei && rm Chr*.norm.50kb.judaei

for i in {1..29}
do
awk '{print "Chr"'${i}'"\t"$1"\t"$2"\t"$4}' Chr${i}.norm.50kb.windows | awk '{if($4==-1)print $1"\t"$2"\t"$3"\t""0"; else print $1"\t"$2"\t"$3"\t"$4}' > Chr${i}.norm.50kb.carmeli
done
cat Chr*.norm.50kb.carmeli > all.norm.50kb.carmeli && rm Chr*.norm.50kb.carmeli


#########################################


#Extraction of SNPs for specified species
gatk CombineGVCFs -R ref/carmeli.fa --variant EC-10A.gvcf --variant EC-1A.gvcf --variant EC-2A.gvcf --variant EC-3A.gvcf --variant EC-4B.gvcf --variant EC-5B.gvcf --variant EC-6A.gvcf --variant EC-7A.gvcf --variant EC-8A.gvcf --variant EC-9A.gvcf -O only.carmeli.gvcf 
gatk GenotypeGVCFs -R ref/carmeli.fa -V ./only.carmeli.gvcf -O ./only.carmeli.vcf
bgzip -f ./only.carmeli.vcf
tabix -p vcf ./only.carmeli.vcf.gz

gatk SelectVariants -select-type SNP -V only.carmeli.vcf.gz -O ./only.carmeli.snp.vcf.gz
gatk VariantFiltration -V ./only.carmeli.snp.vcf.gz --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O ./only.carmeli.snp.HardFilter.vcf.gz && mkdir check/onlycarmeli.HardFilter.done

vcftools --remove-indels --recode --recode-INFO-all --gzvcf ./only.carmeli.snp.HardFilter.vcf.gz  --maf 0.05 --minDP 5 --maxDP 50 --minGQ 20 --hwe 0.001 --max-missing 0.8 --stdout > ./only.carmeli.snp.TiaojianFilter.vcf.gz

bcftools filter -g 5 -O v -o only.carmeli.final.vcf ./only.carmeli.snp.TiaojianFilter.vcf.gz


#########################################


#ROH analysis
plink --vcf only.carmeli.final.vcf --allow-extra-chr --chr-set 29 --out carmeli
plink --bfile carmeli -allow-extra-chr --chr-set 29 --homozyg-window-snp 1000 --homozyg-window-het 10 --homozyg-window-threshold 0.001  --homozyg-kb 100 --out carmeli.ROH


#########################################


#genetic load analysis
#Use of dN/dS to assess genetic load
java -Xmx30g -jar snpEff.jar -v only.carmeli.final.vcf > only.car.anno.vcf

grep -v 'MODIFIER' only.car.anno.vcf > carmeli.anno.have.impact.vcf

grep -v '#' carmeli.anno.have.impact.vcf > carmeli.anno.have.impact.nohead.vcf

grep 'synonymous_variant' carmeli.anno.have.impact.nohead.vcf > carmeli.anno.dS.nohead.vcf
grep 'missense_variant' carmeli.anno.have.impact.nohead.vcf > carmeli.anno.dN.nohead.vcf
grep -v 'synonymous_variant' carmeli.anno.have.impact.nohead.vcf > carmeli.anno.youhai.nohead.vcf

bioawk -c fastx '{print $name, length($seq)}' carmeli.fa > carmeli.chr.length
bedtools makewindows -g carmeli.chr.length -w 500000 > 500k.windows

awk '{print $1"\t"$2-1"\t"$2}' carmeli.anno.dN.nohead.vcf > carmeli.anno.dN.bed
awk '{print $1"\t"$2-1"\t"$2}' carmeli.anno.dS.nohead.vcf > carmeli.anno.dS.bed
#then counted the number of SNPs per chromosome that were annotated as synonymous (dS.bed) and non-synonymous (dN.bed) mutations
