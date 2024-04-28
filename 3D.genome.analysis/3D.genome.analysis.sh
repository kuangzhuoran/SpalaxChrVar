#Take 'S. galili' as an example.

#HiC-Pro pipeline
#Hic-pro v3.1.0
bowtie2-build galili.fa ref
./HiC-Pro-3.1.0/bin/utils/digest_genome.py -r ^GATC -o ref.bed galili.fa
samtools faidx galili.fa
awk '{print $1 "\t" $2 }' galili.fa.fai > galili.sizes
./HiC-Pro-3.1.0/bin/HiC-Pro  -i ./HIC/ref/rawdata -o ./HIC/ref/result -c config-hicpro.txt

#Construct matrix
awk '{print $1"\t"$2"\t"$3"\t"$4-1}' galili.40k.abs.bed > DenseMatrix/galili.40k.bed
awk '{print $1-1"\t"$2-1"\t"$3}' galili.40k.iced.matrix > DenseMatrix/galili.40k.matrix
cd DenseMatrix
python2 ./HiC-Pro_3.1.0/bin/utils/sparseToDense.py -b galili.40k.bed galili.40k.matrix --perchr

#convert matrix to h5
hicConvertFormat -m galili.40k.iced.matrix -o galili.40k.iced --bedFileHicpro galili.40k.abs.bed --inputFormat hicpro --outputFormat h5;
hicCorrectMatrix correct -m galili.40k.iced.h5 --filterThreshold -1.5 5 -o galili.40k.iced.corrected.h5;


#call TAD
for i in {1..29}
do
hicFindTADs -m ../galili.40k.iced.corrected.h5 --outPrefix galili.40k.Chr${i}.TAD -p 20 --chromosomes Chr${i} --correctForMultipleTesting fdr --minDepth 160000 --maxDepth 320000 --step 40000 --thresholdComparisons 0.05
done


#call A/B compartment
for i in {1..29}
do
python2 01.runchangematrix.insulation.py -i ./DenseMatrix/galili.40k_Chr10_dense.matrix -g galili -c Chr${i} -o 01.runchangematrix -s 40000
done

for i in {1..29}
do
perl 02.matrix2compartment.pl -i 01.runchangematrix/galili_Chr${i}_40000.insulation.matrix -o 02.matrix2compartment/Chr${i}  --et
done

awk '{if($3=="gene")print}' galili.gff > galili.gene.gff
for i in {1..29}
do
python 03.matrix2EigenVectors.py -i ./02.matrix2compartment/Chr${i}.zScore.matrix.gz  -r ./galili.gene.gff -v
done


#fithic
./HiC-Pro-3.1.0/bin/utils/hicpro2fithic.py -i galili.10k.iced.matrix -b galili.10k.abs.bed -s galili.10k.iced.matrix.biases
gunzip -c fithic.biases.gz > fithic.biases && gunzip -c fithic.fragmentMappability.gz > fithic.fragmentMappability
awk -F '\t' '!($5=="0")'  fithic.fragmentMappability > fithic.fragmentMappability1
awk NF fithic.biases | awk -F '\t' '!($3=="-1")' > fithic.biases1
wc -l fithic.biases1 fithic.fragmentMappability1 > statistics #查看行数是否一样
gzip fithic.biases1 fithic.fragmentMappability1
python ./fithic/fithic.py -f fithic.fragmentMappability1.gz -i fithic.interactionCounts.gz -t fithic.biases1.gz -o galili.10k -l galili.10k -v -x All -r 10000
cd galili.10k
python2 runfilterpvalue.qvalue.count.py galili.10k.spline_pass1.res10000.significances.txt.gz galili.10k.spline_pass1.res10000.significances.bed

awk 'BEGIN{ FS="\t";OFS="\t" }{if($1=$3)print $1"\t"$2-5000"\t"$2+5000"\t"$3"\t"$4-5000"\t"$4+5000"\t"NR"\t"$5"\t"$6}' galili.10k.spline_pass1.res10000.significances.bed | grep -v 'scaffold'| sort -Vk 1 > galili.loop.bed


#######################################
#get chain file
minimap2 -cx asm5 --cs galili.fa golani.fa > golani2galili.paf
transanno minimap2chain golani2galili.paf --output golani2galili.minimap2.chain

#compartment transformation
cat Chr*.zScore.eigen1.bedGraph | grep -v 'track' | sort -Vk 1 > galili.all.PC1.bed
cat Chr*.zScore.eigen1.bedGraph | grep -v 'track' | sort -Vk 1 > golani.all.PC1.bed

CrossMap.py bed golani2galili.chain golani.All.PC1.bed > golani.trans2.galili.PC1.bed
#然后只提取了和融合相关染色体的信息
awk '{print $6"\t"$7"\t"$8"\t"$9}' golani.trans2.galili.PC1.fusion.bed | sort -Vk 1 | awk '{if($2<$3) print $1"\t"$2"\t"$3"\t"$4; else print $1"\t"$3"\t"$2"\t"$4}' | uniq > tmp1
bedtools intersect -a galili.All.PC1.bed -b tmp1 -F 0.8 -wa -wb > tmp2

#awk '{print $5"\t"$6"\t"$7}' tmp2 | uniq | bedtools sort -i - | bedtools merge -i - > map.region.bed #All Map Length
#awk '{print $3-$2+1}' map.region.bed | awk '{sum+=$1} END {print "All Map Length = ", sum}' #All Map Length

awk '{if($4*$8>0)print}' tmp2 | awk '{print $5"\t"$6"\t"$7"\t"$8}' | uniq | bedtools sort -i - | bedtools merge -i - > stable.region.bed
awk '{print $3-$2+1}' stable.region.bed | awk '{sum+=$1} END {print "Stable Length = ", sum}'  #stable Length

awk '{if($4*$8>0)print}' tmp2 | awk '{if($4>0)print}' | awk '{print $5"\t"$6"\t"$7"\t"$8}' | uniq | bedtools sort -i - | bedtools merge -i - > stableA.region.bed
awk '{print $3-$2+1}' stableA.region.bed | awk '{sum+=$1} END {print "StableA Length = ", sum}'  #stableA Length

awk '{if($4*$8>0)print}' tmp2 | awk '{if($4<0)print}' | awk '{print $5"\t"$6"\t"$7"\t"$8}' | uniq | bedtools sort -i - | bedtools merge -i - > stableB.region.bed
awk '{print $3-$2+1}' stableB.region.bed | awk '{sum+=$1} END {print "StableB Length = ", sum}'  #stableB Length

awk '{if($4*$8<0)print}' tmp2 | awk '{if($4>0)print}' > BtoA.tmp
awk '{print $5"\t"$6"\t"$7}' BtoA.tmp | uniq | bedtools sort -i - | bedtools merge -i - > BtoA.region.bed
awk '{print $3-$2+1}' BtoA.region.bed | awk '{sum+=$1} END {print "B to A = ", sum}'

awk '{if($4*$8<0)print}' tmp2 | awk '{if($4<0)print}' > AtoB.tmp
awk '{print $5"\t"$6"\t"$7}' AtoB.tmp | uniq | bedtools sort -i - | bedtools merge -i - > AtoB.region.bed
awk '{print $3-$2+1}' AtoB.region.bed | awk '{sum+=$1} END {print "A to B = ", sum}'


#shared TADs and boundaries
awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$2"\t"$3"\t"$3}' golani.All.40k.TAD.bed | grep -v 'scaffold' > golani.All.TAD.liftover
wc -l golani.All.40k.TAD.bed
liftOver golani.All.TAD.liftover golani2galili.minimap2.chain golani2galili.body.mapped body.unmapped -minMatch=0.25
grep -v 'scaffold' golani2galili.body.mapped | wc -l
bedtools intersect -a galili.All.40k.TAD.bed -b golani2galili.body.mapped -F 0.8 -f 0.8 -wo | grep -v 'scaffold' > tmp2
wc -l tmp2


awk '{print $1"\t"$2"\t"$3"\t"$1"\t"$2"\t"$3"\t"$3}' golani.All.40k.TADboundaries.bed | grep -v 'scaffold' > golani.All.40k.TADboundaries.liftover 
wc -l golani.All.40k.TADboundaries.bed
liftOver golani.All.40k.TADboundaries.liftover golani2galili.minimap2.chain golani2galili.boundary.mapped boundary.unmapped -minMatch=0.33
grep -v 'scaffold' golani2galili.boundary.mapped | wc -l
bedtools intersect -a golani.All.40k.TADboundaries.bed -b golani2galili.boundary.mapped -wo | grep -v 'scaffold' > tmp2
wc -l tmp2


#Conserved loop
#https://github.com/adadiehl/mapLoopLoci
sed 's/Chr/chrChr/g' golani2galili.minimap2.chain > golani2galili.rename.chain
python2 ./mapLoopLoci.py golani.loop.bed galili.loop.bed golani2galili.rename.chain > golani2galili.out