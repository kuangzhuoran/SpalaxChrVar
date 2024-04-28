#1.Jcvi
gffread carmeli.gff -g carmeli.fa -x carmeli.cds
python -m jcvi.formats.gff bed --type=mRNA --key=ID carmeli.gff -o carmeli.bed
gffread judaei.gff -g judaei.fa -x judaei.cds
python -m jcvi.formats.gff bed --type=mRNA --key=ID judaei.gff -o judaei.bed

python -m jcvi.compara.catalog ortholog --no_strip_names carmeli judaei



#2.NGenomeSyn
nucmer  --mum --mincluster 500 -t 30 galili.fa golanii.fa -p galili_golani && 
delta-filter -1 -i 90 -l 150000 galili_golani.delta > golani_carmeli.filter.delta  && 
show-coords -c -r galili_golani.filter.delta > galili_golani.filter.coords 
perl coords.to.link.pl galili_golani.filter.coords > galili_golani.filter.link

nucmer  --mum --mincluster 500 -t 30 golani.fa carmeli.fa -p golani_carmeli && 
delta-filter -1 -i 90 -l 150000 golani_carmeli.delta >  golani_carmeli.filter.delta  && 
show-coords -c -r golani_carmeli.filter.delta > golani_carmeli.filter.coords 
perl coords.to.link.pl golani_carmeli.filter.coords > golani_carmeli.filter.link

nucmer  --mum --mincluster 500 -t 30 carmeli.fa judaei.fa -p carmeli_judaei &&
delta-filter -1 -i 90 -l 150000 carmeli_judaei.delta >  carmeli_judaei.filter.delta  && 
show-coords -c -r carmeli_judaei.filter.delta > carmeli_judaei.filter.coords
perl coords.to.link.pl carmeli_judaei.filter.coords > carmeli_judaei.filter.link


bioawk -c fastx '{print $name,length($seq)}' galili.fa > galili.length
bioawk -c fastx '{print $name,length($seq)}' golani.fa > golani.length
bioawk -c fastx '{print $name,length($seq)}' carmeli.fa > carmeli.length
bioawk -c fastx '{print $name,length($seq)}' judaei.fa > judaei.length

./bin/NGenomeSyn -InConf NGenomeSyn.config.txt -OutPut OUT
