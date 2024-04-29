biser -o golani --gc-heap 1G --threads 40 golani.softmask.fa
cut -f2,3,4 golani.biser.elem.txt | bedtools sort -i | bedtools merge -i - > golani.biser.elem.merge
awk '{print $3-$2}' golani.biser.elem.merge | awk '{sum += $1};END {print sum}' 