#https://github.com/zhangrengang/Centromics
centromics -l hifi.fq.gz -g ref.fa -pre hifi

#dyad symmetries
#ref:https://doi.org/10.1016/j.cell.2022.06.045
#use EMBOSS
./palindrome -sequence EroSat3.fasta -minpallen 4 -maxpallen 100 -gaplimit 20 -nummismatches 0 -overlap -outfile EroSat3.dyad