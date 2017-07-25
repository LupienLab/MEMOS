#!/bin/bash
#$ -q light.q
#$ -N MOODS

# 1=PFM for TF
# 2=FASTA for regions to be scanned
# 3=output file name
module load python/2.7
python /mnt/work1/software/MOODS/1.9.2/scripts/moods_dna.py -m $1 -s $2 -p 0.0001 > $3
n="$(awk -F' ' '{print NF; exit}' $1)"
n=$(($n-1))
awk -v x=$n -F '[, ]' '{print $1"\t"$2+$5"\t"$2+$5+x"\t"$6}' $3 > $3.bed
