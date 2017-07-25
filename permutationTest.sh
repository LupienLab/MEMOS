#!/bin/bash
#$ -q light.q
#$ -N permTest

# Conducts permutation test by shuffling motif + flanking regions across open chromatin n times & counting overlapping mutations each time
# The proportion of permutations that had more mutations than the actual motifs + flanking regions

# Parameters
# $1=output directory
# $2=motif regions (BED)
# $3=number of flanking bp
# $4=mutations (BED)
# $5=genome file with list of chromosomes & their sizes
# $6=number of permutations

module load bedtools/2.23.0

cd $1
# extend motifs by the specified flanking size
bedtools slop -i $2 -g $5 -b $3 > $1/flank.bed
# count mutations in motifs + flanking regions
mutInMotifs=$(bedtools intersect -wa -a $4 -b $1/flank.bed | sort -u | wc -l)
echo $mutInMotifs > $1/results.txt
count=0
n=$6
# do n permutations
for i in `seq 1 $n`;
do
	if [ -z "$7" ]; then # randomly shuffle motif + flanking regions across genome
        	bedtools shuffle -i $2 -g $5 -excl $2 > $1/shuffle.bed
	else # randomly shuffle motif + flanking regions across specified region
		bedtools shuffle -i $2 -g $5 -excl $2 -incl $7 > $1/shuffle.bed
	fi
	bedtools slop -i $1/shuffle.bed -g $5 -b $3 > $1/shuffleFlank.bed
	# count mutations in those regions
	mutInPerm=$(bedtools intersect -wa -a $4 -b $1/shuffleFlank.bed | sort -u | wc -l)
	# add 1 to count if there's more mutations in the permutated motifs + flanking than the actual motifs + flanking regions
	if [ "$mutInPerm" -gt "$mutInMotifs" ]; then
		count=$((count+1))
        fi
	echo $mutInPerm >> $1/results.txt
done
# calculate p-value
pValue=$(awk "BEGIN {printf \"%.3f\",${count}/${n}}")
echo $pValue >> $1/results.txt
