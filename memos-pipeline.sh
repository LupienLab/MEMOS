#!/bin/bash
#$ -N MEMOS 
#$ -j y 
#$ -cwd 
#$ -q lupiengroup   

# 1=TF
# 2=FASTA for regions to be scanned   
# 3=Mutation file (vcf)  
# 4=Number of permutations for perm test 
# 5=Output directory
# 6=Pvalue (0.0001 recommended)
# 7=Flanking region length
# 8=Background regions for permutation test (if not provided the whole genome is the default)

WDIR="/mnt/work1/users/lupiengroup/CommonScripts/MEMOS"
TFDIR=/mnt/work1/users/lupiengroup/CommonScripts/MEMOS-wrapper/TF-pfms
HG19=$WDIR/hg19.genome
###
module load python/2.7
module load R/3.2.2
###
TF=$1
MUTATION=$3
NUM=$4
OUTDIR=$5
PVAL=$6
FLANK=$7
FASTA=$OUTDIR/tmp.fa
###
mkdir -p $OUTDIR
###
echo "creating a bed file with required columns from the mutation vcf file ..."
grep -F '#' $MUTATION > $OUTDIR/vcfheader.tmp
grep -Fv '#' $MUTATION > $OUTDIR/vcfbody.tmp
awk -F "\t" '{print "chr" $1 "\t" $2 "\t" $2 "\t" $4 ">" $5}' $OUTDIR/vcfbody.tmp > $OUTDIR/mutations.bed
MUT=$OUTDIR/mutations.bed

echo "reformat the fasta file; creating tmp.fa ..."
cp $2 "${OUTDIR}/tmp.fa"
$WDIR/formatFA.sh $FASTA
#FASTA=/mnt/work1/users/lupiengroup/People/parisa/MOODSanalysis/20PCa_H3K27acpeaks_extended500bp_merged.fa
############################# MOODS #############################
echo "started running moods for $TF ..."
for pfm in $TFDIR/$TF/*.pfm
do
	pfmName="${pfm##*/}"
	echo "motif name = $pfmName ..."    
	echo "running MOODS with pvalue = $PVAL ..."
	echo "parameters:  -m $pfm -s $FASTA -p $PVAL > $OUTDIR/${pfmName%.*}-${PVAL}.out "
	python /mnt/work1/software/MOODS/1.9.2/scripts/moods_dna.py -m $pfm -s $FASTA -p $PVAL > $OUTDIR/${pfmName%.*}-${PVAL}.out
	n="$(awk -F' ' '{print NF; exit}' $pfm)" 
	n=$(($n-1)) 
	awk -v x=$n -F '[, ]' '{print $1"\t"$2+$5"\t"$2+$5+x"\t"$6 "\t" $7 "\t" $8}' $OUTDIR/${pfmName%.*}-${PVAL}.out > $OUTDIR/${pfmName%.*}-${PVAL}.bed

####################### Permutation Test #########################
	echo "started running permutation test for $TF ..." 
	mkdir $OUTDIR/${TF}_permutation_${PVAL}
	if [ -z "$8" ]; then
		$WDIR/permutationTest.sh $OUTDIR/${TF}_permutation_${PVAL} $OUTDIR/${pfmName%.*}-${PVAL}.bed $FLANK $MUT $HG19 $NUM	
	else
		 $WDIR/permutationTest.sh $OUTDIR/${TF}_permutation_${PVAL} $OUTDIR/${pfmName%.*}-${PVAL}.bed $FLANK $MUT $HG19 $NUM $8
	fi
	echo "done with the permutation test!"
###################### Creating Plots ############################
	# centerMotifs centers the motifs and returns a file with  postfix .center
	$WDIR/centerMotifs.sh $OUTDIR/${pfmName%.*}-${PVAL}.bed
	Rscript $WDIR/MutationEnrichment_stream.R $OUTDIR $MUT $OUTDIR/${pfmName%.*}-${PVAL}.bed.centers $FLANK $FLANK ${TF}_mutation_enrichment_${FLANK} 1200 1000 $TF 10
	
done 
