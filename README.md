# MEMOS

**M**utation **E**nrichment of **Mo**tif**s**

Developed by [Parisa Mazrooei](https://github.com/mazrooei) and [Tahmid Mehdi](https://github.com/tahmidmehdi)
MEMOS has been used for the analysis of [Mazrooei, Parisa, et al. "Cistrome Partitioning Reveals Convergence of Somatic Mutations and Risk Variants on Master Transcription Regulators in Primary Prostate Tumors." Cancer cell 36.6 (2019): 674-689.](https://www.sciencedirect.com/science/article/abs/pii/S1535610819304799)

# Description

MEMOS computes the enrichment of mutations within motifs and flanking regions of your transcription factor of interest.

Location on the cluster: `/mnt/work1/users/lupiengroup/CommonScripts/MEMOS-wrapper/`

# Installation

Download the repo via `git clone`.

# Usage

```shell
sh memos-pipeline.sh TF FASTA MUT N_PERM OUTDIR PVAL FLANK BED
```

## Parameters

| Parameter | Description |
|-----------|-------------|
| `TF` | Name of your transcription factor of interest. This name will be used to lookup the motif for that transcription factor |
| `FASTA` | FASTA file containing the regions to be scanned |
| `MUT` | VCF file of mutations |
| `N_PERM` | Numer of permutations for permutation test |
| `OUTDIR` | Output directory |
| `PVAL` | P-value for permutation test. 0.001 recommended |
| `FLANK` | Length of flanking region. This will be considered for both left and right flanking regions |
| `BED` | **Optional**. Background regions for permutation test. If not provided the whole genome is the default |

## Notes

1. All paths given should be absolute paths
1. At the moment, VCFs are considered to match [1000GenomeProject](http://www.internationalgenome.org/data) format which means Chromosomes are considered without `chr` at the beginning of the chromosome number
