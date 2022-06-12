# mRNA_circRNA

Analysis of RNA-seq data (mRNA and circRNA) from a rat lumbosacral spinal root avulsion model. Data files include:
# 1. Data Preparation
The data used for analysis can be downloaded from GSE203053, names of the datasets are showed below:
LSRA1.R1.fq; LSRA1.R2.fq; LSRA2.R1.fq; LSRA2.R2.fq; LSRA3.R1.fq; LSRA3.R2.fq
Sham1.R1.fq; Sham1.R2.fq; Sham2.R1.fq; Sham2.R2.fq; Sham3.R1.fq; Sham3.R2.fq
# 2. Preparation of Analysis Software, including installing Miniconda3, setting up environment to run the python and R
# 3. Analysis Procedure, including:
## For mRNA analysis, code for trimmomatic(QC.sh), HISAT2(HISTA2.sh), HTSeq-count(HTseq-count.sh)
## For circRNA analysis, code for BWA and CIRI(step1-BWA+CIRI2-deal, step2-CIRI2_merged.py)
