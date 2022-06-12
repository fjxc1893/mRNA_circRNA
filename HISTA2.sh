while read sample;do
hisat2 -x genome/genome  -1 clean_data/$sample.R1.fq.gz -2 clean_data/$sample.R2.fq.gz -p 10 --rna-strandness RF --fr -S  sam/$sample.sam  2>$sample.hisat2.log
samtools view -uS  sam/$sample.sam | samtools sort -@ 5 -o  sam/$sample.bam 
samtools index sam/$sample.bam
done<sample.list