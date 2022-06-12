while read sample;do
htseq-count -i gene_id  -f bam -s reverse -r name  sam/$sample.bam  genome/gene.gtf > counts/$sample.counts.tx
done<sample.list