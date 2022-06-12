while read sample;do bwa mem -T 19 -t 8 genome.fa clean_data/$sample.R1.fq.gz clean_data/$sample.R2.fq.gz > CIRI/$sample.align.sam;done<sample.list

while read sample;do CIRI2.pl -T 12 -F genome/genome.fa -I CIRI/$sample.align.sam -O CIRI/$sample.CIRI.gtf;done<sample.list