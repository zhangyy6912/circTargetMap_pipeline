hisat2 -x GRCh37.p13.genome.clean.fa -U read1_torRNA_Unmapped_really.fq -S read1.toGenome.hisat2.sam --un ./ --no-softclip --phred33 -p 8 2> read1.toGenome.hisat2.log
sleep 10
hisat2 -x GRCh37.p13.genome.clean.fa -U read2_torRNA_Unmapped_really.fq -S read2.toGenome.hisat2.sam --un ./ --no-softclip --phred33 -p 8 2> read2.toGenome.hisat2.log
sleep 10
perl remove_both_mapping_reads.pl read1_torRNA_Unmapped_really.fq read2_torRNA_Unmapped_really.fq read1.toGenome.hisat2.sam read2.toGenome.hisat2.sam read1.toGenome.unmapped.any.byHisat2.fq read2.toGenome.unmapped.any.byHisat2.fq
bwa mem -t 16 -k 12 -T 15 -o read1.toGenome.bwa.sam GRCh37.p13.genome.clean.fa read1.toGenome.unmapped.any.byHisat2.fq
bwa mem -t 16 -k 12 -T 15 -o read2.toGenome.bwa.sam GRCh37.p13.genome.clean.fa read2.toGenome.unmapped.any.byHisat2.fq

