# read1
perl step1.print_read_blocks.pl read1.toGenome.bwa.sam 1 > out1.read1.toGenome.bwa.blocks
#for BSJ and non-BSJ
perl step2.find_circ_Link_BSJ-supported.pl out1.read1.toGenome.bwa.blocks result2.circBase.spliceSite_circRNA.good.bedpair > out2.read1.circRNA_linked.blocks
perl step2.find_circ_Link_nonBSJ-supported.pl out1.read1.toGenome.bwa.blocks bed_file_of_highly_circularized_circRNAs >> out2.read1.circRNA_linked.blocks

perl step3.remove_potential_intraPairing.pl out2.read1.circRNA_linked.blocks > out3.read1.circRNA_linked.targets 2> out4.read1.circRNA_linked.bed
perl step4.reMap_target_loci.pl out4.read1.circRNA_linked.bed *read1_torRNA_Unmapped_really.fq.gz > out5.read1.target_block.readSeq.fq
bowtie2 --sensitive --end-to-end -p 12 -N 1 -L 17 -x GRCh37.p13.genome.clean.fa -U out5.read1.target_block.readSeq.fq -S out6.read1.target_block.readSeq.reMap.bowtie2.sam 2> out6.read1.target_block.readSeq.reMap.bowtie2.log
hisat2 -p 12 -x GRCh37.p13.genome.clean.fa -U out5.read1.target_block.readSeq.fq -S out6.read1.target_block.readSeq.reMap.hisat2.sam 2> out6.read1.target_block.readSeq.reMap.hisat2.log
perl fq2fa.pl out5.read1.target_block.readSeq.fq > out5.read1.target_block.readSeq.fa
blastn -query out5.read1.target_block.readSeq.fa -task megablast -db human_rRNA.fa -out out6.read1.target_block.readSeq.blast_rRNA.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
perl step5.filter_low_MAPQ.pl out3.read1.circRNA_linked.targets out4.read1.circRNA_linked.bed out6.read1.target_block.readSeq.reMap.bowtie2.sam out6.read1.target_block.readSeq.reMap.hisat2.sam out6.read1.target_block.readSeq.blast_rRNA.out out7.read1.circRNA_linked.targets 2> log.filtered.read1.blocks

# read2
perl step1.print_read_blocks.pl read2.toGenome.bwa.sam 2 > out1.read2.toGenome.bwa.blocks

perl step2.find_circ_Link_BSJ-supported.pl out1.read2.toGenome.bwa.blocks result2.circBase.spliceSite_circRNA.good.bedpair > out2.read1.circRNA_linked.blocks
perl step2.find_circ_Link_nonBSJ-supported.pl out1.read2.toGenome.bwa.blocks bed_file_of_highly_circularized_circRNAs >> out2.read2.circRNA_linked.blocks

perl step3.remove_potential_intraPairing.pl out2.read2.circRNA_linked.blocks > out3.read2.circRNA_linked.targets 2> out4.read2.circRNA_linked.bed
perl step4.reMap_target_loci.pl out4.read2.circRNA_linked.bed *read2_torRNA_Unmapped_really.fq.gz > out5.read2.target_block.readSeq.fq
bowtie2 --sensitive --end-to-end -p 12 -N 1 -L 17 -x GRCh37.p13.genome.clean.fa -U out5.read2.target_block.readSeq.fq -S out6.read2.target_block.readSeq.reMap.bowtie2.sam 2> out6.read2.target_block.readSeq.reMap.bowtie2.log
hisat2 -p 12 -x GRCh37.p13.genome.clean.fa -U out5.read2.target_block.readSeq.fq -S out6.read2.target_block.readSeq.reMap.hisat2.sam 2> out6.read2.target_block.readSeq.reMap.hisat2.log
perl fq2fa.pl out5.read2.target_block.readSeq.fq > out5.read2.target_block.readSeq.fa
blastn -query out5.read2.target_block.readSeq.fa -task megablast -db human_rRNA.fa -out out6.read2.target_block.readSeq.blast_rRNA.out -evalue 1 -word_size 5 -outfmt 7 -num_threads 16
perl step5.filter_low_MAPQ.pl out3.read2.circRNA_linked.targets out4.read2.circRNA_linked.bed out6.read2.target_block.readSeq.reMap.bowtie2.sam out6.read2.target_block.readSeq.reMap.hisat2.sam out6.read2.target_block.readSeq.blast_rRNA.out out7.read2.circRNA_linked.targets 2> log.filtered.read2.blocks


# merge
perl step1.integrate.pl > result1.circRNA_linked_target.bed

bedtools intersect -a result1.circRNA_linked_target.bed -b hg19.gene_element.gencode_noncode_intergenic.bed -wa -wb -f 0.51 > result2.target_to_ele.overlap
perl step2.annotate_target_element.pl result1.circRNA_linked_target.bed result2.target_to_ele.overlap > result3.circRNA_linked_target.anno.table
$ less result3.circRNA_linked_target.anno.table | grep "OnTheSenseStrand" > result3.circRNA_linked_target_sense.anno.table

perl step3.creat_network.pl result3.circRNA_linked_target_sense.anno.table

# file for Monte Carlo simulation
perl generate_input_for_simulation.pl result3.circRNA_linked_target_sense.anno.table
perl add_chr_locus_for_target.pl whole_gene_info.txt circRNA_to_target_read_num.txt interactions_for_simulation.txt


