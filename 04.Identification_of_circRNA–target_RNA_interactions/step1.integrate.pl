#!/usr/bin/perl

my @samples=(
"1.HeLa_rep1",
#"2.HeLa_rep2",
#"3.H1_rep1",
#"4.H1_rep2",
#"5.hNPC_rep1",
#"6.hNPC_rep2",
"7.neuron_rep1",
#"8.neuron_rep2",
#"9.GM12878_rep1",
#"10.GM12878_rep2",
#"11.IMR90_rep1",
#"12.IMR90_rep2",
#"13.HT29_rep1",
#"14.HT29_rep2",
#"23.K562_rep1",
#"24.K562_rep2",
#"25.HepG2_rep1",
#"26.HepG2_rep2",
#"27.addhNPC_rep1",
#"28.addhNPC_rep2",
#"29.addK562_rep1",
#"30.addK562_rep2",
#"31.addGM12878_rep1",
#"32.addGM12878_rep2",
#"33.addH1_rep1",
#"34.addH1_rep2",
#"35.addHepG2_rep1",
#"36.addHepG2_rep2",
#"37.addIMR90_rep1",
#"38.addIMR90_rep2",
#"39.293_rep1",
#"40.293_rep2",
#"41.U87WT_rep1",
#"42.U87WT_rep2"
);

my $index_of_target;
foreach my $d (@samples){
        #remove the circle RNA that from repeats and rRNA regions 
	`bedtools intersect -a ../$d/out7.read1.circRNA_linked.bed -b input0.rRNA_SimpleRep_on_genome.hg19.bed -f 0.8 -v > temp.read1.target.bed`;

	my %circRNA_blocks_read1;
	open(BLOCKSA,"../$d/out7.read1.circRNA_linked.targets") || die;
	while(my $line=<BLOCKSA>){
		chomp $line;
		my @sub=split/\s+/,$line;
		$circRNA_blocks_read1{$sub[2]."\t".$sub[0]."\t".$sub[3]}=$sub[1];
	}
	close BLOCKSA;

	open(OUT,">result0.$d.circRNA_linked.bed") || die;
	open(INA,"temp.read1.target.bed") || die;
	while(my $line=<INA>){
		chomp $line;
		my @sub=split/\s+/,$line;
		$index_of_target++;
		print OUT $line,"\t",$d,"#R1\tindex$index_of_target\n";
		print $line,"\t",$d,"#R1\tindex$index_of_target\t";
	
		if(!exists $circRNA_blocks_read1{$sub[3]."\t".$sub[4]."\t".$sub[5]}){
			die;
		}
		else{
			print "circBlock:";
			print $circRNA_blocks_read1{$sub[3]."\t".$sub[4]."\t".$sub[5]};
			print "\n";
		}
	}


	`bedtools intersect -a ../$d/out7.read2.circRNA_linked.bed -b input0.rRNA_SimpleRep_on_genome.hg19.bed -f 0.8 -v > temp.read2.target.bed`;

	my %circRNA_blocks_read2;
	open(BLOCKSB,"../$d/out7.read2.circRNA_linked.targets") || die;
	while(my $line=<BLOCKSB>){
		chomp $line;
		my @sub=split/\s+/,$line;
		$circRNA_blocks_read2{$sub[2]."\t".$sub[0]."\t".$sub[3]}=$sub[1];
	}
	close BLOCKSB;

	open(INB,"temp.read2.target.bed") || die;
	while(my $line=<INB>){
		chomp $line;
		my @sub=split/\s+/,$line;
		$index_of_target++;
		print OUT $line,"\t",$d,"#R2\tindex$index_of_target\n";
		print $line,"\t",$d,"#R2\tindex$index_of_target\t";
		if(!exists $circRNA_blocks_read2{$sub[3]."\t".$sub[4]."\t".$sub[5]}){
			die;
		}
		else{
			print "circBlock:";
			print $circRNA_blocks_read2{$sub[3]."\t".$sub[4]."\t".$sub[5]};
			print "\n";
		}
	}
	close INA;
	close INB;
	close OUT;
	
	`rm -rf temp.read1.target.bed temp.read2.target.bed`;

}
