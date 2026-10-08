#!/usr/bin/perl
die "perl $0 out4.read1.circRNA_linked.bed read.fq\n"if(@ARGV != 2);
my $target_of_circ_table=shift;
my $input_read_fq=shift;

my %should_remap;
open(TCT,$target_of_circ_table) || die;
while(my $line=<TCT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my @block_info=split/:/,$sub[5];
	$should_remap{$sub[4]}{$block_info[1].":".$block_info[2]}=1;
	#print $sub[4],"\t",$block_info[1].":".$block_info[2],"\taa\n";
}

open(IRF,"gzip -c -d $input_read_fq|") || die;
while(my $name=<IRF>){
	my $seq=<IRF>;
	my $symbol=<IRF>;
	my $qual=<IRF>;
	chomp $name;
	chomp $seq;
	chomp $qual;
	$name=~s/^@//;
	if(exists $should_remap{$name}){
		my @blocks=keys %{$should_remap{$name}};
		foreach my $se (@blocks){
			my ($start,$end)=split/:/,$se;
			my $subSeq=substr($seq,$start,$end-$start);
			my $subQual=substr($qual,$start,$end-$start);
			print "@",$name,":$start:$end\n";
			print $subSeq,"\n";
			print "+\n";
			print $subQual,"\n";
		}
	}
}

	
