#!/usr/bin/perl
die "perl $0 out2.read*.circRNA_linked.blocks\n" if(@ARGV != 1);
my $circRNA_linked_targetBlocks=shift;

open(CLT,$circRNA_linked_targetBlocks) || die;
while(my $line=<CLT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my @circ_info=split/:/,$sub[2];
	my @target_info=split/:/,$sub[3];
        my $type=$sub[4];

	my ($circ_chr,$circ_start,$circ_end)=split/_/,$circ_info[2];
	my $circ_strand=$circ_info[4];	#+ or -

	my $target_chr=$target_info[3];
	my $target_start=$target_info[4];
	my $target_end=$target_info[5];
	my $target_strand_symbol=$target_info[7];
	my $target_strand;
	if($target_strand_symbol eq "Minus"){
		$target_strand="-";
	}
	elsif($target_strand_symbol eq "Plus"){
		$target_strand="+";
	}
	else{
		die;
	}

	if($circ_strand ne $target_strand){
		print $line,"\n";
		warn $target_chr,"\t",$target_start,"\t",$target_end,"\t",$sub[2],"\t",$sub[0],"\t",$sub[3],"\n";
		next;
	}
	if($target_chr ne $circ_chr){
		print $line,"\n";
		warn $target_chr,"\t",$target_start,"\t",$target_end,"\t",$sub[2],"\t",$sub[0],"\t",$sub[3],"\n";
		next;
	}	
	if($target_start > $circ_end){
		print $line,"\n";
		warn $target_chr,"\t",$target_start,"\t",$target_end,"\t",$sub[2],"\t",$sub[0],"\t",$sub[3],"\n";
		next;
	}
	if($target_end < $circ_start){
		print $line,"\n";
		warn $target_chr,"\t",$target_start,"\t",$target_end,"\t",$sub[2],"\t",$sub[0],"\t",$sub[3],"\n";
		next;
	}
}
