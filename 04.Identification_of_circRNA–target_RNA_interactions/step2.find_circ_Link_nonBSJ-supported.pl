#!/usr/bin/perl
use strict;
die "perl $0 input.blocks circBedPair\n" if(@ARGV != 2);
my $blocks=shift;
my $circ_splice_sites_bedpair=shift;

my $slop=10;
open(TMP,">chimeric_pairBlocks.bedpair") || die;
open(BK,$blocks) || die;
while(my $line=<BK>){
	chomp $line;;
	my @sub=split/\s+/,$line;
	my $readID=shift @sub;
	foreach (0..$#sub-1){
		my @info_this=split/:/,$sub[$_];
		my @info_next=split/:/,$sub[$_+1];
	
		my $order_this=$info_this[0];	
		my $chr_this=$info_this[3];
		my $start_this=$info_this[4];
		my $end_this=$info_this[5];
		my $strand_this=$info_this[6];
		my $RNAstrand_this=$info_this[7];
		
		my $order_next=$info_next[0];
		my $chr_next=$info_next[3];
		my $start_next=$info_next[4];
		my $end_next=$info_next[5];
		my $strand_next=$info_next[6];
		my $RNAstrand_next=$info_next[7];

		if($end_this < $start_this or $end_next < $start_next){	#for checking
			die;
		}
		
		if($chr_this ne $chr_next){
			next;
		}
		if($strand_this ne $strand_next){
			next;
		}
		if($RNAstrand_this ne $RNAstrand_next){	#same mapping strand must from same RNA strand
			die;
		}
		if($strand_this eq "+" and $end_next < $start_this){#chimeric
			print TMP $chr_this,"\t",$end_this-$slop,"\t",$end_this+$slop,"\t";
			print TMP $chr_next,"\t",$start_next-$slop,"\t",$start_next+$slop,"\t";
			print TMP $readID,"($order_this:$order_next)\t255\t";
			if($RNAstrand_this eq "Plus"){
				print TMP "+\t+\n";
			}
			elsif($RNAstrand_this eq "Minus"){
				print TMP "-\t-\n";
			}
			else{
				die;
			}
		}
		elsif($strand_this eq "-" and $start_next > $end_this){
			print TMP $chr_next,"\t",$end_next-$slop,"\t",$end_next+$slop,"\t";
			print TMP $chr_this,"\t",$start_this-$slop,"\t",$start_this+$slop,"\t";
			print TMP $readID,"($order_this:$order_next)\t255\t";
			if($RNAstrand_this eq "Plus"){
				print TMP "+\t+\n";
			}
			elsif($RNAstrand_this eq "Minus"){
				print TMP "-\t-\n";
			}
			else{
				die;
			}
		}
		else{
			next;
		}
	}
}
close BK;
close TMP;

`bedtools pairtopair -a chimeric_pairBlocks.bedpair -b $circ_splice_sites_bedpair > overlap.chimeric_to_circRNA.list`;	#we require strand

my %is_circRNA_junction;
open(OLP,"overlap.chimeric_to_circRNA.list") || die;
while(my $line=<OLP>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($sub[6]=~/(.+)\((\d+):(\d+)\)$/){
		my $readID=$1;
		my $this_block_id=$2;
		my $next_block_id=$3;
		my $circ_abundance=(split/:/,$sub[16])[-1];	#not use
		$is_circRNA_junction{$readID}{$this_block_id."\t".$sub[16].":".$sub[18]}=1;
	}
	else{
		die;
	}
}
close LOP;

open(BKT,$blocks) || die;
while(my $line=<BKT>){
	chomp $line;
        my @sub=split/\s+/,$line;
        my $readID=shift @sub;

	
	if(exists $is_circRNA_junction{$readID}){
		my %candidate_target_block;
		foreach (0..$#sub){
			my $index=(split/:/,$sub[$_])[0];
			if(exists $candidate_target_block{$index}){
				die;
			}
			$candidate_target_block{$index}=$sub[$_];
		}

		foreach my $index_circ (sort {$a<=>$b} keys %{$is_circRNA_junction{$readID}}){
			my ($p,$this_circ)=split/\t/,$index_circ;

			my $previous_block_index=$p-1;
			if(exists $candidate_target_block{$previous_block_index}){
				print $readID,"\t",$p,":",$p+1,"\t",$this_circ,"\t",$candidate_target_block{$previous_block_index},"\n";
			}
			my $next_block_index=$p+2;
			if(exists $candidate_target_block{$next_block_index}){
				print $readID,"\t",$p,":",$p+1,"\t",$this_circ,"\t",$candidate_target_block{$next_block_index},"\n";
			}
		}
	}
}

sleep(5);
`rm -rf chimeric_pairBlocks.bedpair overlap.chimeric_to_circRNA.list`;

