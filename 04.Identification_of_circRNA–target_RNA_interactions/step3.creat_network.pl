#!/usr/bin/perl
use strict;
die "perl $0 result3.circRNA_linked_target.anno.table\n" if(@ARGV != 1);
my $circRNA_target_table=shift;

my %circ_to_target;
my %circ_node;
my %target_node;

my %circ_chr;
my %target_chr;

open(CTT,$circRNA_target_table) || die;
while(my $line=<CTT>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $circ=$sub[3];
	if($sub[9] eq "Intergenic"){
		next;
	}
	my @candidate_target=split/\|/,$sub[9];
	my %non_redundant_target;
	foreach my $sym_ele (@candidate_target){
		my @info=split/#/,$sym_ele;
		if(@info != 2){
			warn $sym_ele,"\taa\n";
			die;
		}
		my $symbol=$info[0];
		my $target=$sub[0].":".$symbol.":".$sub[-1];
		if(exists $non_redundant_target{$symbol}){
			next;
		}
		$non_redundant_target{$symbol}=1;

		$circ_to_target{$circ}{$target}{$sub[4]."\t".$sub[6]}=1;
		$target_node{$target}{$sub[4]."\t".$sub[6]}=1;
		$target_chr{$target}=$sub[0];
	}
	$circ_node{$circ}{$sub[4]."\t".$sub[6]}=1;
	if($circ=~/:(chr.+)_\d+_\d+:/){
		$circ_chr{$circ}=$1;
	}	
	else{
		die;
	}
}

open(NET,">result4.circRNA_to_target.network.txt") || die;
open(NODE,">result5.circRNA_and_target.node.txt") || die;

print NET "circRNA\tTargetRNA\tStrength\tIntraChrorInterChr\n";
foreach my $circ (keys %circ_to_target){
	foreach my $target (keys %{$circ_to_target{$circ}}){
		print NET $circ,"\t",$target,"\t";
		my @all_read=keys %{$circ_to_target{$circ}{$target}};
		print NET $#all_read+1,"\t";
		if($circ_chr{$circ} eq $target_chr{$target}){
			print NET "IntraChr\n";
		}
		else{
			print NET "InterChr\n";
		}
	}
}

my $index;
print NODE "symbol\tcircORtarget\tchr\toutReads\tinReads\treName\n";
foreach my $circ (keys %circ_node){
	print NODE $circ,"\tcircRNA\t";
	print NODE $circ_chr{$circ},"\t";
	my @all_read=keys %{$circ_node{$circ}};
	print NODE $#all_read+1,"\t";
	print NODE "0\t";
	$index++;
	print NODE "circ$index\n";
}

foreach my $target (keys %target_node){
	print NODE $target,"\ttargetRNA\t";
	print NODE $target_chr{$target},"\t";
	print NODE "0\t";
	my @all_read=keys %{$target_node{$target}};
	print NODE $#all_read+1,"\t";
	$index++;
	print NODE "target$index\n";
}







