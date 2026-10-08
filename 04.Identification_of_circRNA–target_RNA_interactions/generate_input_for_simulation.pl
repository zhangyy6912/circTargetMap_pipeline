#!/usr/bin/perl
use strict;
die "perl $0 result3.circRNA_linked_target.anno.table\n" if(@ARGV != 1);
my $circRNA_target_table=shift;

my %circ_to_target;

my %target_chr;
my %circ;

open(CTT,$circRNA_target_table) || die;
while(my $line=<CTT>){
	chomp $line;
	my @sub=split/\s+/,$line;
        #get the circRNA information
	my $circ=$sub[3];
        if($circ=~/(\S+):\d+:(chr.+)_(\d+)_(\d+):\d+:(\S+)/){
    	    $circ{$circ}{'name'}=$1;
    	    $circ{$circ}{'chr'}=$2;
    	    $circ{$circ}{'start'}=$3;
    	    $circ{$circ}{'end'}=$4;
            $circ{$circ}{'strand'}=$5;
            $circ{$circ}{'type'}='circRNA';
        }
        else{
                die;
        }

        #get the target information
	my $target_info=$sub[5];
        my $target_chr;
        my $target_start;
        my $target_end;
        my $target_strand;
        if($target_info=~/\d+:\d+:\d+:(chr.+):(\d+):(\d+):\S+:(\S+)/){
    	    $target_chr=$1;
    	    $target_start=$2;
    	    $target_end=$3;
            if($4 eq "Plus"){
                $target_strand="+";
            }
            elsif($4 eq "Minus"){
                $target_strand="-";
            }
            else{
                die;
            }
        }
        else{
                die;
        }
        #if the circRNA binding site located in the circRNA region, remove this read
        if(($target_chr == $circ{$circ}{'chr'}) and ($target_strand == $circ{$circ}{'strand'}) and (int($target_start) >= int($circ{$circ}{'start'})) and (int($target_end) <= int($circ{$circ}{'end'}))){
            next;
        }
	if($sub[9] eq "Intergenic" or $sub[9] eq "Unknown"){
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
		my $target=$symbol;
		if(exists $non_redundant_target{$symbol}){
			next;
		}
		$non_redundant_target{$symbol}=1;
                if($sub[-1] eq "OnTheSenseStrand"){
		    $circ_to_target{$circ}{$target}{$sub[4]."\t".$sub[6]}=1;
		    $target_chr{$target}=$sub[0];
                } 
	}
        #hsa-TUBA1A_0003:65660:chr12_49523025_49579773:0:-
        #print $circ, "\n";
}
close CTT;

open(NET,">circRNA_to_target_read_num.txt") || die;

print NET "circRNA_chr\tcircRNA_start\tcircRNA_end\tcircRNA_type\tcircRNA\tcircRNA_strand\tTargetRNA\tStrength\ttargetChr\n";
foreach my $circ (keys %circ_to_target){
        #print $circ, "\n";
	foreach my $target (keys %{$circ_to_target{$circ}}){
		print NET $circ{$circ}{'chr'}, "\t", $circ{$circ}{'start'},"\t", $circ{$circ}{'end'},"\t", $circ{$circ}{'type'},"\t", $circ{$circ}{'name'},"\t", $circ{$circ}{'strand'},"\t";
                print NET $target,"\t";
		my @all_read=keys %{$circ_to_target{$circ}{$target}};
		print NET $#all_read+1,"\t";
                print NET $target_chr{$target},"\n";
	};
}

close NET;
