#!/usr/bin/perl
die "perl $0 in.sam_by_bwa 1/2 \n" if(@ARGV != 2);
my $input_sam=shift;
my $read1_or_read2=shift;

my %reads_input_to_bwa;
my %reads_blocks;
my %read_len;
my %reads_seq;
my %reads_qual;
open(IS,$input_sam) || die;
while(my $line=<IS>){
	chomp $line;
	if($line=~/^@/){
		next;
	}
	my @sam_info=split/\s+/,$line;
	my $read_id=$sam_info[0];
	#if($read_id ne "622867"){
	#if($read_id ne "622880"){
	#	next;	
	#}
	$reads_input_to_bwa{$read_id}=1;
	if(exists $reads_seq{$read_id}){
	}
	else{
		if($sam_info[1] eq "0" or $sam_info[1] eq "2048"){
			$read_len{$read_id}=length($sam_info[9]);
			$reads_seq{$read_id}=$sam_info[9];
			$reads_qual{$read_id}=$sam_info[10];
		}
		elsif($sam_info[1] eq "16" or $sam_info[1] eq "2064"){
			my $tmp_seq=$sam_info[9];
			$tmp_seq=reverse $tmp_seq;
			$tmp_seq=~tr/ATCG/TAGC/;
			my $tmp_qual=$sam_info[10];
			$tmp_qual=reverse $tmp_qual;
			$read_len{$read_id}=length($sam_info[9]);
			$reads_seq{$read_id}=$tmp_seq;
			$reads_qual{$read_id}=$tmp_qual;
		}
		elsif($sam_info[1] eq "4"){
			next;
		}
	}
	my ($start_in_read_this,$end_in_read_this,$chr_name_this,$start_in_chr_this,$end_in_chr_this,$strand_this)=obtain_read_blocks($line);
	if(exists $reads_blocks{$read_id}){
		my @tmp=@{$reads_blocks{$read_id}};
		if($#tmp > 3){	#only store four fragment; more fragment is more complex
			next;
		}
		push (@tmp,[$start_in_read_this,$end_in_read_this,$chr_name_this,$start_in_chr_this,$end_in_chr_this,$strand_this]);
		@{$reads_blocks{$read_id}}=@tmp;
	}
	else{
		@{$reads_blocks{$read_id}}=[$start_in_read_this,$end_in_read_this,$chr_name_this,$start_in_chr_this,$end_in_chr_this,$strand_this];
	}
}


my $pet_index;
foreach my $read_id (keys %reads_blocks){
	my @raw_blocks=@{$reads_blocks{$read_id}};
	@raw_blocks=sort {$a->[0] <=> $b->[0]} @raw_blocks;
	my @blocks;

	foreach my $i (0..$#raw_blocks){
		my $conflict_i;
		foreach my $j (0..$#raw_blocks){
			if($j==$i){
				next;
			}
			$conflict_i=$conflict_i+count_conflict_bases($raw_blocks[$i][0],$raw_blocks[$i][1],$raw_blocks[$j][0],$raw_blocks[$j][1]);
		}
		#print $raw_blocks[$i][0],"\t",$raw_blocks[$i][1],"\t",$conflict_i,"\n";
		if($conflict_i > 0.25*($raw_blocks[$i][1]-$raw_blocks[$i][0])){
			next;
		}
		else{
			push (@blocks,[$raw_blocks[$i][0],$raw_blocks[$i][1],$raw_blocks[$i][2],$raw_blocks[$i][3],$raw_blocks[$i][4],$raw_blocks[$i][5]]);
			#$effective_len=$effective_len+$raw_blocks[$i][1]-$raw_blocks[$i][0];
		}
	}

	my %effective_mapped;
	foreach my $i (0..$#blocks){
		foreach ($blocks[$i][0]+1..$blocks[$i][1]){
			$effective_mapped{$_}++;
		}
	}

	my $effective_len;
	foreach (keys %effective_mapped){
		if($effective_mapped{$_} == 1){
			$effective_len++;
		}
	}

	if($effective_len < 0.9*$read_len{$read_id}){
		print UNMAP $read_id,"\n";
		next;
	}

	if($#blocks < 2){	#specific for circRNAJunc-to-RNA interactions
		next;
	}

	print $read_id,"\t";
	foreach my $i (0..$#blocks){
		my $this_start=$blocks[$i][0];
		my $this_end=$blocks[$i][1];
		my $this_chr=$blocks[$i][2];
		my $this_chr_start=$blocks[$i][3];
		my $this_chr_end=$blocks[$i][4];
		my $this_strand=$blocks[$i][5];
		my $plus_or_minus_this;
		if(($this_strand eq "+" and $read1_or_read2 eq "2") or ($this_strand eq "-" and $read1_or_read2 eq "1")){
			$plus_or_minus_this="Plus";
		}
		elsif(($this_strand eq "-" and $read1_or_read2 eq "2") or ($this_strand eq "+" and $read1_or_read2 eq "1")){
			$plus_or_minus_this="Minus";
		}

		print $i,":",$this_start,":",$this_end,":",$this_chr,":",$this_chr_start,":",$this_chr_end,":",$this_strand,":",$plus_or_minus_this,"\t";
	}
	print "\n";
}


sub count_conflict_bases{
	my @four=@_;
	my $min=$four[0] < $four[2] ? $four[2] : $four[0];
	my $max=$four[1] < $four[3] ? $four[1] : $four[3];
	if($max > $min){
		return $max-$min;
	}
	return 0;
}
		
		


sub obtain_read_blocks{
	my $sam_line=shift;
	my @sub=split/\s+/,$sam_line;
	my $strand;
	if($sub[1] eq "0" or $sub[1] eq "2048"){
		$strand="+";
	}
	elsif($sub[1] eq "16" or $sub[1] eq "2064"){
		$strand="-";
	}
	else{
		die;
	}
	my $cigar=$sub[5];
	$cigar=~s/H/S/g;
	my $start_in_read;
	my $end_in_read;
	my $chr_name=$sub[2];
	my $start_in_chr=$sub[3];
	my $end_in_chr;
	if($strand eq "+"){
		if($cigar=~/^(\d+)S(\d+)M$/){
			$start_in_read=$1;
			$end_in_read=$1+$2;
			$end_in_chr=$start_in_chr+$2;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)I(\d+)M$/){
			$start_in_read=$1;
			$end_in_read=$1+$2+$3+$4;
			$end_in_chr=$start_in_chr+$2+$4;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)D(\d+)M$/){
			$start_in_read=$1;
			$end_in_read=$1+$2+$4;
			$end_in_chr=$start_in_chr+$2+$3+$4;
		}
		elsif($cigar=~/^(\d+)M(\d+)S$/){
			$start_in_read=0;
			$end_in_read=$1;
			$end_in_chr=$start_in_chr+$1;
		}
		elsif($cigar=~/^(\d+)M(\d+)I(\d+)M(\d+)S$/){
			$start_in_read=0;
			$end_in_read=$1+$2+$3;
			$end_in_chr=$start_in_chr+$1+$3;			
		}
		elsif($cigar=~/^(\d+)M(\d+)D(\d+)M(\d+)S$/){
			$start_in_read=0;
			$end_in_read=$1+$3;
			$end_in_chr=$start_in_chr+$1+$2+$3;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)S$/){
			$start_in_read=$1;
			$end_in_read=$1+$2;
			$end_in_chr=$start_in_chr+$2;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)D(\d+)M(\d+)S$/){
			$start_in_read=$1;
			$end_in_read=$1+$2+$4;
			$end_in_chr=$start_in_chr+$2+$3+$4;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)I(\d+)M(\d+)S$/){
			$start_in_read=$1;
			$end_in_read=$1+$2+$3+$4;
			$end_in_chr=$start_in_chr+$2+$4;
		}
		else{
			#warn $cigar,"\t";
			next;
		}
	}
	elsif($strand eq "-"){
		if($cigar=~/^(\d+)S(\d+)M$/){
			$start_in_read=0;
			$end_in_read=$2;
			$end_in_chr=$start_in_chr+$2;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)D(\d+)M$/){
			$start_in_read=0;
			$end_in_read=$2+$4;
			$end_in_chr=$start_in_chr+$2+$3+$4;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)I(\d+)M$/){
			$start_in_read=0;
			$end_in_read=$2+$3+$4;
			$end_in_chr=$start_in_chr+$2+$4;
		}
		elsif($cigar=~/^(\d+)M(\d+)S$/){
			$start_in_read=$2;
			$end_in_read=$2+$1;
			$end_in_chr=$start_in_chr+$1;
		}
		elsif($cigar=~/^(\d+)M(\d+)D(\d+)M(\d+)S$/){
			$start_in_read=$4;
			$end_in_read=$4+$3+$1;
			$end_in_chr=$start_in_chr+$1+$2+$3;
		}
		elsif($cigar=~/^(\d+)M(\d+)I(\d+)M(\d+)S$/){
			$start_in_read=$4;
			$end_in_read=$4+$3+$2+$1;
			$end_in_chr=$start_in_chr+$1+$3;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)S$/){
			$start_in_read=$3;
			$end_in_read=$3+$2;
			$end_in_chr=$start_in_chr+$2;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)D(\d+)M(\d+)S$/){
			$start_in_read=$5;
			$end_in_read=$5+$4+$2;
			$end_in_chr=$start_in_chr+$2+$3+$4;
		}
		elsif($cigar=~/^(\d+)S(\d+)M(\d+)I(\d+)M(\d+)S$/){
			$start_in_read=$5;
			$end_in_read=$5+$4+$3+$2;
			$end_in_chr=$start_in_chr+$2+$4;
		}
		else{
			#warn $sam_line,"\n";
			#warn $cigar,"\t";
			next;
		}
	}
	else{
		die;
	}

	return  ($start_in_read,$end_in_read,$chr_name,$start_in_chr,$end_in_chr,$strand);	
}









