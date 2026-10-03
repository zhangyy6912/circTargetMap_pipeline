#!/usr/bin/perl
die "perl $0 input.read star.align.sam star.chimeric.sam\n" if(@ARGV != 3);
my $input_to_star_fastq=shift;
my $star_align_sam=shift;
my $star_chimeric_sam=shift;

my %is_star_chimeric;
open(SCS,$star_chimeric_sam) || die;
while(my $line=<SCS>){
	chomp $line;
	if($line=~/^@/){
		next;
	}
	my $id=(split/\s+/,$line)[0];
	$id=~s/^Head_//;
	$id=~s/^Tail_//;
	$is_star_chimeric{$id}=1;
}

my %is_fully_aligned;
open(SAS,$star_align_sam) || die;
while(my $line=<SAS>){
	chomp $line;
	if($line=~/^@/){
		next;
	}
	my @sub=split/\s+/,$line;
	if($sub[1] > 255){
		next;
	}
	my $cigar=$sub[5];
	my $mapped_len;
	while($cigar=~/(\d+)M/g){
		$mapped_len+=$1;
	}
	if(length($sub[9])-$mapped_len < 15){	#chimeric fragment > 15
		$is_fully_aligned{$sub[0]}=1;
	}
	#warn $line,"\n";
	#warn $mapped_len,"\taa\n";
	#exit;
}

open(ISF,$input_to_star_fastq) || die;
while(my $id_line=<ISF>){
	my $seq=<ISF>;
	my $symbol=<ISF>;
	my $qual=<ISF>;
	my $short_id=(split/\s+/,$id_line)[0];
	$short_id=~s/^@//;
	if(exists $is_star_chimeric{$short_id} or exists $is_fully_aligned{$short_id}){
		next;
	}
	else{
		print $id_line,$seq,$symbol,$qual;
	}
}
	
