#!/usr/bin/perl
die "perl $0 input.fastq unmapped.id\n" if(@ARGV != 2);
my $input_fastq=shift;
my $unmapped_id=shift;

my %is_unmapped;
open(UI,$unmapped_id) || die;
while(my $line=<UI>){
	chomp $line;
	$is_unmapped{$line}=1;
}

open(INF,$input_fastq) || die;
while(my $id=<INF>){
	my $seq=<INF>;
	my $symbol=<INF>;
	my $qual=<INF>;
	my $simple_id=(split/\s+/,$id)[0];
	$simple_id=~s/^@//;
	if(exists $is_unmapped{$simple_id}){
		print $id,$seq,$symbol,$qual;
	}
}
