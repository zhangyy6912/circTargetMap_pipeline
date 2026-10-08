#!/usr/bin/perl
die "perl $0 in.fa\n" if(@ARGV != 1);
my $in_fa=shift;

my $seq_id=0;

$/=">";
open(IN,$in_fa) || die;
<IN>;
while(my $block=<IN>){
	chomp $block;
	my @lines=split/\n/,$block;
	my $first=shift @lines;
	my $id=(split/\s+/,$first)[0];
	my $seq=join"",@lines;
	$seq_id++;
	print "@",$id,".$seq_id\n";
	print $seq,"\n";
	print "+\n";
	print "F" x length($seq);
	print "\n";
}

	
	

