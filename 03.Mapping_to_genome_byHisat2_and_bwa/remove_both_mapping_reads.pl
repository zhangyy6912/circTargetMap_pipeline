#!/usr/bin/perl
die "perl $0 read1.fq read2.fq read1_hisat2.sam read2_hisat2.sam read1_filter.fq read2_filter.fq \n" if(@ARGV != 6);

my $input_read1=shift;
my $input_read2=shift;
my $input_sam1=shift;
my $input_sam2=shift;
my $output_read1=shift;
my $output_read2=shift;

#load in the read1 sam file
if ($input_sam1 =~ /.gz$/) {
    open(IS1, "gunzip -c $input_sam1 |") || die "can't open pipe to $input_sam1";
}
else {
    open(IS1, $input_sam1) || die "can't open $input_sam1";
}

%read1_map={};
while (my $line = <IS1>){
    chomp $line;
    if($line=~/^@/){
        next;
    }
    my @sam_info=split/\s+/,$line;
    my $read_id=$sam_info[0];
    if($sam_info[1] eq "4"){
        next;
    }
    $read1_map{$read_id}=1;
} 
close IS1;

#load in the read2 sam file
if ($input_sam2 =~ /.gz$/) {
    open(IS2, "gunzip -c $input_sam2 |") || die "can't open pipe to $input_sam2";
}
else {
    open(IS2, $input_sam2) || die "can't open $input_sam2";
}
%read2_map={};
while (my $line = <IS2>){
    chomp $line;
    if($line=~/^@/){
        next;
    }
    my @sam_info=split/\s+/,$line;
    my $read_id=$sam_info[0];
    if($sam_info[1] eq "4"){
        next;
    }
    $read2_map{$read_id}=1;
} 
close IS2;

#load in the read1 and remove the both mapped reads
if ($input_read1 =~ /.gz$/) {
    open(IR1, "gunzip -c $input_read1 |") || die "can't open pipe to $input_read1";
}
else {
    open(IR1, $input_read1) || die "can't open $input_read1";
}
#@ST-E00205:676:HTJL3CCXY:4:2207:3194:10468
#AGATTGTGCCACTGCACTCCAGCCTGGGCGACAGTGCGAGACTCTGTCTC
#+
#AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
open (OR1, ">$output_read1") || die;
while(my $line=<IR1>){
    chomp $line;
    my $id=$line;
    $id=~s/@//g;
    #print $line, "\n";
    my $seq=<IR1>;
    my $sym=<IR1>;
    my $qual=<IR1>;
    if(exists $read1_map{$id} and exists $read2_map{$id}){
        next; #discard the reads the both mapped 
    }
    else{
        #write down read1
        print OR1 "\@$id\n";
        print OR1 "$seq";
        print OR1 "$sym";
        print OR1 "$qual";
    }
}
close IR1;
close OR1;

#load in the read2 and output the read
open (OR2, ">$output_read2") || die;
if ($input_read2 =~ /.gz$/) {
    open(IR2, "gunzip -c $input_read2 |") || die "can't open pipe to $input_read2";
    print "ok\n";
}
else {
    open(IR2, $input_read2) || die "can't open $input_read2";
}
while (my $line =<IR2>){
    chomp $line;
    my $id=$line;
    $id=~s/@//g;
    $seq=<IR2>;
    $sym=<IR2>;
    $qual=<IR2>;
    if(exists $read1_map{$id} and exists $read2_map{$id}){
        next; #discard the reads the both mapped   
    }
    else{
        #write down read2
        print OR2 "\@$id\n";
        print OR2 "$seq";
        print OR2 "$sym";
        print OR2 "$qual";
    }
}
close IR2;
close OR2;

