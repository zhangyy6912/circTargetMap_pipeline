#!/usr/bin/perl

use strict;
die "perl $0 gene_info  circRNA2target  circRNA2target_edit\n" if(@ARGV != 3);
my $gene_info=shift;
my $circRNA2target=shift;
my $circRNA2target_edit=shift;

my %gene_info;
open(GN,$gene_info) || die;
while (my $line=<GN>){
    chomp $line;
    #chr1    11868   14412   genecode        DDX11L1 pseudogene      +
    my @sub=split/\s+/,$line;
    my $gene=$sub[4];
    $gene_info{$gene}{'chr'}=$sub[0];
    $gene_info{$gene}{'start'}=$sub[1];
    $gene_info{$gene}{'end'}=$sub[2];
    $gene_info{$gene}{'type'}=$sub[5];
    $gene_info{$gene}{'strand'}=$sub[6];
}

open (OUT, ">$circRNA2target_edit") || die;
open (IN, $circRNA2target) || die;
while(my $line=<IN>){
    chomp $line;
    #chrX	139865341	139866825	circRNA	CDR1as	-  POLE    1       chr12
    my @sub=split/\s+/,$line;
    my $target=$sub[6];
    my $target_chr=$sub[8];
    if (exists $gene_info{$target} and $target_chr eq $gene_info{$target}{'chr'}){
        print OUT "$sub[0]\t$sub[1]\t$sub[2]\t$sub[3]\t$sub[4]\t$sub[5]\t$gene_info{$target}{'chr'}\t$gene_info{$target}{'start'}\t$gene_info{$target}{'end'}\t$gene_info{$target}{'type'}\t$target\t$gene_info{$target}{'strand'}\t$sub[7]\n";
    }
}
close OUT;
