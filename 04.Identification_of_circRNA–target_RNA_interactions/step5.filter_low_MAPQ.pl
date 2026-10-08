#!/usr/bin/perl
die "perl $0 out3.read1.circRNA_linked.targets out4.read1.circRNA_linked.bed out6.target_block.readSeq.reMap.sam out6.target_block.readSeq.reMap.sam out6.read1.target_block.readSeq.blast_rRNA.out outprefix\n" if(@ARGV != 6);
my $circRNA_linked_table=shift;
my $circRNA_linked_bed=shift;
my $circRNA_target_reMap_sam_byBowtie2=shift;
my $circRNA_target_reMap_sam_byHisat2=shift;
my $circRNA_target_blast_out=shift;
my $outprefix=shift;

my %read_rRNA_fraction;
open(CTBO,$circRNA_target_blast_out) || die;
while(my $line=<CTBO>){
        chomp $line;
        if($line=~/^#/){
                next;
        }
        my @sub=split/\s+/,$line;
        $sub[0]=~s/^@//;
        if(exists $read_rRNA_fraction{$sub[0]}){
                next;
        }
        else{
                $read_rRNA_fraction{$sub[0]}=$sub[3];
        }
}

my %read_block_MAPQ_bowtie2;
open(CTRS,$circRNA_target_reMap_sam_byBowtie2) || die;
while(my $line=<CTRS>){
	chomp $line;
	if($line=~/^@/){
		next;
	}
	my @sub=split/\s+/,$line;
	if(exists $read_block_MAPQ_bowtie2{$sub[0]}){
		die;
	}
	
	my $map_strand;
	if($sub[1] eq "0"){
		$map_strand="+";
	}
	elsif($sub[1] eq "16"){
		$map_strand="-";
	}
	else{
		$map_strand="NA";
	}

	$read_block_MAPQ_bowtie2{$sub[0]}=$sub[4]."\t".$sub[2]."\t".$sub[3]."\t".$map_strand;
}

my %read_block_MAPQ_hisat2;
open(CTRSHI,$circRNA_target_reMap_sam_byHisat2) || die;
while(my $line=<CTRSHI>){
        chomp $line;
        if($line=~/^@/){
                next;
        }
        my @sub=split/\s+/,$line;
        if($sub[1] > 255){      #not primary
                next;
        }
        if(exists $read_block_MAPQ_hisat2{$sub[0]}){
                die;
        }
        my $map_strand;
        if($sub[1] eq "0"){
                $map_strand="+";
        }
        elsif($sub[1] eq "16"){
                $map_strand="-";
        }
        else{
                $map_strand="NA";
        }
        $read_block_MAPQ_hisat2{$sub[0]}=$sub[4]."\t".$sub[2]."\t".$sub[3]."\t".$map_strand;
}


open(OUTT,">out7.$outprefix.targets") || die;
open(OUTB,">out7.$outprefix.bed") || die;
open(CT,$circRNA_linked_table) || die;
open(CB,$circRNA_linked_bed) || die;
while(my $line_b=<CB>){
	my $line_t=<CT>;
	chomp $line_b;
	chomp $line_t;
	my @sub_b=split/\s+/,$line_b;
	my @sub_t=split/\s+/,$line_t;
	
	if($sub_b[3] ne $sub_t[2]){	#check whether these two files are matched
		die;
	}

	my @block_info=split/:/,$sub_b[5];
	my $key=$sub_b[4].":".$block_info[1].":".$block_info[2];
	my $block_tag_len=$block_info[2]-$block_info[1];

        my $is_bad_in_bowtie2;
        if(!exists $read_block_MAPQ_bowtie2{$key}){
                $is_bad_in_bowtie2=1;
        }
        else{
                my ($mapq,$map_chr,$map_loci,$map_strand)=split/\s+/,$read_block_MAPQ_bowtie2{$key};
                if($mapq < 20){
                        $is_bad_in_bowtie2=1;
                }
                if($map_chr ne $block_info[3]){
                        $is_bad_in_bowtie2=1;
                }
                if($map_strand ne $block_info[6]){
                        $is_bad_in_bowtie2=1;
                }
                if($map_loci ne $block_info[4]){
                        $is_bad_in_bowtie2=1;
                }
        }

        my $is_bad_in_hisat2;
        if(!exists $read_block_MAPQ_hisat2{$key}){
                $is_bad_in_hisat2=1;
        }
        else{
                my ($mapq,$map_chr,$map_loci,$map_strand)=split/\s+/,$read_block_MAPQ_hisat2{$key};
                if($mapq < 20){
                        $is_bad_in_hisat2=1;
                }
                if($map_chr ne $block_info[3]){
                        $is_bad_in_hisat2=1;
                }
                if($map_strand ne $block_info[6]){
                        $is_bad_in_hisat2=1;
                }
                if($map_loci ne $block_info[4]){
                        $is_bad_in_hisat2=1;
                }
        }

        if($is_bad_in_bowtie2 and $is_bad_in_hisat2){
                warn $line_b,"\t","not uniquely mapped by both bowtie2 and hisat2\n";
                next;
	}

        my $rRNA_fraction=$read_rRNA_fraction{$key}/$block_tag_len;
        if($rRNA_fraction > 0.8){
                warn $line_b,"\t",$rRNA_fraction,"\tfrom rRNA; maybe\n";
                next;
        }

	print OUTT $line_t,"\n";
	print OUTB $line_b,"\n";
}

close OUTT;
close OUTB;

