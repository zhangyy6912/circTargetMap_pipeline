#!/usr/bin/perl
use strict;
die "perl $0 input.blocks.read1 input.blocks.read2 circBedPair\n" if(@ARGV != 3);
my $blocks_read1=shift;
my $blocks_read2=shift;
my $circ_splice_sites_bedpair=shift;

my $slop=10;
my %read1_circRNA=&get_circRNA_junction($blocks_read1);
my %read2_circRNA=&get_circRNA_junction($blocks_read2);
my %read1_target=&get_target($blocks_read1);
my %read2_target=&get_target($blocks_read2);

#get the circRNA and target pair
foreach my $key (keys %read1_circRNA){ #get the target of read1
    if(exists $read2_target{$key}){
            print $key, "\t", "$read1_circRNA{$key}{'block'}", "\t", $read1_circRNA{$key}{'circRNA'}, "\t", $read2_target{$key}, "\t", "1-2\n";
    }
}

foreach my $key (keys %read2_circRNA){ #get the target of read1
    #print $key, "\n";
    if(exists $read1_target{$key}){
            print $key, "\t", $read2_circRNA{$key}{'block'}, "\t", $read2_circRNA{$key}{'circRNA'}, "\t", $read1_target{$key}, "\t", "2-1\n";
    }
}

sub get_circRNA_junction() {
    my ($blocks)=@_;
    #print $blocks, "\n";
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
    		if($strand_this eq "+" and $end_next < $start_this){#chimeric, the coordination here is chr localization not read
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

    #get the chimeric reads that overlap with circRNAs  
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
    		$is_circRNA_junction{$readID}{'block'}=$this_block_id.":".$next_block_id;
                $is_circRNA_junction{$readID}{'circRNA'}=$sub[16].":".$sub[18];
    	}
    	else{
    		die;
    	}
    }
    close LOP;
    return %is_circRNA_junction;
    `rm chimeric_pairBlocks.bedpair overlap.chimeric_to_circRNA.list`;
}

sub get_target{
    my %is_target={};
    my ($blocks)=@_;
    #print $blocks, "\n";
    open(BKT,$blocks) || die;
    while(my $line=<BKT>){
    	chomp $line;
        my @sub=split/\s+/,$line;
        my $readID=shift @sub;
        if($#sub+1 ==1){#only one block
            $is_target{$readID}=$sub[0];
        }
        elsif($#sub+1 ==2){
    		my @info_this=split/:/,$sub[0];
    		my @info_next=split/:/,$sub[1];
    	
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
    		
    		if(($chr_this eq $chr_next) and ($strand_this ne $strand_next) and ($RNAstrand_this ne $RNAstrand_next)){
    		if($strand_this eq "+" and $start_next > $end_this){#normal junction
                        $is_target{$readID}=$sub[0];
    		}
    		elsif($strand_this eq "-" and $end_next < $start_this){
                        $is_target{$readID}=$sub[0];
    		}
    		else{
    			next;
    		}
             }
    	}
    }
    close BK;
    return %is_target;
}


