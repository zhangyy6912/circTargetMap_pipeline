#!/usr/bin/perl
die "perl $0 in.sam_by_bwa 1/2 outdir\n" if(@ARGV != 3);
my $input_sam=shift;
my $read1_or_read2=shift;
my $outdir=shift;
$outdir=~s/\/$//;

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

open(UNMAP,">$outdir/unmapped.$read1_or_read2.readID.list") || die;
foreach (keys %reads_input_to_bwa){
	if(exists $reads_blocks{$_}){
	}
	else{
		print UNMAP $_,"\n";
	}
}

open(FIRST,">$outdir/read$read1_or_read2.firstBlcok.list") || die;
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

	#print the first block to mate with its mates 
	my $first_start=$blocks[0][0];
	my $first_end=$blocks[0][1];
	my $first_chr=$blocks[0][2];
	my $first_chr_start=$blocks[0][3];
	my $first_chr_end=$blocks[0][4];
	my $first_strand=$blocks[0][5];
	
	if($first_strand eq "+"){
		print FIRST $read_id,"\t0\t";
		my $offset=$first_chr_end-$first_chr_start-($first_end-$first_start);
		if($offset==0){
			print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start,"M\t*\t0\t0\t";
		}
		elsif($offset > 0){#deletion
			print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start,"M",$offset,"D\t*\t0\t0\t";
		}
		elsif($offset < 0){#insertion
			print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
		}
		my $first_seq=substr($reads_seq{$read_id},$first_start,$first_end-$first_start);
		my $first_qual=substr($reads_qual{$read_id},$first_start,$first_end-$first_start);
		print FIRST $first_seq,"\t",$first_qual,"\n";
	}
	elsif($first_strand eq "-"){
		print FIRST $read_id,"\t16\t";
		my $offset=$first_chr_end-$first_chr_start-($first_end-$first_start);
		if($offset==0){
			print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start,"M\t*\t0\t0\t";
		}
		elsif($offset > 0){#deletion
			print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start,"M",$offset,"D\t*\t0\t0\t";
		}
		elsif($offset < 0){
			print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
		}
		#print FIRST $first_chr,"\t",$first_chr_start,"\t255\t",$first_end-$first_start,"M\t*\t0\t0\t";
		my $first_seq=substr($reads_seq{$read_id},$first_start,$first_end-$first_start);
		my $first_qual=substr($reads_qual{$read_id},$first_start,$first_end-$first_start);
		$first_seq=~tr/ATCG/TAGC/;
		$first_seq=reverse $first_seq;
		$first_qual=reverse $first_qual;
		print FIRST $first_seq,"\t",$first_qual,"\n";
	}
	else{
		die;
	}
	#print over for the first block

	if($#blocks < 1){
		next;
	}

	foreach my $i (0..$#blocks-1){
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

		my $next_start=$blocks[$i+1][0];
		my $next_end=$blocks[$i+1][1];
		my $next_chr=$blocks[$i+1][2];
		my $next_chr_start=$blocks[$i+1][3];
		my $next_chr_end=$blocks[$i+1][4];
		my $next_strand=$blocks[$i+1][5];
		my $plus_or_minus_next;
                if(($next_strand eq "+" and $read1_or_read2 eq "2") or ($next_strand eq "-" and $read1_or_read2 eq "1")){
                        $plus_or_minus_next="Plus";
                }
                elsif(($next_strand eq "-" and $read1_or_read2 eq "2") or ($next_strand eq "+" and $read1_or_read2 eq "1")){
                        $plus_or_minus_next="Minus";
                }

		my $space=$next_start-$this_end;

		my $impossible_situation;
		if($this_end-$this_start < $next_end-$next_start){
			my $please_revise_next;
			my $please_revise_this;
			if($this_strand eq "+"){
				$this_end=$this_end+$space;
				$this_chr_end=$this_chr_end+$space;
			}
			elsif($this_strand eq "-"){
				if($this_chr_start-$space >= 1){
					$this_end=$this_end+$space;
					$this_chr_start=$this_chr_start-$space;
				}
				else{
					$please_revise_next=1;
				}
			}
			else{
				die;
			}
			if($please_revise_next){
	                        if($next_strand eq "+"){
	                                if($next_chr_start-$space >=1){
	                                        $next_start=$next_start-$space;
	                                        $next_chr_start=$next_chr_start-$space;
	                                }
	                                else{
	                                        $please_revise_this=1;
	                                }
	                        }
	                        elsif($next_strand eq "-"){
	                                $next_start=$next_start-$space;
	                                $next_chr_end=$next_chr_end+$space;
	                        }
	                        else{
	                                die;
	                        }
			}
			if($please_revise_this){
				$impossible_situation=1;
			}				
		}
		else{
			my $please_revise_next;
			my $please_revise_this;
			if($next_strand eq "+"){
				if($next_chr_start-$space >=1){
					$next_start=$next_start-$space;
					$next_chr_start=$next_chr_start-$space;
				}
				else{
					$please_revise_this=1;
				}
			}
			elsif($next_strand eq "-"){
				$next_start=$next_start-$space;
				$next_chr_end=$next_chr_end+$space;
			}
			else{
				die;
			}
			if($please_revise_this){
	                        if($this_strand eq "+"){
	                                $this_end=$this_end+$space;
	                                $this_chr_end=$this_chr_end+$space;
	                        }
	                        elsif($this_strand eq "-"){
	                                if($this_chr_start-$space >= 1){
	                                        $this_end=$this_end+$space;
	                                        $this_chr_start=$this_chr_start-$space;
	                                }
	                                else{
	                                        $please_revise_next=1;
	                                }
	                        }
	                        else{
	                                die;
	                        }
			}
                        if($please_revise_next){
				$impossible_situation=1;
			}
		}

		if($impossible_situation){
			warn $read_id,"\thave impossible mapping fragment\n";
			next;
		}

		my $is_singleton;	
		if($this_strand eq $next_strand and $this_chr eq $next_chr){
			if($this_strand eq "+" and $next_chr_start > $this_chr_end){	#singleton
				$is_singleton=1;
				#next;
			}
			elsif($this_strand eq "-" and $next_chr_end < $this_chr_start){	#singleton
				$is_singleton=1;
				#next;
			}
		}

		if($is_singleton){
                        #print $read_id,"\tgapped\n";
                        #print $this_chr,"\t",$this_chr_start,"\t",$this_chr_end,"\t",$this_strand,"\n";
                        #print $next_chr,"\t",$next_chr_start,"\t",$next_chr_end,"\t",$next_strand,"\n";

			$pet_index++;
			if($this_strand eq "+" and $next_strand eq "+"){
				print "Part_",$pet_index,"_",$plus_or_minus_this,"_",$read_id,"\t0\t";
				my $offset=$this_chr_end-$this_chr_start-($this_end-$this_start);
				if($offset==0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start,"M\t*\t0\t0\t";
				}
				elsif($offset > 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start,"M",$offset,"D\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
				}
				my $this_seq=substr($reads_seq{$read_id},$this_start,$this_end-$this_start);
				my $this_qual=substr($reads_qual{$read_id},$this_start,$this_end-$this_start);
				print $this_seq,"\t",$this_qual,"\n";

				print "Part_",$pet_index,"_",$plus_or_minus_next,"_",$read_id,"\t0\t";
				my $offset=$next_chr_end-$next_chr_start-($next_end-$next_start);
				if($offset==0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start,"M\t*\t0\t0\t";
				}
				elsif($offset > 0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start,"M",$offset,"D\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
				}
				my $next_seq=substr($reads_seq{$read_id},$next_start,$next_end-$next_start);
				my $next_qual=substr($reads_qual{$read_id},$next_start,$next_end-$next_start);
				print $next_seq,"\t",$next_qual,"\n";
			}
			elsif($this_strand eq "-" and $next_strand eq "-"){
				print "Part_",$pet_index,"_",$plus_or_minus_next,"_",$read_id,"\t16\t";
				my $offset=$next_chr_end-$next_chr_start-($next_end-$next_start);
				if($offset==0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start,"M\t*\t0\t0\t";
				}
				elsif($offset > 0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start,"M",$offset,"D\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
				}
				my $next_seq=substr($reads_seq{$read_id},$next_start,$next_end-$next_start);
				my $next_qual=substr($reads_qual{$read_id},$next_start,$next_end-$next_start);
				$next_seq=~tr/ATCG/TAGC/;
				$next_seq=reverse $next_seq;
				$next_qual=reverse $next_qual;
				print $next_seq,"\t",$next_qual,"\n";

				print "Part_",$pet_index,"_",$plus_or_minus_this,"_",$read_id,"\t16\t";
				my $offset=$this_chr_end-$this_chr_start-($this_end-$this_start);
				if($offset == 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start,"M\t*\t0\t0\t";
				}
				elsif($offset > 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start,"M",$offset,"D\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
				}
				my $this_seq=substr($reads_seq{$read_id},$this_start,$this_end-$this_start);
				my $this_qual=substr($reads_qual{$read_id},$this_start,$this_end-$this_start);
				$this_seq=~tr/ATCG/TAGC/;
				$this_seq=reverse $this_seq;
				$this_qual=reverse $this_qual;
				print $this_seq,"\t",$this_qual,"\n";
			}
			else{
				die;
			}

		}
		else{
			#print $read_id,"\tChimeric\n";
			#print $this_chr,"\t",$this_chr_start,"\t",$this_chr_end,"\t",$this_strand,"\n";
			#print $next_chr,"\t",$next_chr_start,"\t",$next_chr_end,"\t",$next_strand,"\n";
			$pet_index++;

			my $this_seq=substr($reads_seq{$read_id},$this_start,$this_end-$this_start);
			my $this_qual=substr($reads_qual{$read_id},$this_start,$this_end-$this_start);
			my $next_seq=substr($reads_seq{$read_id},$next_start,$next_end-$next_start);
			my $next_qual=substr($reads_qual{$read_id},$next_start,$next_end-$next_start);

			if($this_strand eq "+"){
				print "ChimericWhole_",$pet_index,"_",$plus_or_minus_this,"_Head_",$read_id,"\t0\t";
				my $offset=$this_chr_end-$this_chr_start-($this_end-$this_start);
				if($offset==0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start,"M",$next_end-$next_start,"S\t*\t0\t0\t";
				}
				elsif($offset > 0){#deletion
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start,"M",$offset,"D",$next_end-$next_start,"S\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$this_end-$this_start+$offset,"M",0-$offset,"I",$next_end-$next_start,"S\t*\t0\t0\t";
				}	
				print $this_seq,$next_seq,"\t";
				print $this_qual,$next_qual,"\n";
				
			}
			elsif($this_strand eq "-"){
				print "ChimericWhole_",$pet_index,"_",$plus_or_minus_this,"_Head_",$read_id,"\t16\t";
				my $offset=$this_chr_end-$this_chr_start-($this_end-$this_start);
				if($offset==0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$next_end-$next_start,"S",$this_end-$this_start,"M\t*\t0\t0\t";
				}
				elsif($offset > 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$next_end-$next_start,"S",$this_end-$this_start,"M",$offset,"D\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $this_chr,"\t",$this_chr_start,"\t255\t",$next_end-$next_start,"S",$this_end-$this_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
				}
				my $tmp_seq=$this_seq.$next_seq;
				my $tmp_qual=$this_qual.$next_qual;
				$tmp_seq=~tr/ATCG/TAGC/;
				$tmp_seq=reverse $tmp_seq;
				$tmp_qual=reverse $tmp_qual;
				print $tmp_seq,"\t",$tmp_qual,"\n";
			}
			else{
				die;
			}
			if($next_strand eq "+"){
				print "ChimericWhole_",$pet_index,"_",$plus_or_minus_next,"_Tail_",$read_id,"\t0\t";
				my $offset=$next_chr_end-$next_chr_start-($next_end-$next_start);
				if($offset==0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$this_end-$this_start,"S",$next_end-$next_start,"M\t*\t0\t0\t";
				}
				elsif($offset > 0){#deletion
					print $next_chr,"\t",$next_chr_start,"\t255\t",$this_end-$this_start,"S",$next_end-$next_start,"M",$offset,"D\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$this_end-$this_start,"S",$next_end-$next_start+$offset,"M",0-$offset,"I\t*\t0\t0\t";
				}
				print $this_seq,$next_seq,"\t";
				print $this_qual,$next_qual,"\n";
				
			}
			elsif($next_strand eq "-"){
				print "ChimericWhole_",$pet_index,"_",$plus_or_minus_next,"_Tail_",$read_id,"\t16\t";
			 	my $offset=$next_chr_end-$next_chr_start-($next_end-$next_start);
				if($offset==0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start,"M",$this_end-$this_start,"S\t*\t0\t0\t";
				}
				elsif($offset > 0){#deletion
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start,"M",$offset,"D",$this_end-$this_start,"S\t*\t0\t0\t";
				}
				elsif($offset < 0){
					print $next_chr,"\t",$next_chr_start,"\t255\t",$next_end-$next_start+$offset,"M",0-$offset,"I",$this_end-$this_start,"S\t*\t0\t0\t";
				}
                                my $tmp_seq=$this_seq.$next_seq;
                                my $tmp_qual=$this_qual.$next_qual;
                                $tmp_seq=~tr/ATCG/TAGC/;
                                $tmp_seq=reverse $tmp_seq;
				$tmp_qual=reverse $tmp_qual;
                                print $tmp_seq,"\t",$tmp_qual,"\n";
	
			}
			else{
				die;
			}
		}
	}
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









