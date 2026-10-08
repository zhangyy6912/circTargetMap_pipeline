#!/usr/bin/perl
die "perl $0 result1.circRNA_linked_target.bed result2.target_to_ele.overlap\n" if(@ARGV != 2);
my $target_bed=shift;
my $target_ele_overlap=shift;

my %target_anno;
open(TEO,$target_ele_overlap) || die;
while(my $line=<TEO>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $target_strand_symbol=(split/:/,$sub[5])[-1];
	my $target_strand;
	if($target_strand_symbol eq "Plus"){
		$target_strand="+";
	}
	elsif($target_strand_symbol eq "Minus"){
		$target_strand="-";
	}
	else{
		die;
	}

	my $index=$sub[7];

	my $gene_symbol=(split/_/,$sub[12],3)[-1];
	my $ele;
	if($gene_symbol =~ /Intergenic/){
		$ele="Intergenic";
		$target_anno{$index}{"Intergenic"}=1;
	}
	else{
		my $sense_or_anti;
		if($target_strand eq $sub[14]){
			$sense_or_anti="sense";
		}
		else{
			$sense_or_anti="anti";
		}			
		
		$ele=$gene_symbol."#".$sub[13];
		$target_anno{$index}{$sense_or_anti}{$ele}=1;
	}
}

open(TB,$target_bed) || die;
while(my $line=<TB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my $index=$sub[7];
	if(exists $target_anno{$index}){
		if(exists $target_anno{$index}{"sense"}){
			my @all_candidate=keys %{$target_anno{$index}{"sense"}};
			$gene_attr=join"|",@all_candidate;
			print $line,"\t",$gene_attr,"\tOnTheSenseStrand\n";
		}
		elsif(exists $target_anno{$index}{"anti"}){
			my @all_candidate=keys %{$target_anno{$index}{"anti"}};
			$gene_attr=join"|",@all_candidate;
			print $line,"\t",$gene_attr,"\tOnTheAntiStrand\n";
		}
		elsif(exists $target_anno{$index}{"Intergenic"}){
			print $line,"\tIntergenic\tIntergenic\n";
		}
		else{
			die;
		}
	}
	else{
		print $line,"\tUnknown\tUnknown\n";
	}
}
	
	

