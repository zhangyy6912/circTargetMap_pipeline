#!/usr/bin/python3
import sys,os

work_path = os.getcwd()
# script_path = "/data2/zhangyy/01.RIC-seq-tools/3.analysis"

species = "human"
ref_version = "hg38"

rRNA_STAR_index_dir = "/data0/database/1.human/rRNA_reference_index"
rRNA_genome_ref_fasta = "/data0/database/1.human/rRNA_reference_index/rRNA_reference.fa"
chrM_STAR_index_dir = "/data0/database/1.human/chrM_reference_index"
chrM_genome_ref_fasta = "/data0/database/1.human/chrM_reference_index/chrM.fa"
genome_STAR_index_dir = "/data0/database/1.human/GRCh38.p13_cleangenome_index"
genome_genome_ref_fasta = "/data0/database/1.human/GRCh38.p13_cleangenome_index/hg38.clean.fa"

threads = 60  # default if server is fully available, use 100 or more.

samples = (
"85.32-RIC-1",
"86.32-RIC-2",

)


###change the number to use the following step
match5 = 1



### 5.split_rRNA_and_otherReads
if match5:
    os.system("mkdir 5.split_rRNA_and_otherReads")
    os.system(f"cd {work_path}/5.split_rRNA_and_otherReads")
    outfile = open(rf"{work_path}/5.split_rRNA_and_otherReads/work.sh", "w")
    for i in samples:
        os.system(f"mkdir {work_path}/5.split_rRNA_and_otherReads/{i}")
        os.system(f"ln -s ../../../2.trim/{i}/read1.clean.rmDup.rmPoly.fq {work_path}/5.split_rRNA_and_otherReads/{i}")
        os.system(f"ln -s ../../../2.trim/{i}/read2.clean.rmDup.rmPoly.fq {work_path}/5.split_rRNA_and_otherReads/{i}")

        print(f"STAR --runMode alignReads --genomeDir {rRNA_STAR_index_dir} --readFilesIn ./{i}/read1.clean.rmDup.rmPoly.fq --outFileNamePrefix ./{i}/{i}_read1_torRNA_ --outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN {threads} --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5", file=outfile)
        print(f"STAR --runMode alignReads --genomeDir {rRNA_STAR_index_dir} --readFilesIn ./{i}/read2.clean.rmDup.rmPoly.fq --outFileNamePrefix ./{i}/{i}_read2_torRNA_ --outReadsUnmapped Fastx --outFilterMultimapNmax 100 --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 --limitIObufferSize 1500000000 --runThreadN {threads} --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5", file=outfile)
        print(f"perl select_not_fully_aligned.pl ./{i}/read1.clean.rmDup.rmPoly.fq ./{i}/{i}_read1_torRNA_Aligned.out.sam ./{i}/{i}_read1_torRNA_Chimeric.out.sam > ./{i}/{i}_read1_torRNA_Unmapped_really.fq", file=outfile)
        print(f"perl select_not_fully_aligned.pl ./{i}/read2.clean.rmDup.rmPoly.fq ./{i}/{i}_read2_torRNA_Aligned.out.sam ./{i}/{i}_read2_torRNA_Chimeric.out.sam > ./{i}/{i}_read2_torRNA_Unmapped_really.fq", file=outfile)

        os.system(f"mkdir {work_path}/5.split_rRNA_and_otherReads/{i}/z1.second_round_byBWA")
        print(f"bwa mem -t {threads} -k 12 -T 15 -o ./{i}/z1.second_round_byBWA/read1_futher_by_bwa.sam {rRNA_genome_ref_fasta} ./{i}/{i}_read1_torRNA_Unmapped_really.fq", file=outfile)
        print(f"bwa mem -t {threads} -k 12 -T 15 -o ./{i}/z1.second_round_byBWA/read2_futher_by_bwa.sam {rRNA_genome_ref_fasta} ./{i}/{i}_read2_torRNA_Unmapped_really.fq", file=outfile)
        print(f"perl collect_chimeric_ligation_from_sam.pl ./{i}/z1.second_round_byBWA/read1_futher_by_bwa.sam 1 ./{i}/z1.second_round_byBWA > ./{i}/z1.second_round_byBWA/out1.read1.chimeric.sam", file=outfile)
        print(f"perl collect_chimeric_ligation_from_sam.pl ./{i}/z1.second_round_byBWA/read2_futher_by_bwa.sam 2 ./{i}/z1.second_round_byBWA > ./{i}/z1.second_round_byBWA/out1.read2.chimeric.sam", file=outfile)
        print(f"perl get_final_unmapped_reads.pl ./{i}/{i}_read1_torRNA_Unmapped_really.fq ./{i}/z1.second_round_byBWA/unmapped.1.readID.list > ./{i}/z1.second_round_byBWA/read1_torRNA_Unmapped_really.FutherByBwa.fq", file=outfile)
        print(f"perl get_final_unmapped_reads.pl ./{i}/{i}_read2_torRNA_Unmapped_really.fq ./{i}/z1.second_round_byBWA/unmapped.2.readID.list > ./{i}/z1.second_round_byBWA/read2_torRNA_Unmapped_really.FutherByBwa.fq", file=outfile)

    print(f"cd {work_path}/5.split_rRNA_and_otherReads")
    print("bash work.sh")

