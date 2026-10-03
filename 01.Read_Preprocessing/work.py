import sys,os

samples=(
"1.A549-1",
"2.A549-2",
"3.HEK293T-1",
"4.HEK293T-2",
)


for i in samples:
	cmd3 = "mkdir %s"%(i)
	cmd4 = "gzip -c -d ../1.rawdata/%s/read1.fq.gz > ./%s/read1.fq"%(i,i)
	cmd5 = "gzip -c -d ../1.rawdata/%s/read2.fq.gz > ./%s/read2.fq"%(i,i)
	cmd6 = "fastqc --extract -o ./%s -t 10 ./%s/read1.fq ./%s/read2.fq"%(i,i,i)
	cmd7 = "java -jar trimmomatic-0.39.jar PE -phred33 -threads 10 -quiet ./%s/read1.fq ./%s/read2.fq ./%s/read1.clean.pair.fq ./%s/read1.clean.unpair.fq ./%s/read2.clean.pair.fq ./%s/read2.clean.unpair.fq ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:7:8:true LEADING:25 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30"%(i,i,i,i,i,i)
	os.system(cmd3)
	os.system(cmd4)
	os.system(cmd5)
	os.system(cmd6)
	os.system(cmd7)

	cmd8 = "perl remove_duplicated_reads.pl ./%s/read1.clean.pair.fq ./%s/read2.clean.pair.fq ./%s/read1.clean.pair.rmDup.fq ./%s/read2.clean.pair.rmDup.fq"%(i,i,i,i)
	cmd9 = "perl remove_duplicated_reads_SE.pl ./%s/read1.clean.unpair.fq ./%s/read1.clean.unpair.rmDup.fq"%(i,i)
	cmd10 = "perl remove_duplicated_reads_SE.pl ./%s/read2.clean.unpair.fq ./%s/read2.clean.unpair.rmDup.fq"%(i,i)
	cmd11 = "cat ./%s/read1.clean.pair.rmDup.fq ./%s/read1.clean.unpair.rmDup.fq > ./%s/read1.clean.rmDup.fq"%(i,i,i)
	cmd12 = "cat ./%s/read2.clean.pair.rmDup.fq ./%s/read2.clean.unpair.rmDup.fq > ./%s/read2.clean.rmDup.fq"%(i,i,i)
	os.system(cmd8)
	os.system(cmd9)
	os.system(cmd10)
	os.system(cmd11)
	os.system(cmd12)

	cmd14 = "cutadapt -j 10 -b \"A{100}\" -b \"C{100}\" -b \"G{100}\" -b \"T{100}\" -n 3 --minimum-length=30 -e 0.1 -o ./%s/read1.clean.rmDup.rmPoly.fq ./%s/read1.clean.rmDup.fq"%(i,i)
	cmd15 = "cutadapt -j 10 -b \"A{100}\" -b \"C{100}\" -b \"G{100}\" -b \"T{100}\" -n 3 --minimum-length=30 -e 0.1 -o ./%s/read2.clean.rmDup.rmPoly.fq ./%s/read2.clean.rmDup.fq"%(i,i)
	cmd16 = "fastqc -t 10 --extract -o ./%s ./%s/read1.clean.rmDup.rmPoly.fq ./%s/read2.clean.rmDup.rmPoly.fq"%(i,i,i)
	os.system(cmd14)
	os.system(cmd15)
	os.system(cmd16)
