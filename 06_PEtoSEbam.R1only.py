#PURPOSE: change bam flags from PE to SE for downsampling and counting
#Based on this script: https://bioinformatics.stackexchange.com/questions/9228/convert-paired-end-bam-into-a-single-end-bam-and-keep-all-the-reads

import sys, pysam, random
#NOTE: using an old pysam (0.7.6) so some variable names have changed; see https://pysam.readthedocs.io/en/latest/release.html Release 0.8.1

if(len(sys.argv) != 2):
        sys.exit("Enter name of merged bamfile to convert to SE and remove repeated reads from")

print "Reading in file..."
outname = '.'.join((sys.argv[1]).strip().split(".")[0:4])
print outname

paired =         1
proper_pair =    2
mate_unmapped =  8
mate_reverse =   32
first_in_pair =  64
second_in_pair = 128
unpaired_f = 0
unpaired_r = 16

bam_in = pysam.Samfile(sys.argv[1], "rb")
bam_out = pysam.Samfile(outname+".SE_R1orMateUnmappedOnly.bam", "wb", template=bam_in)
#bam_in = pysam.Samfile("LLWW_SRR2060837.trimmomatic.DomPwkCombo.DOMallelesONLY.bam", "rb")
#bam_out = pysam.Samfile("LLWW_SRR2060837.trimmomatic.DomPwkCombo.DOMallelesONLY.SE.bam", "wb", template=bam_in)
#bam_in = pysam.Samfile("temp.bam", "rb")
#bam_out = pysam.Samfile("temp.SE.bam", "wb", template=bam_in)
for line in bam_in:
    #Only keep R1, convert to SE flags
#    if line.flag & (first_in_pair | mate_unmapped | unpaired_f | unpaired_r):
     if (first_in_pair | mate_unmapped | unpaired_f | unpaired_r):
	if(unpaired_f):
		print line
        line.pos += line.tlen
	line.rnext = 0
        line.pnext = 0
        line.flag &= ~(paired + proper_pair + mate_unmapped + mate_reverse + first_in_pair + second_in_pair)
	bam_out.write(line)

bam_in.close()
bam_out.close()
