import subprocess
import pysam
import sys
import pandas

COORDINATES = "./SNP/all.snps.sorted.bed"

#proc = subprocess.run(["which", "samtools"], stdout=subprocess.PIPE)

#if(len(proc.stdout)==0):
#    print("Samtools not in path. Exiting")
#    exit(1)
#else:
#    print("Found Samtools")


if len(sys.argv) <= 2:
    print("please provide at least two exome bam files as arguments.")
    exit(2)
else:
    print("bam provided " + str(len(sys.argv[1:])) + " bam files: ")
    print('\n'.join('{}: {}'.format(*k) for k in enumerate(sys.argv[1:])))

#read bed file

coord = pandas.read_csv(COORDINATES, sep='\t',header=None)
out = coord[[3]]


out.columns=['SNP']
print(out)

#create dataframe from bed file


for file in sys.argv[1:]:
    samfile = pysam.Samfile(file, "rb")
    if 'chr' in samfile.references[1]:
        print("hg19 assumed: chr detected in the reference of ", file)
        chr = 'chr'
    else:
        print("GRCh38 assumed: chr not detected in reference of ", file)
        chr = ''

    for index,row in coord.iterrows():
        #print(row)

        for pileupcolumn in samfile.pileup(chr+str(row[0]), row[1], row[2]):
            for pileupread in pileupcolumn.pileups:
                if pileupcolumn.pos == row[1]:
                    #print(row[3])
                    #print("query pos")
                    if pileupread.query_position != None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                    #print(base)

        #add column to dataframe



#plot graph
