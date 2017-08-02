import subprocess
import pysam
import sys

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
    print("bam files provided: ", sys.argv[1:])


#check if chr or not chr

#samtools mpileup (in memory?)

#plot pairs
