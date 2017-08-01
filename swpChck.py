import subprocess

proc = subprocess.run(["which", "samtools"], stdout=subprocess.PIPE)

if(len(proc.stdout)==0):
    print("Samtools not in path. Exiting")
    exit(1)
else:
    print("Found Samtools")


#check if bam files have correct coordinate

#read bam files (what number?)

#check if chr or not chr

#samtools mpileup (in memory?)

#plot pairs
