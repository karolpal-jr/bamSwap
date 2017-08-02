import subprocess
import pysam
import sys
import pandas
import matplotlib.pyplot as plt
import pylab

COORDINATES = "./snp/all.snps.sorted.txt"

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total:
        print()



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
index = coord[4].tolist()


columns=[sys.argv[1:]]
df = pandas.DataFrame(index=index,columns=columns)

full = len(df)*len(sys.argv[1:])
inc = 0


printProgressBar(inc, full, prefix='Progress:', suffix='Complete', length=50)
for file in sys.argv[1:]:
    samfile = pysam.Samfile(file, "rb")
    if 'chr' in samfile.references[1]:
        print("hg19 assumed: chr detected in the reference of ", file)
        chr = 'chr'
    else:
        print("GRCh37 assumed: chr not detected in reference of ", file)
        chr = ''

    for index,row in coord.iterrows():
        inc = inc + 1
        #print(row)
        ref = 0
        all = 0
        for pileupcolumn in samfile.pileup(chr+str(row[0]), row[1], row[2]):
            for pileupread in pileupcolumn.pileups:
                if pileupcolumn.pos == row[1]: # 0 based coordinates ??
                    if pileupread.query_position != None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        all = all + 1
                        if row[3]== base.upper():
                            ref = ref + 1
                        #print(base)
        printProgressBar(inc, full, prefix='Progress:', suffix='Complete', length=50)
        if all > 0:
            df.loc[[row[4]],[file]] = ref/all
        else:
            df.loc[[row[4]], [file]] = 0

        #add column to dataframe



#plot graph

plt.plot(df[sys.argv[1]].tolist(), df[sys.argv[2]].tolist(),"o")
pylab.savefig("out.pdf")