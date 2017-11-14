import subprocess
import pysam
import sys
import pandas
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pylab
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("bam1")
parser.add_argument("bam2")
parser.add_argument("-n","--name",help="name of the output plot")
parser.add_argument("-s","--species",help="currently onlhy 'huma' and 'dog' are recognized",default='human',choices=['human','dog'])
args = parser.parse_args()

#print(args.name)

if args.species == 'human':
    COORDINATES = "./snp/all.snps.sorted.txt"
elif args.species == 'dog':
    COORDINATES = "./snp/all.dog_snps_canfam1.0.bed.txt"

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



#read bed file

coord = pandas.read_csv(COORDINATES, sep='\t',header=None)
index = coord[4].tolist()

columns=[args.bam1, args.bam2]
df = pandas.DataFrame(index=index,columns=columns)

full = 2*len(coord.index)
inc = 0


printProgressBar(inc, full, prefix='Progress:', suffix='Complete', length=50)
for file in [args.bam1,args.bam2]:
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
            df.loc[[row[4]],[file]] = (all-ref)/all
        else:
            df.loc[[row[4]], [file]] = -0.1

        #add column to dataframe



#plot graph
b1 = args.bam1.split("/")
b2 = args.bam2.split("/")


if(args.name != None):
    plt.title(args.name)
plt.xlabel(b1[len(b1)-1].replace(".bam", ""))
plt.ylabel(b2[len(b2)-1].replace(".bam", ""))
plt.plot(df[args.bam1].tolist(), df[args.bam2].tolist(), "o", mfc='None',)
if(args.name != None):
    fn = args.name
else:
    fn = 'out'
pylab.savefig(fn+".pdf")