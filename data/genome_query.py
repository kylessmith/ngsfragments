import genomepy
from ngsfragments import fragments
from collections import Counter


frags = fragments("/Users/ksmith10/Desktop/CSF/AL328_Swift_43222_3_18_2015dup.bam", verbose=True)

g = genomepy.Genome("hg19", genome_dir="./hg19")

k = 2
kmers = {
    "AA":0,
    "AT":0,
    "AG":0,
    "AC":0,
    "TA":0,
    "TT":0,
    "TG":0,
    "TC":0,
    "GA":0,
    "GT":0,
    "GG":0,
    "GC":0,
    "CA":0,
    "CT":0,
    "CG":0,
    "CC":0
}
kmers = Counter()
for chrom in peaks:
    print(chrom)
    for i in peaks[chrom]:
        s = g[chrom][i.start:i.end]
        for i in range(len(s.seq)-k):
            kmers[s.seq[i:i+k].upper()] += 1