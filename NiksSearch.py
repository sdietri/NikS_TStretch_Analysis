import re
import seaborn as sns
import pandas as pd
import glob

from subprocess import call
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

class searchRes:
    def __init__(self):
        self.SRR = ""
        self.info = ""
        self.hitSeq = ""
        self.nikS_seq = ""
        self.promreg = ""
        self.T_stretchLen = 0
    

def searchSRR(SRR, fastaDirectory, nikSstartFasta):
    myRes = searchRes()
    handle = Entrez.efetch(db="sra", id=SRR)
    Info = handle.read()
    m = re.search("<TITLE>(.+)</TITLE>", Info)
    myRes.info = m.group(1)
    myRes.SRR = SRR
    fasta = fastaDirectory + SRR + ".fas"
    f = open(fasta,"r")
    lines = f.readlines()
    seqs = {}
    seqName = ""
    seq = ""
    for line in lines:
        if line.startswith(">"):
            if not seqName == "":
                seqs[seqName] = seq;
                seq = ""
            seqName = line.split(" ")[0].strip().replace(">","")
        else:
            seq = seq + line.strip()
    seqs[seqName] = seq
    call(["makeblastdb", "-in", fasta, "-out", "blastdb", "-dbtype","nucl"])
    call(["blastn", "-query", nikSstartFasta, "-db", "blastdb", "-outfmt","6","-out","blast.xml"])
    #parse blastres
    f = open("blast.xml","r")
    lines = f.readlines()
    for line in lines:
        splitL = line.split("\t")
        if splitL[6] == "1" and int(splitL[3])>90:
            start = int(splitL[8])-1
            stop = int(splitL[9])
            if stop < start:
                start = len(seqs[splitL[1]]) - start
                stop = len(seqs[splitL[1]])- stop
                my_dna = Seq(seqs[splitL[1]], generic_dna)
                seq = str(my_dna.reverse_complement())
                if start < 40:
                    continue
                nikS_seq = seq[start:stop]
                promreg = seq[start-40:start]
                m = re.search("[T]{6,}", promreg)
                if not m == None:
                    myRes.nikS_seq = nikS_seq
                    myRes.promreg = promreg
                    myRes.hitSeq = splitL[1]
                    myRes.T_stretchLen = len(m.group(0))
                    break
                    
            else: 
                if start < 40:
                    continue
                nikS_seq = seqs[splitL[1]][start:stop]
                promreg = seqs[splitL[1]][start-40:start]
                m = re.search("[T]{6,}", promreg)
                if not m == None:
                    myRes.nikS_seq = nikS_seq
                    myRes.promreg = promreg
                    myRes.hitSeq = splitL[1]
                    myRes.T_stretchLen = len(m.group(0))
                    break
    return myRes


def main(args):
    Entrez.email = args.email
    
    hits = []
    for name in glob.glob(args.sraFastaDir + 'SRR*.fas'):
        hits.append(searchSRR(name.split("/")[-1].replace(".fas",""), args.sraFastaDir, args.genestartfasta))

    o = open("NikS_col.csv","w")
    o.write("SRR\tDescription\tTop hit Sequence\tNikS sequence\tPromoter Sequence\tT-Stretch length\n")
    for curRes in hits:
        info = [curRes.SRR, curRes.info, curRes.hitSeq, curRes.nikS_seq, curRes.promreg, str(curRes.T_stretchLen)]
        o.write("\t".join(info) + "\n")
    o.flush()
    o.close()

    f = open("NikS_col.csv", "r")
    patients = []
    counts = {}
    lines = f.readlines()
    for line in lines[1:]:
        splitL = line.split("\t")
        curPat = splitL[1].split(":")[1].split(",")[0][1:]
        locus = splitL[1].split(":")[1].split(",")[1].replace(" ","")  
        if "colony" in locus:
            locus = locus.split("colony")[0]   
        locus = locus[0]   
        Tlen = int(splitL[5])
        if Tlen == 0:
            continue
        if not curPat in patients:
            patients.append(curPat)
            counts[curPat] = {"c":{}, "a":{},"f":{}}
        if not Tlen in counts[curPat][locus]:
            counts[curPat][locus][Tlen]=1
        else:
            counts[curPat][locus][Tlen]+=1

    pats = []
    locus = []
    lens = []
    TLcounts = []
    for curPat in patients:
        for curLocus in counts[curPat]:
            for curLen in counts[curPat][curLocus]:
                curCount = counts[curPat][curLocus][curLen]
                pats.append(curPat)
                locus.append(curLocus)
                lens.append(curLen)
                TLcounts.append(curCount)
    df = pd.DataFrame(list(zip(pats, locus, lens, TLcounts)),
                     columns=['Patient', 'Locus', 'Tlen', 'Count'])
    df = df.sort_values(by=['Patient'])

    sns.set_style("whitegrid")
    #pal = sns.cubehelix_palette(12, start=.5, rot=-.75)
    pal = reversed(sns.color_palette("RdYlGn", 12))
    g = sns.catplot(data=df, x='Locus',y='Tlen',col="Patient", hue='Count',palette=pal,
                    height=5,aspect=.15,kind="strip", dodge=False, jitter=False,sharey=True,
                    order=['c','a','f'])
    axes = g.axes.flatten()
    for i in range(len(axes)):
        curTitle = axes[i].title.get_text()
        curTitle = curTitle.split("=")[1].replace(" patient ","").replace(" (timepoint 2)","_2")
        axes[i].set_title(curTitle)
        axes[i].set_ylim([0,27])
    g.fig.tight_layout()
    g.set_axis_labels(x_var="", y_var="T Stretch length")
    g.despine(left=False)

    g.fig.subplots_adjust(wspace=0, hspace=0.2)
    g._legend.set_bbox_to_anchor((1.04, 0.5, 0, 0))

    g.savefig("NikS_dot.png")

    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("email", help="Email to use as access identification for Entrez services")
    parser.add_argument("sraFastaDir", help="Directory containing the fasta converted SRR reads")
    parser.add_argument("genestartfasta", help="FASTA with the gene start sequence")
    args = parser.parse_args()
    
    main(args)
