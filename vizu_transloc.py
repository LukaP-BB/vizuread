"""
Script de création de graphes de translocation
"""

from dataclasses import dataclass
from pathlib import Path
import logging

from matplotlib import pyplot as plt

from vizuread import Read, get_reads_from, plot_region



def format_big_number(i) : 
    """
    Takes an int or string represnting an int and returns a string with a space separating every 3 numbers
    """
    try : 
        i = int(i)
    except ValueError : 
        raise Exception("This function accepts only ints or strings parsable into an int")

    s = ""
    i = str(i)
    for j, c in enumerate(reversed(i)) : 
        if j%3==0 : 
            s = c + " " + s 
        else : 
            s = c + s
    return s

@dataclass
class Link() :
    c1 : str
    p1 : int
    c2 : str
    p2 : int
    col : str

    def to_string(self) :
        return f"{self.c1}\t{self.p1}\t{self.p1}\t{self.c2}\t{self.p2}\t{self.p2}\tclass={self.col}\n"


@dataclass
class Transloc() :
    f : Path
    c1 : str
    start1 : int
    end1 : int
    c2 : str
    start2 : int
    end2 : int

    def get_reads(self, padding=0, samtools_options="-F 2") :
        reads1 = get_reads_from(self.f, (self.c1, int(self.start1)-padding, int(self.end1)+padding), samtools_command="satmools", samtools_options=samtools_options)
        reads2 = get_reads_from(self.f, (self.c2, int(self.start2)-padding, int(self.end2)+padding), samtools_command="satmools", samtools_options=samtools_options)
        
        reads1 = [r for r in reads1 if r.receiver_chr == self.c2]
        reads2 = [r for r in reads2 if r.receiver_chr == self.c1]

        return reads1, reads2

    def plot(self, save_to=""):
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,5))

        reads1, reads2 = self.get_reads()

        if len(reads1) == 0 or len(reads2) == 0 :
            logging.error(f"Pas assez de reads retrouvés pour {self.f}")
            return

        plot_region(ax=ax[0], reads=reads1)
        plot_region(ax=ax[1], reads=reads2)
        
        ax[0].set_title(f"{self.c1}:{format_big_number(self.start1)}-{format_big_number(self.end1)}")
        ax[1].set_title(f"{self.c2}:{format_big_number(self.start2)}-{format_big_number(self.end2)}")
        plt.suptitle(f"Patient : {self.f.name.split('_')[0]}")
        if save_to == "" :
            plt.show()
            plt.close()
            plt.cla()
        else :
            plt.savefig(save_to)
            plt.close()
            plt.cla()

def links_from_reads(r1:list[Read], r2:list[Read], file_path:Path|str) :
    links = [
        Link(r.chr, r.pos, r.receiver_chr, r.pos_receiver, r.is_forward)
        for r in r1
    ]
    with open(file_path, "w+") as fp :
        for l in links :
            fp.write(l.to_string())

def kario_from_reads(reads:list[Read], file_path):
    c1 = reads[0].chr
    positions1 = [r.pos for r in reads]
    min1, max1 = min(positions1), max(positions1)
    c2 = reads[0].receiver_chr
    positions2 = [r.pos_receiver for r in reads]
    min2, max2 = min(positions2), max(positions2)

    with open(file_path, "w+") as fp :
        fp.write(f"chr - {c1} {c1.replace('chr', '')} {min1-200} {max1+200} {c1}\n")
        fp.write(f"chr - {c2} {c2.replace('chr', '')} {min2-200} {max2+200} {c2}\n")




if __name__ == "__main__" :
    
    nom = "T32943"
    nom = "T32924"
    nom = "T33028"
    nom = "T32886"
    bam = f"/home/luka/Projet/routine/hebdo/{nom}_realigned.fixed.recal.bam"

    reads1 = get_reads_from(bam, position="chr14", samtools_options="-F 2")
    reads1 = [r for r in reads1 if r.receiver_chr == "chr11"]
    
    reads2 = get_reads_from(bam, position="chr11", samtools_options="-F 2")
    reads2 = [r for r in reads2 if r.receiver_chr == "chr14"]

    fig, axes = plt.subplots(1, 2, figsize=(10,5))

    plot_region(reads=reads1, ax=axes[0])
    plot_region(reads=reads2, ax=axes[1])

    plt.show()

    # links_from_reads(reads1, reads2, Path(f"circos/links/{nom}.links.txt"))                
    # kario_from_reads(reads1, Path(f"circos/links/{nom}.links.kario"))
