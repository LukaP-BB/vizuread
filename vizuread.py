import re
from statistics import mean
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import cm


class ReadException(Exception):
    """base exception for the read errors"""

class MalformedEntry(ReadException):
    """raised when a read doesn't contain the expected fields"""

class EmptyEntry(ReadException):
    """raised when a blank line is passed in place of a read"""

class BadCigar(ReadException) :
    """raised when the cigar string could not be parsed"""

MAPQ_COLORS_PROPERLY_PAIRED = cm.get_cmap("Blues", 256)
MAPQ_COLORS = cm.get_cmap("Oranges", 256)
READ_SPACING = 20

def get_mean_qual(seq: str):
    """returns the mean quality of a phred Q string"""
    scores = [ord(c)-33 for c in seq]
    return round(mean(scores), 1)

def is_forward(flags:int):
    """returns if a read is forward or reverse based on the presence of the x10 flag in the SAM flags"""
    s = [ # crée une liste de 12 ints (0 ou 1) qui correspondent chacun à un flag
        int(_) for _ in 
        list(reversed(format(flags, "b").rjust(12, "0")))
    ]
    # le 5ème bit correspond au flag "read reverse strand"
    return s[4] == 0

def parse_cigar(cigar:str) :
    """returns a list of segments to be plotted as arrows"""
    matches = re.findall(r"(\d+)(\w)", cigar)
    if len(matches) == 0 and cigar != "*" :
        raise BadCigar(f"cigar '{cigar}' could not be parsed")
    
    adds_to_ref_span = {
        "M" : True,
        "I" : False,
        "D" : True,
        "N" : True,
        "S" : False,
        "H" : False,
    }
    segments = []
    reff_span = 0
    plot_length = 0
    for i, match in enumerate(matches) :
        length, operation = int(match[0]), match[1]

        if adds_to_ref_span[operation] :
            reff_span += length

        if operation != "I" :
            plot_length += length

        segments.append((operation, length))
    return segments, reff_span, plot_length

class Read():
    """
    Class for a read extracted with samtools
    """
    def __init__(self, e:str) -> None:
        """
        Initialize a Read object from a line obtained with the command `samtools view file.bam | cut -f 2,3,4,5,6,7,8,10,11`
        """
        fields = e.strip().split()
        
        expected = 9
        if len(fields) == 0 :
            raise EmptyEntry()
        if len(fields) != expected :
            raise MalformedEntry(f"Unexpected number of fields in a read (got {len(fields)}, expected {expected}) : \n{e}")

        try :
            self.flag = int(fields[0])
            self.chr = fields[1]
            self.pos = int(fields[2])
            self.mapQ = int(fields[3])
            self.cigar = fields[4]
            self.receiver_chr = fields[5]
            self.pos_receiver = int(fields[6])
            self.seq = fields[7]
            self.qual = fields[8]
        except ValueError :
            # happens if one str could not be converted to int
            infos = f"flag={fields[0]}, pos={fields[2]}, mapQ={fields[3]}, pos_mate={fields[6]}"
            raise MalformedEntry(
                f"Could not convert one of the following to an integer : {infos}\nRead : '{e}'"
            )

        self.length = len(self.seq)
        self.mean_qual = get_mean_qual(self.qual)
        self.is_forward = is_forward(self.flag)
        self.is_properly_paired = (self.receiver_chr == "=")

        self.segments, self.ref_span, self.plot_len = parse_cigar(self.cigar)

        if self.is_forward :
            self.start = self.pos
            self.end = self.pos+self.length
        else :
            self.start = self.pos-self.length
            self.end = self.pos


    def plot(self, ax:plt.axes, ypos:int, **kwargs) :
        """
        Plots an arrow representing a read on a plt.ax at the given y position.
        
        Kwargs are passed to plt.arrow https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.arrow.html 
        """
        if self.is_properly_paired :
            color = MAPQ_COLORS_PROPERLY_PAIRED(self.mapQ)
        else :
            color = MAPQ_COLORS(self.mapQ)

        cursor_pos = self.pos

        # defining colors and widths for the different types of segments
        colors = {"H" : "black", "S" : "grey", "D" : "red", "I" : "green"}
        widths = {"D" : 0.3, "I" : 0.9}

        for i, segment in enumerate(self.segments) :
            operation, length = segment
            # we use the MAPQ color only for "M" segments, aka matches
            # for other segments, we use defined colors
            s_color = colors.get(operation, color)
            # same logic for the arrow width
            width = widths.get(operation, 0.8)
            head_width = width if operation != "D" else 0.3

            # for insertions, we make sure they appear on top and we center them
            if operation == "I" : 
                z_order = 2 
                drift = (length/2) if self.is_forward else (-length/2)
            else : 
                z_order = 1
                drift = 0

            if self.is_forward :
                head_length = 4 if i+1 == len(self.segments) else 0 # tracing the arrow head only for the tip of the read
                ax.arrow(
                    x=cursor_pos+drift, y=ypos, dx=length, 
                    dy=0, width=width, head_width=head_width, 
                    head_length=head_length, length_includes_head=True, 
                    color=s_color, zorder=z_order,
                    **kwargs
                )
            else :
                head_length = 4 if i==0 else 0 # tracing the arrow head only for the tip of the read
                ax.arrow(
                    x=cursor_pos+length+drift, y=ypos, dx=-length, 
                    dy=0, width=width, head_width=head_width, 
                    length_includes_head=True, head_length=head_length, 
                    color=s_color, zorder=z_order,
                    **kwargs
                )   

            # for insertions, the reference span stay the same
            if operation != "I" : cursor_pos += length 

    def overlap(self, o:"Read") :
        """
        Returns if 2 reads are overlapping with each other
        """
        if o.start <= self.end <= o.end :
            return True
        if self.start <= o.end <= self.end :
            return True

        return False

    def __str__(self) -> str:
        return f"{self.chr}:{self.start}-{self.end} {self.cigar}"

def parse_position(pos:str) :
    """
    Parse a genomic position from a string. Chromosomes can be named "chr" or simply called by number
    Comas, tabs and whitespaces are deleted, allowing input strings like :
    ```verb
    - chr11:36270167-36270242
    - chr11:36,270,167-36,270,242
    - 11:36 270 167 - 36 270 242
    - ...
    """
    pos = pos.replace(",", "").replace(" ", "").replace("\t", "")
    regex = r"(?P<chrom>(chr)?(\d+|X|Y)):(?P<start>\d+)-(?P<end>\d+)"
    match = re.match(regex, pos)
    try :
        return match.group("chrom"), match.group("start"), match.group("end")
    except AttributeError :
        raise ValueError(f"Could not parse the string '{pos}' into valid positions with the regex '{regex}'")


def get_reads_from(bam_file, position, samtools_command="samtools", samtools_options=""):
    """
    Returns a generator of reads from a given region. 
    If samtools isn't in your path, you can overwrite the default samtools_command kwarg by an appropriate one.

    You can use the helper function parse_position() to input a string similar to what you could give to IGV or UCSC genome browser.

    You can use the samtools_options kwarg to specify filtering options to the `samtools view` command, eg 
    ```py
    position = parse_position("chr11:36,270,167-36 270 242")
    reads = get_reads_from(f, *position, samtools_command="samtools")

    # OR convert to list for multiple iterations over the reads : 
    list_of_reads = list(get_reads_from(...))
    last_read = list_of_reads[-1]
    ```
    """
    if isinstance(position, tuple) or isinstance(position, list) :
        try :
            chrom, start, end = position
        except Exception :
            raise Exception("The tuple or list passed as the position argument must have 3 elements, no more, no less. Three shall be the number thou shalt count, and the number of the counting shall be three. Four shalt thou not count, neither count thou two, excepting that thou then proceed to three. Five is right out")
    elif isinstance(position, str) :
        chrom, start, end = parse_position(position)
    else :
        raise Exception("position argument expected either a tuple or a string")

    samtools = f"{samtools_command} view {samtools_options} {bam_file} {chrom}:{start}-{end}"
    print(samtools)
    sam = sp.Popen(samtools.split() , stdout=sp.PIPE, stderr=sp.PIPE)
    cut_command = "cut -f 2,3,4,5,6,7,8,10,11"
    cut = sp.Popen(cut_command.split(), stdin=sam.stdout, stdout=sp.PIPE)
    for line in cut.stdout.readlines() :
        yield Read(line.decode()) 


def plot_region(
    bam_file:str=None, region:str=None, ax:plt.Axes=None, 
    reads=[], samtools_command="samtools", samtools_options="", 
    piling="spaced", **kwargs) :
    """
    Plots reads from a specific region on a matplotlib ax. Returns the list of Read objects.
    
    By default, the reads will be piled up with an algorithm trying to optimize the space they use on the plot.
    This algorithm can be changed with the `piling` kwarg and accepts :
    - "compact" : the reads will take the minimum amount of space possible
    - "spaced"  : same as compact, but with some spacing on the left and the right of each read
    - "seq"     : the reads are placed on the bottom of the graph again only if there is a break between the reads
    - None      : each read corresponds to a line on the graph

    A list of reads can be directly passed. In that case, every other arguments except `ax` will be ignored.
    This is useful if you wanna retrieve a list of reads and perform custom operations on them before plotting them. 

    If samtools isn't in your path, you can overwrite the default samtools_command kwarg by an appropriate one.

    You can filter reads by passing additional options to the samtools command. See `samtools view --help` for all available options. 
    Options that make samtools write it's output to a file instead of stdout are to avoid.

    Flags explanation for the -f, -F and -G options : https://broadinstitute.github.io/picard/explain-flags.html

    Additional kwargs are passed to plt.arrow https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.arrow.html 

    ```
    # create sublots
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,5))
    # define regions
    region1 = (chr1, 10000, 20000)      # define a tuple
    region2 = "chr2:4587639-5789456"    # or use a string à la IGV
    # plot the regions on the given axes
    plot_region(bam_file, region1, ax[0])
    plot_region(bam_file, region2, ax[1], samtools_options="-F 2") # use the flag -F 2 to exclude all reads that are properly paired.
    # do some additional manipulations to your axes
    ax[1].set_title("some title")
    plt.show()
    """

    if ax is None : raise Exception("ax must be defined")

    if reads == [] :
        if bam_file is None : raise Exception(f"bam_file must be defined")

        # consuming the generator into a list so we can return it
        reads = list(get_reads_from(bam_file, region, samtools_command=samtools_command, samtools_options=samtools_options))
    
    if piling is None :
        for i, r in enumerate(reads) :
            r.plot(ax, i, **kwargs)

    elif piling in {"compact", "spaced"} :
        if piling == "spaced" :
            padding = READ_SPACING
        else :
            padding = 0
        rightmosts = [] # this will keep track of rightmost positions
        i = 0
        for r in reads :
            for j, right_pos in enumerate(rightmosts) :
                if right_pos + padding < r.start :
                    rightmosts[j] = r.start+r.plot_len
                    r.plot(ax, j)
                    i = j
                    break
                else :
                    i = j+1
            else :
                # if the loop didn't break
                rightmosts.append(r.start+r.plot_len)
                r.plot(ax, i)
                        
    elif piling == "seq" :
        rightmost = 0
        i = 0
        for r in reads :
            if r.start > rightmost :
                i = 0
                r.plot(ax, i)
            else :
                r.plot(ax, i)

            if r.start + r.plot_len > rightmost :
                rightmost = r.start + r.plot_len
            i += 1

    else :
        raise Exception("piling argument has to be one of [None, 'compact', 'seq', 'spaced']")        

    return reads

def plot_transloc(
    f = "T30989_realigned.fixed.recal.bam",
    c1 = "chr11",
    c2 = "chr14",
    s = "69,638,162",
    e = "69,639,433"
    ) :

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,5))
    reads = plot_region(f, (c1, s, e), ax[0], samtools_options="-F 2")
    p = [r.pos_receiver for r in reads if r.receiver_chr == c2]
    mi, ma = min(p), max(p)

    reads2 = get_reads_from(f, c2, mi, ma, flags="-F 2")
    plot_region(ax=ax[1], reads=reads2)
    ax[0].set_title(f"Reads {c1}")
    ax[1].set_title(f"Reads {c2}")    
    plt.suptitle(f"Reads attestant d'une translocation {c1}-{c2}")
    plt.show()
