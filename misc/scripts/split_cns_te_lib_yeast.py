import argparse
import pysam
import re

# from ragtag_utilities.utilities import reverse_complement

"""
Like bedtools getfasta, but use the gff ID attribute as the FASTA header and
always force strandedness. 
"""


def main():
    parser = argparse.ArgumentParser(description="Get fasta sequences from a GFF file")
    parser.add_argument("fasta", metavar="<sequences.fasta>", type=str, help="The full-length consensus TE library.")
    parser.add_argument("gff", metavar="<genes.gff>", type=str, help="The gff file for full-length TE library.")
    
    args = parser.parse_args()
    fasta_file = args.fasta
    gff_file = args.gff
    
    inputseq = pysam.FastaFile(fasta_file)
    
    superfamily = {}
    internal_region = {}
    ltr_region = {}

    # Read the gff file
    with open(gff_file, "r") as f:
        for line in f:
            if line.startswith("##sequence-region"):
                fmt, tefamily, start, end = line.rstrip().split()
                if re.match(r"ty3", tefamily.lower()):
                    superfamily[tefamily] = "Gypsy"
                else:
                    superfamily[tefamily] = "Copia"
                internal_region[tefamily] = {"start": None, "end": None}
                ltr_region[tefamily] = {"start": None, "end": None}

            if not line.startswith("#"):
                tefamily, source, feature, start, end, score, strand, fname, attributes = line.rstrip().split("\t")
                start, end = int(start)-1, int(end)
                
                #     raise ValueError("Need an ID attribute for each gff line.")
                if (feature == "long_terminal_repeat") and (start == 0):
                    internal_region[tefamily]["start"] = end
                    ltr_region[tefamily]["start"] = start
                    ltr_region[tefamily]["end"] = end

                if (feature == "long_terminal_repeat") and (start > 0):
                    internal_region[tefamily]["end"] = start

    # print internal coding region
    for i in internal_region.keys():
        repeatmasker_id = i + "-LTR#LTR/" + superfamily[i]
        print(">" + repeatmasker_id)
        print(inputseq.fetch(i, ltr_region[i]["start"], ltr_region[i]["end"]))
        repeatmasker_id = i + "-I#LTR/" + superfamily[i]
        print(">" + repeatmasker_id)
        print(inputseq.fetch(i, internal_region[i]["start"], internal_region[i]["end"]))
if __name__ == "__main__":
    main()
