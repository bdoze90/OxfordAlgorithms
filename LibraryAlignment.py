"""File for identifying short sequences from Nanopore library sequencing data."""

from Bio import pairwise2

alignment = pairwise2.align.localxx("GACTACCG","GACTTACCGTAGCCCAGTAGTAGACTTACCGTGGATAGACTAGTCGACTTACCGTATTTTGACTACCG")
print(alignment)
