"""Code for sorting an embedded library of reads (non-Illumina/Oxford barcoded) that are across inherently barcoded sequences (distinct genes that
can be identified by alignment score to a section)"""

class ReadSorter:

    def __init__(self, read_directory, ref_read_file):
        

    def parse_fastqs(self, read_directory):
        for f in read_directory:
            if f.endswith("fastq"):
