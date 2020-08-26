"""File takes a single or set of FASTQ files and counts the times a sequence from a list appears.  Can set a barcode
sequence as well to determine if the sequence is appearing in the right place, as well as a minimum sequence length.
Optional: Use of blast alignment algorithm to determine allow for partial matches."""

import os,re
from Bio import pairwise2

class LibSeqCount:

    def __init__(self):
        self.MySeqList = list()
        self.barcode = str()
        self.readrange = (100, 20000)  # length min and max of the read
        self.results = list()
        self.alignment_ref = "GTTTTAGAGCTAGAAATAGCAAG"

    def import_seq_info(self):
        f = open('/home/trinhlab/Desktop/SaccLibGroup1.txt')
        for line in f:
            self.MySeqList.append(line[:-1].split("\t")[1])

    def import_fastq_data(self, file_directory, barcode, multi=False, use_alignment=False):
        os.chdir(file_directory + barcode)
        filestring = str()
        for fastq_file in os.listdir(os.curdir):
            if fastq_file.endswith(".fastq"):
                f = open(fastq_file)
                for line in f:
                    filestring += line
                f.close()
        f = open(file_directory + "scaffold_alignment_librarycalloutput" + barcode + ".txt", 'w')
        for sequence in self.MySeqList:
            sequence = sequence.upper()
            count = 0
            if use_alignment:
                seqcount = [m.start() for m in re.finditer(sequence, filestring)]
                for match in seqcount:
                    grnaseq = filestring[match+20:match+23]
                    if pairwise2.align.globalxx(grnaseq, self.alignment_ref, score_only=True) > 19:
                        count += 1
                revseqcount = [m.start() for m in re.finditer(self.revcom(sequence), filestring)]
                for match in revseqcount:
                    grnaseq = filestring[match-23:match]
                    if pairwise2.align.globalxx(grnaseq, self.alignment_ref, score_only=True) > 19:
                        count += 1
            else:
                seqcount = re.finditer(sequence, filestring)
                for match in seqcount:
                    count += 1
                revseqcount = re.finditer(self.revcom(sequence), filestring)
                for match in revseqcount:
                    count += 1
            f.write(sequence + "," + str(count) + "\n")
            print("Sequence " + sequence + " completed.")
        f.close()

    def revcom(self, sequence):
        revseq = ""
        change = {'A':'T',
                  'T':'A',
                  'C':'G',
                  'G':'C'}
        for nt in sequence:
            rnt = change[nt]
            revseq = rnt + revseq
        return revseq


L = LibSeqCount()
L.import_seq_info()
barcodes = ["04","05","06"] #"02","03","04","05","06","07","08","09","10","12"]
for num in barcodes:
    L.import_fastq_data("/home/trinhlab/Documents/SaccLibraryAmpliconFastQ/Barcodes/", "barcode" + num, use_alignment=True)
