"""File takes a single or set of FASTQ files and counts the times a sequence from a list appears.  Can set a barcode
sequence as well to determine if the sequence is appearing in the right place, as well as a minimum sequence length.
Optional: Use of blast alignment algorithm to determine allow for partial matches."""

import os,re
from Bio import pairwise2
import multiprocessing as mp

class LibSeqCount:

    def __init__(self):
        self.MySeqList = list()
        self.barcode = str()
        self.readrange = (100, 2500)  # length min and max of the read
        self.results = list()
        self.alignment_ref = str()  # This is the gRNA scaffold that is used to control for false positives

    def import_seq_info(self,seqlistfile):
        f = open(seqlistfile)
        for line in f:
            #self.MySeqList.append(line[:-1].split("\t")[1])
            self.MySeqList.append(line[:-1])

    # Function for importing a fastq file and then identifying hits
    def import_fastq_data(self, file_directory, barcode, multi=False, use_alignment=False):
        os.chdir(file_directory + barcode)
        filestring = str()
        for fastq_file in os.listdir(os.curdir):
            if fastq_file.endswith(".fastq"):
                f = open(fastq_file)
                for line in f:
                    filestring += line
                f.close()
        f = open(file_directory + "alignment_librarycalloutput_20ntstrict" + barcode + ".txt", 'w')
        for sequence in self.MySeqList:
            sequence = sequence.upper()
            count = 0
            if use_alignment:
                seqcount = [m.start() for m in re.finditer(self.alignment_ref,filestring)]
                for match in seqcount:
                    grnaseq = filestring[match-20:match]
                    if pairwise2.align.globalxx(grnaseq,sequence,score_only=True) > 15:
                        count += 1
                revseqcount = [m.start() for m in re.finditer(self.revcom(self.alignment_ref),filestring)]
                for match in revseqcount:
                    grnaseq = filestring[match + len(self.alignment_ref):match + len(self.alignment_ref) + 20]
                    if pairwise2.align.globalxx(grnaseq,sequence,score_only=True) > 15:
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
L.import_seq_info("/Volumes/Lexar/CalebILseqs.txt")
#  L.import_fastq_data("/home/trinhlab/Desktop/pCasSAlib/EcoLibExtract/fastq_pass/", "", use_alignment=True)
barcodes = ["01","02","03","04","05","06"]
for num in barcodes:
    L.import_fastq_data("/Volumes/Lexar/fastq_pass/", "barcode" + num, use_alignment=True)
