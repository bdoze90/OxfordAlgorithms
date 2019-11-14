"""File takes a single or set of FASTQ files and counts the times a sequence from a list appears.  Can set a barcode
sequence as well to determine if the sequence is appearing in the right place, as well as a minimum sequence length.
Optional: Use of blast alignment algorithm to determine allow for partial matches."""

import os,re

class LibSeqCount:

    def __init__(self):
        self.MySeqList = list()
        self.barcode = str()
        self.readrange = (100,20000)  # length min and max of the read
        self.results = list()

    def import_seq_info(self):
        f = open('/Users/brianmendoza/Desktop/SaureusVipLib1.txt')
        for line in f:
            self.MySeqList.append(line[:-1])

    def import_fastq_data(self,file_directory, multi=False):
        f = open(file_directory)
        filestring = str()
        for line in f:
            filestring += line
        f.close()
        for sequence in self.MySeqList:
            seqcount = re.finditer(sequence, filestring)
            count = 0
            for match in seqcount:
                count += 1
            print(sequence, count)
            self.results.append((sequence, count))



L = LibSeqCount()
L.import_seq_info()
L.import_fastq_data("/Users/brianmendoza/Desktop/pCasSAlib/EcoLibExtract/pCasSAEcoLibExtractAll.fastq")