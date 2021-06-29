"""File for identifying short sequences from Nanopore library sequencing data."""

from Bio import pairwise2
import subprocess

import os
import csv


def run_centrifuge_sorting():

    cent_exe = '/Users/brianmendoza/Documents/centrifuge-1.0.3-beta/centrifuge'
    ind_path = '/Users/brianmendoza/Documents/centrifuge-1.0.3-beta/indices/myexe' ###Change to appropriate index path
    master_dir = '/Volumes/Lexar/IL_RCA_6' ###Please change to the root directory of fastq files (code recursively goes through barcoded folders if present

    for barcode_dir in sorted(os.listdir(master_dir)): ###MASTER LOOP THROUGH EVERY DIRECTORY IN SEQUENCING FOLDER
        barcode_dir = master_dir + '/' + barcode_dir
        print(barcode_dir)
        master_fastq = barcode_dir + '/master.fastq'
        """ Concatenate all fastq files to a single large file """
        cat_fastq = 'cat ' + barcode_dir + '/*.fastq > ' + master_fastq
        os.system(cat_fastq)

        """  Run centrifuge """
        cline = cent_exe + ' -q -x ' + ind_path + ' ' + master_fastq +  ' > ' + os.getcwd() + '/' + barcode_dir.split('/')[-1] + '_master.tsv'
        os.system(cline)
        os.remove(master_fastq)


def concatenate_and_cross_reference(barcodes):

    gene_path = "/Volumes/Lexar/IL_RCA_Gene_6/"
    marker_path = "/Volumes/Lexar/IL_RCA_UraLeu_6/"
    fullseq_path = "/Volumes/Lexar/IL_RCA_Full_3/"

    DataFrame = dict()

    for barcode in barcodes:
        BarcodeDataFrame = dict()
        # Establish the table with marker_path
        f = open(marker_path + "barcode" + barcode + "_master.tsv")
        for line in f:
            if line.startswith("readID"):  #throw away first line
                r = 1
            else:
                readID = f.readline()[:-1].split("\t")[0]
                data = f.readline()[:-1].split("\t")[1:]
                BarcodeDataFrame[readID] = data

        f.close()

        # Compare reads with those that match the gene as well
        y = open(gene_path + "barcode" + barcode + "_master.tsv")
        for line in y:
            if line.startswith("readID"):  # throw away first line
                r = 1
            else:
                readID = y.readline()[:-1].split("\t")[0]
                data = y.readline()[:-1].split("\t")[1:]
                if readID in BarcodeDataFrame:
                    for item in data:
                        BarcodeDataFrame[readID].append(item)
        y.close()
        DataFrame[barcode] = BarcodeDataFrame

    # print out relevant data:
    with open("/Volumes/Lexar/IL_RCA_Data_6.csv", 'w') as outfile:
        outwrite = csv.writer(outfile, delimiter=',', quotechar='"')
        for barcode in DataFrame:
            count = 0
            for read in DataFrame[barcode]:
                if len(DataFrame[barcode][read]) > 7:
                    count += 1
                    outwrite.writerow([barcode,DataFrame[barcode][read][0],DataFrame[barcode][read][7]])
            print(count)

    outfile.close()




#run_centrifuge_sorting()
concatenate_and_cross_reference(["18","22","08"])