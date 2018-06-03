#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 11:20:35 2018
@author: sofieolundvillumsen
"""

######################################################################
############################## Import ################################
######################################################################
from Bio import SeqIO
from itertools import chain #Sorts first by the first element, then the next element 
import os

######################################################################
########################### Variables ################################
######################################################################
file = os.path.abspath("Dataset_Uniprot_AnimalToxins_All.fa")
outputfile = 'Dataset_Sorted_AnimalToxins_All.fa'

######################################################################
########################### Funktioner ###############################
######################################################################
#Function for reading and sorting data set
def ReadSeq(filename):
    Data = []
    Doublet = []
    progress = [3000,6000,9000,12000]
    count = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        count = count +1
        if count in progress: 
            print(count)
        if seq_record.id not in chain (*Data): #Doubles are sorted out in proportion to ID and sequence
            if seq_record.seq not in chain (*Data):
                Data.append((seq_record.id,seq_record.seq))
            else:
                Doublet.append((seq_record.id,seq_record.seq))
        else:
            Doublet.append((seq_record.id,seq_record.seq))
    return Data, Doublet

#Sorted data is saved into a new file, which is saved in the working directory 
def FastaOutput(SortedData,filenameOut):
    ofile = open(filenameOut, "w")
    for i in range(len(SortedData)):
        ofile.write(">" + str(SortedData[i][0]) + "\n" + str(SortedData[i][1]) + "\n")
    ofile.close()
    print('File created')


######################################################################
############################### Main #################################
######################################################################
SortedData,Doublet = ReadSeq(file)
FastaOutput(SortedData, outputfile)