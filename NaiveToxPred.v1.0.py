#!/usr/bin/env python3                                                                                                                                                          
# -*- coding: utf-8 -*-

"""                                                                                                                                                                               
Created on Wed Mar 14 15:56:02 2018                                                                                                                                         
@author: sofieolundvillumsen                                                                                                                                                         
"""

##############################################################################################################################################                                               
################################################################### IMPORT ###################################################################                                                  
##############################################################################################################################################                                 
from Bio import SeqIO
import numpy as np
from collections import OrderedDict
import argparse
import os
import pandas as pd
import os.path
import time
import shutil
import sys


##############################################################################################################################################                                               
################################################################ VARIABLES ###################################################################                                                   
############################################################################################################################################## 
Version = 1.0
Outdir_Web = 'services/NaiveToxPred-' + str(Version) + '/tmp/NaiveToxPred.' + str(os.getpid())
Outdir_UNIX = "NaiveToxPred." + str(os.getpid())

#DataSetPos_Web = '/usr/cbs/bio/src/naivetoxpred-1.0/dataset_pos.fa'
DataSetPos_Web = '/usr/cbs/bio/src/naivetoxpred-1.0/VenomToxinsData.fa'  

#DataSetPos_UNIX = os.path.abspath("dataset_pos.fa") 
DataSetPos_UNIX = os.path.abspath("VenomToxinsData.fa") 

##############################################################################################################################################
################################################################## ARGUMENTS #################################################################
##############################################################################################################################################
# This section is for handling the arguments.
parser = argparse.ArgumentParser(description="Prediction of venom toxins based on protein sequences in fasta format.")

#Required argument                                                                                                                                                
parser.add_argument("--InputFile",
                    help="Fasta file containing query sequences. This argument is the only optional argument which is required.", required=True) 

#Optional arguments                                                                                                                                                
parser.add_argument("--Evalue",
                    help="Set desired e-value treshold. Default = '1e-5' when running BLASTp.",
                    type=float,
                    default=1e-5)
parser.add_argument("--SignalpStatus",
                    help="If input file contains both mature- and precursor protein sequences, then SignalpStatus should be False. Default = 'True'.",
                    default='True',
                    required=False)
parser.add_argument("--SignalpFile",
                    help="A General Feature Format (GFF) file from a SignalP run. The file is created by following command: signalp -n SignalpDataOut.txt <name on inputfile>.",
                    required=False)
parser.add_argument("--BLASTFile",
                    help="BLAST file obtained in the following way: 1. BLASTp database is created: makeblastdb -in dataset_pos.fa -dbtype prot and 2. Run BLASTp with created database: -de blastp -query <name on inputfile> -db dataset_pos.fa -out BLASTDataOut.txt -outfmt 6 -evalue 0.00001 -max_target_seqs 1 -num_threads 4.",
                    required=False)
parser.add_argument("--BenchmarkFile",
                    help="Benchmark file including identifier of the protein sequences, as in the input file, and their prediction assigned either T(=Toxin) or N(=NonToxin), seperated from the sequence name with a tap. As follows: <Protein sequence identifier>    T/N.",
                    required=False)
parser.add_argument("--Web",
                    help="This indicates whether the script is run from a UNIX terminal or from the webserver. When Web = True, there will be an output to the screen (stdout) and Web = False files will be created and saved in a directory. Default = 'False'.",
                    required=False,
                    default=False)
parser.add_argument("--Directory",
                    help="Directory containing the created files. Default = 'NaiveToxPred.ProcessID'.",
                    required=False)
parser.add_argument("--Verbose",
                    help="When Verbose = True, then the BLASTp and SignalP system calls are printed to the screen. Default = 'False'.",
                    default=False,
                    required=False)

args = parser.parse_args()


 
##############################################################################################################################################                                                        
################################################################## FUNCTIONS #################################################################                                                        
##############################################################################################################################################
#Load functions===============================================================================================================================
NamesPos = []
NamesInput = []

#This function loads and handles the positive data set, the InputFile, the SignalpFile and BLASTFile respectively.
def LoadSeq(PosDataSet,InputFile, SignalpDataSet, BLASTDataSet):
    for record in SeqIO.parse(PosDataSet, "fasta"):
        NamesPos.append(record.id) 
    for record in SeqIO.parse(InputFile, "fasta"):
        NamesInput.append(record.id)
    MatrixSignalpOutput = np.genfromtxt(SignalpDataSet,dtype='str') #Converts the file created with the signalp server to a convinent matrix.                               
    NamesSignalp = []
    if (np.size(MatrixSignalpOutput) == 9): 
        MatrixSignalpOutput.shape=9
        for i in range(np.size(MatrixSignalpOutput)): 
            NamesSignalp.append(MatrixSignalpOutput[0])
    else:
        NamesSignalp.append(MatrixSignalpOutput[:,0])
        NamesSignalp = [list(x) for x in NamesSignalp] #Converts the names from being stored in an array to be stored in a list.                                                           
        NamesSignalp = [item for sublist in NamesSignalp for item in sublist]
    MatrixBLASTOutput = np.genfromtxt(BLASTDataSet,dtype='str')
    NamesBLAST = []
    if (np.size(MatrixBLASTOutput) < 13):
        for i in range(np.size(MatrixBLASTOutput)):
            NamesBLAST.append(MatrixBLASTOutput[0])
    else:
        NamesBLAST.append(MatrixBLASTOutput[:,0])
        NamesBLAST = [list(x) for x in NamesBLAST] 
        NamesBLAST = [item for sublist in NamesBLAST for item in sublist]
    return NamesPos, NamesInput, MatrixSignalpOutput, NamesSignalp, MatrixBLASTOutput, NamesBLAST

#This function does the same as the LoadSeq function, except it doesn't handle a SignalpFile. This function is used when the seqeunces in the InputFile doesn't contain signal peptides.
def LoadSeqNoSigP(PosDataSet,InputFile,BLASTDataSet):
    for record in SeqIO.parse(PosDataSet, "fasta"):
        NamesPos.append(record.id)
    for record in SeqIO.parse(InputFile, "fasta"):
        NamesInput.append(record.id)  
    MatrixBLASTOutput = np.genfromtxt(BLASTDataSet,dtype='str')
    NamesBLAST = []
    if (np.size(MatrixBLASTOutput) < 13):
        for i in range(np.size(MatrixBLASTOutput)):
            NamesBLAST.append(MatrixBLASTOutput[0])
    else:
        NamesBLAST.append(MatrixBLASTOutput[:,0])
        NamesBLAST = [list(x) for x in NamesBLAST]
        NamesBLAST = [item for sublist in NamesBLAST for item in sublist] 
    return NamesPos, NamesInput, MatrixBLASTOutput, NamesBLAST



#Prediction function===========================================================================================================================                                                  
#This function sorts the MatrixBLASTOutput by evalues. E-values which is above a given treshold is assigned NonToxin and evalues under the cut-off value are assigned Toxin.
def EvalueCutOff(NamesInput, NamesBLAST, MatrixBLASTOutput): 
    NonToxin1 = list(set(NamesInput)-set(NamesBLAST)) #The names in the BLASTFile are subtract from the names in the InputFile. These names are NonToxins.
    Toxin1 = []

    if (np.size(MatrixBLASTOutput) < 13):
        if MatrixBLASTOutput[10].astype(np.float) < args.Evalue:
            Toxin1.append(MatrixBLASTOutput[0])
        else:
            NonToxin1.append(MatrixBLASTOutput[0])
    else:
        for i in range(len(MatrixBLASTOutput[:,10])): #The e-values which are found in column 11 are sorted at a given cut-off value.
            if MatrixBLASTOutput[i,10].astype(np.float) < args.Evalue: 
                Toxin1.append(MatrixBLASTOutput[i,0])
            else:           
                NonToxin1.append(MatrixBLASTOutput[i,0])       

    if args.SignalpStatus == 'True': #If the sequences have signal peptides, then they are sorted additionally with help from the signalp server.
        Toxin2=[]
        for i in range(len(Toxin1)):
            if Toxin1[i] in NamesSignalp:
                Toxin2.append(Toxin1[i])                
            else:
                NonToxin1.append(Toxin1[i])  
        Toxin2 = list(OrderedDict.fromkeys(Toxin2)) #This command filters out any dublicates there should be after the BLAST search.
        Toxin = list(set(Toxin2)-set(NonToxin1))
        NonToxin1 = list(OrderedDict.fromkeys(NonToxin1))         
        NonToxin = list(set(NonToxin1)-set(Toxin2))
    else:
        Toxin = list(OrderedDict.fromkeys(Toxin1))                                   
        NonToxin1 = list(OrderedDict.fromkeys(NonToxin1))                             
        NonToxin = list(set(NonToxin1)-set(Toxin))
    return Toxin, NonToxin



#Output functions===============================================================================================================================                                                    
BlastYN = []
ToxinYN = []

#This function creates a organized output for the user, this function creates output with information from both BLASTp and signalp.
def UserOutput(NamesInput, NamesSignalp, MatrixSignalpOutput, NamesBLAST, MatrixBLASTOutput, Toxin):
    #The lists are sorted, to obtain same order in the columns.
    NamesInput.sort()
    NamesSignalp.sort()
    np.sort(MatrixSignalpOutput)
    Toxin.sort()
    SignalpYN = []
    for i in range(len(NamesInput)):
        if NamesInput[i] in NamesSignalp:
            SignalpYN.append('Yes')
        else:
            SignalpYN.append('No')            
    SignalpStartN = []
    for i in range(len(SignalpYN)):
            if SignalpYN[i]=='No':
                SignalpStartN.append('-')
            else:
                SignalpStartN.append('Yes')

    SignalpEndN = []
    for i in range(len(SignalpYN)):
        if SignalpYN[i]=='No':
            SignalpEndN.append('-')
        else:
            SignalpEndN.append('Yes')
    
    SignalpScoreN = []
    for i in range(len(SignalpYN)):
        if SignalpYN[i]=='No':
            SignalpScoreN.append('-')
        else:
            SignalpScoreN.append('Yes')
                

    if (np.size(MatrixSignalpOutput) == 9):
        SignalpStartN = MatrixSignalpOutput[3]
        SignalpEndN = MatrixSignalpOutput[4]
        SignalpScoreN = MatrixSignalpOutput[5]
        
    else:
        k = 0
        MatrixSignalpStart = MatrixSignalpOutput[:,3]
        for i in range(len(SignalpStartN)):
            #In the list SignalpStartN, the positions assigned "Yes" are being replaced by the start position of the signal peptide.                                       
            if SignalpStartN[i]=='Yes':
                for j in range(len(MatrixSignalpStart)):
                    if k < len(MatrixSignalpStart) and i < len(SignalpStartN): #This is necessary, otherwise the loop is going to run too many times.           
                        SignalpStartN[i] = MatrixSignalpStart[k]
                        i = i + 1
                        k = k + 1
                        break
        l = 0
        MatrixSignalpEnd = MatrixSignalpOutput[:,4]
        for i in range(len(SignalpEndN)):
            if SignalpEndN[i]=='Yes':
                for j in range(len(MatrixSignalpEnd)):
                    if l < len(MatrixSignalpEnd) and i < len(SignalpEndN):
                        SignalpEndN[i] = MatrixSignalpEnd[l]
                        i = i + 1
                        l = l + 1
                        break
        m = 0
        MatrixSignalpScore = MatrixSignalpOutput[:,5]
        for i in range(len(SignalpScoreN)):
            if SignalpScoreN[i]=='Yes':
                for j in range(len(MatrixSignalpScore)):
                    if m < len(MatrixSignalpScore) and i < len(SignalpScoreN):
                        SignalpScoreN[i] = MatrixSignalpScore[m]
                        i = i + 1
                        m = m + 1
                        break

    for i in range(len(NamesInput)):
        if NamesInput[i] in NamesBLAST:
            BlastYN.append('Yes')
        else:
            BlastYN.append('No')        
            
    EvalueN = []
    for i in range(len(BlastYN)):
            if BlastYN[i]=='No':
                EvalueN.append('-')
            else:
                EvalueN.append('Yes')       

    #The following will remove all the dublicate rows from the matrix, leaving the best BLAST hit in MatrixBLASTOutputSorted. This is the work of the unique function. 
    if (np.size(MatrixBLASTOutput) < 13):
        MatrixBLASTOutputSorted = MatrixBLASTOutput
        np.sort(MatrixBLASTOutputSorted)
        EvalueN = MatrixBLASTOutputSorted[10]

    else:
        _, indices = np.unique(MatrixBLASTOutput[:, 0], return_index=True) #The uniqe function returns the indices of the first occurences of the unique values, aka. the best BLAST score.
        MatrixBLASTOutputSorted = MatrixBLASTOutput[indices, :]
        np.sort(MatrixBLASTOutputSorted)
        MatrixEvalue = MatrixBLASTOutputSorted[:,10]

        c = 0
        for i in range(len(EvalueN)):
            if EvalueN[i]=='Yes':
                for j in range(len(MatrixEvalue)):
                    if c < len(MatrixEvalue) and i < len(EvalueN):
                        EvalueN[i] = MatrixEvalue[c]
                        i = i + 1
                        c = c + 1
                        break

    for i in range(len(NamesInput)):
        if NamesInput[i] in Toxin:
            ToxinYN.append('Yes')
        else:
            ToxinYN.append('No')            
    data = {'Identifier':  NamesInput,'SPstart': SignalpStartN,'SPend': SignalpEndN,'SPscore': SignalpScoreN,'SP': SignalpYN,'BLAST': BlastYN,'Evalue': EvalueN,'Toxin': ToxinYN}
    MatrixFinalOutput = pd.DataFrame(data, index=[list(range(0,len(NamesInput)))], columns=['Identifier','SPstart','SPend','SPscore','SP','BLAST','Evalue','Toxin'])
    Prediction = "\nPredicted {:d} toxins out of {:d} input sequences.\n".format(len(Toxin),len(NamesInput))
    print(Prediction)
    MatrixFinalOutput.to_csv(OutdirPrediction, header=True, index = False, sep='\t') #The output is saved as a .txt file.                                                                     
    return MatrixFinalOutput


#This function does the same as the UserOutput function, the difference is that it will not include a signalp search.
def UserOutputNoSigP(NamesInput ,MatrixBLASTOutput, Toxin):
    NamesInput.sort()
    Toxin.sort()
    for i in range(len(NamesInput)):
        if NamesInput[i] in NamesBLAST:
            BlastYN.append('Yes')
        else:
            BlastYN.append('No')        

    EvalueN = []
    for i in range(len(BlastYN)):
            if BlastYN[i]=='No':
                EvalueN.append('-')
            else:
                EvalueN.append('Yes')    

    
    if (np.size(MatrixBLASTOutput) < 36):
        MatrixBLASTOutputSorted = MatrixBLASTOutput
        np.sort(MatrixBLASTOutputSorted)
        EvalueN = MatrixBLASTOutputSorted[10]
        
    else:
        _, indices = np.unique(MatrixBLASTOutput[:, 0], return_index=True)
        MatrixBLASTOutputSorted = MatrixBLASTOutput[indices, :]
        np.sort(MatrixBLASTOutputSorted)
        MatrixEvalue = MatrixBLASTOutputSorted[:,10]

        c = 0
        for i in range(len(EvalueN)):
            if EvalueN[i]=='Yes':
                for j in range(len(MatrixEvalue)):
                    if c < len(MatrixEvalue) and i < len(EvalueN):
                        EvalueN[i] = MatrixEvalue[c]
                        i = i + 1
                        c = c + 1
                        break   

    for i in range(len(NamesInput)):
        if NamesInput[i] in Toxin:
            ToxinYN.append('Yes')
        else:
            ToxinYN.append('No')                
    data = {'Identifier':  NamesInput,'BLAST': BlastYN,'Evalue': EvalueN,'Toxin': ToxinYN}
    MatrixFinalOutput = pd.DataFrame(data, index=[list(range(0,len(NamesInput)))], columns=['Identifier','BLAST','Evalue','Toxin'])
    print ("\nPredicted {:d} toxins out of {:d} input sequences.\n".format(len(Toxin),len(NamesInput)))
    MatrixFinalOutput.to_csv(OutdirPrediction, header=True, index = False, sep='\t') #The output is saved as a .txt file.
    
    return MatrixFinalOutput #NoSignalpOutput                                                                                                                           


#Performance function===============================================================================================================================                                                    
#This function is used if the user have a file with some premade assumptions according to Toxin and NonToxin classification.
def Performance(BenchmarkDataSet,Toxin, NonToxin, NamesPos):
    MatrixBenchmarkDataSet = np.genfromtxt(BenchmarkDataSet,dtype='str') 
    FalsePos = list(set(Toxin)-set(NamesPos))
    TruePos = list(set(Toxin)-set(FalsePos))
    TrueNeg = list(set(NonToxin)-set(NamesPos))
    FalseNeg = list(set(NonToxin)-set(TrueNeg))    
    #If the dataset only contains either toxins or non-toxins and they are classified right. 
    if len(TruePos)==0 and len(FalseNeg)==0 or len(TrueNeg)==0 and len(FalsePos)==0 or len(TrueNeg)==0 and len(FalseNeg)==0 or len(TruePos)==0 and len(FalsePos)==0 or len(TruePos)==0 and len(FalsePos)==0 and len(FalseNeg)==0:
        Accuracy = (len(TrueNeg) + len(TruePos))/(len(TrueNeg)+len(TruePos)+len(FalsePos)+len(FalseNeg))
        MatrixPerformanceOutput = "##By evalue {0:.2e}: TP = {1:d}, FP = {2:d}, TN = {3:d}, FN = {4:d}, Accuracy = {5:.4f}".format(args.Evalue,  len(TruePos),len(FalsePos),len(TrueNeg),len(FalseNeg), Accuracy)
    else: 
        Accuracy = (len(TrueNeg) + len(TruePos))/(len(TrueNeg)+len(TruePos)+len(FalsePos)+len(FalseNeg))
        Specificity = len(TrueNeg)/(len(TrueNeg+FalsePos)) #SPEC
        Sensitivity = len(TruePos)/(len(TruePos+FalseNeg)) #SENS
        BalancedAccuracy = (Specificity+Sensitivity)/2 #B.ACC
        NegativePredictedValue = len(TrueNeg)/(len(TrueNeg+FalseNeg)) #NPV
        PositivePredictedValue = len(TruePos)/(len(TruePos+FalsePos)) #PPV
        FScore = 2*len(TruePos)/(2*len(TruePos+FalsePos+FalseNeg)) #F1
        MatthewsCorCoef = (len(TruePos)*len(TrueNeg)-len(FalsePos)*len(FalseNeg))/(len(TruePos+FalsePos)*len(TruePos+FalseNeg)*len(TrueNeg+FalsePos)*len(TrueNeg+FalseNeg))**(.5) #MCC    
        MatrixPerformanceOutput = "##When the evalue is {0:.2e}: TP = {1:d}, FP = {2:d}, TN = {3:d}, FN = {4:d}, Accuracy = {5:.2f}, PPV = {6:.2f}, NPV = {7:.2f}, Sensitivity = {8:.2f}, Specificity = {9:.2f}, F1 = {10:.2f}, MCC = {11:.4f}".format(args.Evalue, len(TruePos),len(FalsePos),len(TrueNeg),len(FalseNeg),Accuracy,PositivePredictedValue,NegativePredictedValue, Sensitivity,Specificity, FScore, MatthewsCorCoef)        
    file = open(OutdirPerformance,'w')
    file.write(MatrixPerformanceOutput)
    file.close()
    return MatrixPerformanceOutput


#BLAST and SignalP functions=====================================================================================================================                                         
#This function runs BLASTp on the InputFile.
def RunBLAST(InputFile):
    if args.Web is False:
        MakeDataBase = "makeblastdb -in " + DataSetPos +  " -dbtype prot > /dev/null" #/dev/null hides irrelevant screen output.
        os.system(MakeDataBase)                                                                                          
        BLASTCmd = "blastp -query " + args.InputFile + " -db " + DataSetPos +  " -out " + OutdirBLAST + " -outfmt 6 -evalue 0.00001 -max_target_seqs 1 -num_threads 4"
    else:
        BLASTCmd = "blastall -p blastp -i " + args.InputFile + " -d " + DataSetPos + " -o " + OutdirBLAST + " -m 8 -e 0.00001 -b 1" #-b 1, takes only the best hit
    os.system(BLASTCmd)
    #The BLAST file is made and the wait command is used when the file is being created                                                           
    while not os.path.exists(OutdirBLAST):
        time.sleep(30)
    if os.path.isfile(OutdirBLAST):
        BLASTCmdOut = os.path.abspath(OutdirBLAST)
    else:
        raise ValueError("%s isn't a file!" % os.path.abspath(OutdirBLAST))    
    return BLASTCmdOut


#This functions runs Signalp on the InputFile
def RunSignalp(InputFile):
    SignalpCmd = "signalp -n " + OutdirSigP + " " + args.InputFile + " > /dev/null"
    os.system(SignalpCmd)
    
    #while not os.path.exists(OutdirSigP):
     #   time.sleep(10)
    if os.path.isfile(OutdirSigP):
        SignalpCmdOut = os.path.abspath(OutdirSigP)
    else:
        print('test')
        SignalpCmdOut=None
        args.SignalpStatus='False'
    return SignalpCmdOut



##############################################################################################################################################
#################################################################### MAIN ####################################################################
##############################################################################################################################################
#Check file==================================================================================================================================                                                  
if not os.path.isfile(args.InputFile):
    print("Input file not found")
    sys.exit(0)





#Web?=========================================================================================================================================                                                   
if args.Web is not False:
    Outdir = Outdir_Web
    DataSetPos = DataSetPos_Web
else:
    if args.Directory is None:
        Outdir = Outdir_UNIX
    else:
        Outdir = args.Directory
    DataSetPos = DataSetPos_UNIX

OutdirBLAST = Outdir + '/BLASTDataOut.txt'
OutdirPrediction = Outdir  + '/Prediction.txt'

MakedirCmd = 'mkdir -p ' + Outdir
os.system(MakedirCmd)


#Calling functions============================================================================================================================                                               
if args.BLASTFile is None:
    BLASTCmd = RunBLAST(args.InputFile)
else:
    BLASTCmd = args.BLASTFile

if os.stat(BLASTCmd).st_size == 0:
    print("No BLAST hits, thus no predicted toxins.")
    sys.exit(0)


if args.SignalpFile is None:
    OutdirSigP = Outdir + '/SignalpDataOut.txt'
    SignalpCmd = RunSignalp(args.InputFile)
else:
    SignalpCmd = args.SignalpFile


if args.SignalpStatus == 'True':
    NamesPos,NamesInput,MatrixSignalpOutput,NamesSignalp,MatrixBLASTOutput, NamesBLAST = LoadSeq(DataSetPos,args.InputFile,SignalpCmd,BLASTCmd)
    Toxin, NonToxin = EvalueCutOff(NamesInput,NamesBLAST,MatrixBLASTOutput)
    MatrixFinalOutput = UserOutput(NamesInput,NamesSignalp,MatrixSignalpOutput,NamesBLAST,MatrixBLASTOutput,Toxin)
else:
    NamesPos,NamesInput,MatrixBLASTOutput,NamesBLAST = LoadSeqNoSigP(DataSetPos,args.InputFile,BLASTCmd)
    Toxin,NonToxin = EvalueCutOff(NamesInput,NamesBLAST,MatrixBLASTOutput)
    MatrixFinalOutput = UserOutputNoSigP(NamesInput,MatrixBLASTOutput,Toxin)



if args.BenchmarkFile is not None:
    OutdirPerformance = Outdir  + '/Performance.txt'
    MatrixPerformanceOutput = Performance(args.BenchmarkFile,Toxin,NonToxin,NamesPos)

if args.Web is False:
    if args.BenchmarkFile is not None: #For when PerformanceLoop.sh is run
        print(MatrixPerformanceOutput)
else:
    print(MatrixFinalOutput)


#Verbose=======================================================================================================================================    
if args.Verbose is not False:
    print ("\nMaking directory: " + (MakedirCmd))
    print ("Output BLAST file: " +  OutdirBLAST)
    print ("Output prediction file: " + OutdirPrediction)
    if args.SignalpStatus == 'True':
        print ("Output SignalP file: " + OutdirSigP)
    if args.BenchmarkFile is not None:
        print ("Output performance file: " + OutdirPerformance)
    print ("SignalpStatus = " + str(args.SignalpStatus))

    BLASTCmd_Web = "blastall -p blastp -i " + args.InputFile + " -d " + DataSetPos + " -o " + OutdirBLAST + " -m 8 -e 0.00001 -b 1"
    BLASTCmd_UNIX = "blastp -query " + args.InputFile + " -db " + DataSetPos+  " -out " + OutdirBLAST + " -outfmt 6 -evalue 0.00001 -max_target_seqs 1 -num_threads 4"
    if args.Web is not False:
        print("\nBLASTp system call for web mode:\n" + BLASTCmd_Web)
    else:
        print("\nBLASTp system call for UNIX mode:\n" + BLASTCmd_UNIX)
    if args.SignalpStatus == 'True':
        SignalpCmd = "signalp -n " + OutdirSigP + " " + args.InputFile
        print("\nSignalp system call:\n" + SignalpCmd)
