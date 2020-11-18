#!/usr/bin/python3
#Import Packages:
import argparse
import pysam
import FLAME_FUNC.FLAME_FUNC as FF #Break this down into each part based on the command.

#Input command:
parser = argparse.ArgumentParser(description = "FLAME: Full Length Adjacency Matrix Enumeration")
##Obligatory Inputs:
parser.add_argument("-I", dest = "INPUT", help = "Input file")
parser.add_argument("-GTF", dest = "GTF", help = "Reference File in GTF format")
parser.add_argument("-G", dest = "GENE", help = "Targeted Gene")
parser.add_argument("--range", dest = "RANGE", help = "Variance Range", default = 20)
parser.add_argument("-O", dest = "OUTPUT", help = "Output Prefix", default = "Flame")
##Optional Inputs:
parser.add_argument("-R", dest = "REF", help = "Reference File in Fasta format", default = 0)
#parser.add_argument("-THRESHOLD?", dest = "THRESHOLD", help = "Threshold for the Frequency Analysis, the lower threshold the more comprehensive but the longer processing power required", default = 0.01)?
parser.add_argument("-B", dest = "SAM", help = "Shortread Sequence", default = 0)
parser.add_argument("--verbose", dest = "VERBOSE", help = "Verbose output", action="store_true")
#ADD A FUNCTION THAT "COLLAPSES" THE QUANTIFICATION READS?
##Parser Funciton:
args = parser.parse_args()
print("\n-------------------------------------------------------------------------------------------")
print("\n-----------\tFLAME: Full Length Adjacency Matrix Enumeration\t\t\t-----------\n")
print("-------------------------------------------------------------------------------------------")
print("\n-----------\tInitiating FLAME\t\t\t\t\t\t-----------")
print("Input:\t\t{}\n\
GTF:\t\t{}\n\
Gene:\t\t{}\n\
Range:\t\t{}\n\
Output:\t\t{}-[Suffix]".format(
    args.INPUT,
    args.GTF,
    args.GENE,
    args.RANGE,
    args.OUTPUT
))
##The optional choices:
###The Optional Removal of both.
if args.REF == 0 and args.SAM == 0:
    print("\n")
###Load the Reference if it is specified to be used for the Adjacent Splice Site Signal (S3) Detection.
if args.REF != 0:
    print("Reference:\t{}\n".format(args.REF), end = "")
    Flame_Main_REFFILE = open(args.REF, "r")
    #Remove the header of the reference as well as make it one continuous string.
    Flame_Main_REF = ""
    for Flame_Main_REFLINE in Flame_Main_REFFILE:
        if ">" in Flame_Main_REFLINE:
            pass
        else:
            Flame_Main_REF += Flame_Main_REFLINE.rstrip()
###Load the Shortread if it specified to be used for the confirmation of splice signal using short read RNA-seq.
if args.SAM != 0:
    print("Shortread:\t{}\n".format(args.SAM), end = "")
    Flame_Main_SHORTREAD1 = str(args.SAM)
    if Flame_Main_SHORTREAD1.endswith(".sam") or Flame_Main_SHORTREAD1.endswith(".SAM"):
        Flame_Main_SHORTREAD2 = pysam.AlignmentFile("%s" %args.SAM, "r")
        Flame_Main_SHORTREAD3 = pysam.AlignmentFile("%s" %args.SAM, "r")
    elif Flame_Main_SHORTREAD1.endswith(".bam") or Flame_Main_SHORTREAD1.endswith(".BAM"):
        Flame_Main_SHORTREAD2 = pysam.AlignmentFile("%s" %args.SAM, "rb")
        Flame_Main_SHORTREAD3 = pysam.AlignmentFile("%s" %args.SAM, "rb")
    else:
        print("Shortread input file format not recognized ([File].bam/[File].BAM/[File].sam/[File].SAM)")
###Setting the threshold for frequency required:
#if args.THRESHOLD:
#    print("Threshold:\t{}\n".format(args.THRESHOLD), end = "")
###Output a Verbose output that includes the raw "Faulty Ranges" and the Reference Translate Table.
if args.VERBOSE:
    print("Verbose:\tOn\n")
print("-------------------------------------------------------------------------------------------\n")

#Variables:
##Input Files and References:
Flame_Main_INPUT1 = open(args.INPUT, "r")
Flame_Main_INPUT2 = open(args.INPUT, "r")
Flame_Main_GTFFILE = open(args.GTF, "r") #Storage[List] of the GTF file.
Flame_Main_RANGESIZE = int(args.RANGE) #Storage(Integer) of the Windowsize.
Flame_Main_GENENAME = args.GENE #Storage("String") of the Specific Gene.
Flame_Main_REFERENCELIST = [] #Storage[Nested List] of the reference.
Flame_Main_FREQUENCYWINDOWSIZE = 2 #FIX SO THAT THIS WINDOWSIZE IS FLEXIBLE
##FLAME: Printing Files:
Flame_Main_Counter1 = 0
Flame_Main_Counter2 = 0
Flame_Main_Counter3_1 = False
Flame_Main_Counter3_2 = False
Flame_Main_Counter4 = 0
Flame_Main_Counter5 = 0

print("-----------\tInitiate Creation of Reference\t\t\t\t\t-----------")
#Function to transform the GTF-file into an efficient python nested list: Exon[Name, Start, Length, Stop].
#Note that this reference creator is not aware of metadata such as which chromosome (In multichromosomal organism).
Flame_Main_REFERENCELIST = FF.CREATEREFFUNC(Flame_Main_GTFFILE,
                                            args.GENE) #((The GTF Input file), (The name of the Gene))

print("-----------\tInitiate Filter Function\t\t\t\t\t-----------")
#Function to classify the reads as either "Correct" or "Faulty" and store them as 
Flame_Main_CORRECREADS, Flame_Main_FAULTYREADS = FF.FILTERFUNC(Flame_Main_INPUT1,
                                                               Flame_Main_INPUT2,
                                                               Flame_Main_REFERENCELIST,
                                                               Flame_Main_RANGESIZE) #Input, Reference, Rangesize, CorrectOutput, 

##Print out Correct:
if Flame_Main_CORRECREADS != []:
    Flame_Main_OUTPUT = open("%s.Correct.bed" %args.OUTPUT, "w+")
    for Flame_Main_COUNT in Flame_Main_CORRECREADS:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
    Flame_Main_OUTPUT.close()
else:
    pass

##Print out Faulty:
if Flame_Main_FAULTYREADS != []:
    Flame_Main_OUTPUT = open("%s.Faulty.bed" %args.OUTPUT, "w+")
    for Flame_Main_COUNT in Flame_Main_FAULTYREADS:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
    Flame_Main_OUTPUT.close()
else:
    pass   

#Translate Function. They have built in Print function. Also make a function that will make either Correct/Faulty reads = 0 to avoid downstream computation, if one chooses.
Flame_Main_TRANSLATECORREC, Flame_Main_TRANSLATEFAULTY = FF.TRANSLATEFUNC(Flame_Main_REFERENCELIST,
                                                                          Flame_Main_RANGESIZE,
                                                                          Flame_Main_CORRECREADS,
                                                                          Flame_Main_FAULTYREADS)
if args.VERBOSE:
    ##Print out Correct:
    if Flame_Main_TRANSLATECORREC != []:
        Flame_Main_OUTPUT = open("%s.CorrectTrans.txt" %args.OUTPUT, "w+")
        for Flame_Main_COUNT in Flame_Main_TRANSLATECORREC:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
        Flame_Main_OUTPUT.close()
    else:
        pass
    
if args.VERBOSE:
    ##Print out Faulty:
    if Flame_Main_TRANSLATEFAULTY != []:
        Flame_Main_OUTPUT = open("%s.FaultyTrans.txt" %args.OUTPUT, "w+")
        for Flame_Main_COUNT in Flame_Main_TRANSLATEFAULTY:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
        Flame_Main_OUTPUT.close()
    else:
        pass

#Correct Adjacency Matrix Function:
print("-----------\tInitiate Creation of Empty Adjacency Matrix\t\t\t-----------")
Flame_Main_ADJMTX1 = FF.EMPTYADJMTXFUNC(Flame_Main_REFERENCELIST)
print("-----------\tInitiate Creation of Adjacency Matrix, Correct\t\t\t-----------")
Flame_Main_ADJMTX1 = FF.CORRECTADJMTXFUNC(Flame_Main_TRANSLATECORREC,
                                          Flame_Main_REFERENCELIST,
                                          Flame_Main_ADJMTX1)
##Print out Adjacency Matrix, Correct:
###Print the Column- and Rowheaders.
Flame_Main_OUTPUT = open("%s.AdjacencyCorrect.tsv" %args.OUTPUT, "w+")
for Flame_Main_COUNT in Flame_Main_ADJMTX1:
    if Flame_Main_Counter1 <= (len(Flame_Main_REFERENCELIST)-1):
        Flame_Main_OUTPUT.write("\t" + str(Flame_Main_REFERENCELIST[Flame_Main_Counter1][0]))
        Flame_Main_Counter1 += 1
    elif Flame_Main_Counter1 > (len(Flame_Main_REFERENCELIST)-1):
        Flame_Main_OUTPUT.write("\t" + "END" + "\t")
Flame_Main_OUTPUT.write("\n")
###Filling the Adjacency Matrix itself.
for Flame_Main_COUNT1 in Flame_Main_ADJMTX1:
    if Flame_Main_Counter2 <= (len(Flame_Main_REFERENCELIST)-1):
        Flame_Main_OUTPUT.write(str(Flame_Main_REFERENCELIST[Flame_Main_Counter2][0]) + "\t")
        for Flame_Main_COUNT2 in Flame_Main_COUNT1:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT2) + "\t")
        Flame_Main_OUTPUT.write("\n")
        Flame_Main_Counter2 += 1
    elif Flame_Main_Counter2 > (len(Flame_Main_REFERENCELIST)-1):
        pass
Flame_Main_OUTPUT.close()


print("-----------\tInitiate Singling of Faulty Exons\t\t\t\t-----------")
Flame_Main_POTENTIALS = FF.FAULTYSEPERATORFUNC(Flame_Main_TRANSLATEFAULTY,
                                               Flame_Main_REFERENCELIST)

print("-----------\tInitiate Novel Splice Site Detection, Part1: Frequency\t\t-----------")
Flame_Main_GENEREFERENCE = FF.FREQUENCYSITEFUNC(Flame_Main_POTENTIALS,
                                                Flame_Main_REFERENCELIST,
                                                Flame_Main_RANGESIZE,
                                                Flame_Main_FREQUENCYWINDOWSIZE)


print("-----------\tInitiate Novel Splice Site Detection, Part2: Threshold\t\t-----------")
#Prepare GTF reference:
Flame_Main_SPLICECANDIDATES = []
#Command itself:
Flame_Main_SPLICECANDIDATES = FF.FREQUENCYTHRESHFUNC(Flame_Main_GENEREFERENCE,
                                                     0.01, #args.THREHSOLD, #MAKE CHANGEABLE (?)
                                                     Flame_Main_FAULTYREADS)

if args.REF != 0:
    print("-----------\tInitiate Novel Splice Site Detection, Part3: Splice Signal\t-----------")
    Flame_Main_SPLICECANDIDATES = FF.SPLICESIGNALFUNC(Flame_Main_SPLICECANDIDATES,
                                                      Flame_Main_REF)
    Flame_Main_Counter3_1 = True
else:
    pass


if args.SAM != 0:
    print("-----------\tInitiate Novel Splice Site Detection, Part4: Shortread\t\t-----------")
    #Prepare a tmp dictionary:
    Flame_Main_SPLICESITECOUNT = {}
    #Command itself:
    Flame_Main_SPLICECANDIDATES, Flame_Main_SPLICESITECOUNT = FF.SHORTREADFUNC(Flame_Main_SHORTREAD2,
                                                                               Flame_Main_SHORTREAD3,
                                                                               Flame_Main_SPLICECANDIDATES,
                                                                               Flame_Main_REFERENCELIST)
    #The Raw Quantification of Shortread Splice Sites:
    Flame_Main_OUTPUT = open("%s.ShortreadSplice.tsv" %args.OUTPUT, "w+") #Make this an verbose option.
    Flame_Main_OUTPUT.write("Genomic Position" +
                            "\t" +
                            "Number" +
                            "\n")
    for k, v in Flame_Main_SPLICESITECOUNT.items():
        Flame_Main_OUTPUT.write(str(k) +
                                "\t" +
                                str(v) +
                                "\n")
    Flame_Main_OUTPUT.close()
    Flame_Main_Counter3_2 = True
else:
    pass


##Print out Potential Splice Sites:
Flame_Main_OUTPUT = open("%s.PotentialSplice.tsv" %args.OUTPUT, "w+")
if Flame_Main_Counter3_1 == False and Flame_Main_Counter3_2 == False:
    Flame_Main_OUTPUT.write("Gene Position" +
                            "\t" +
                            "Supporting Faulty Reads, Absolute" +
                            "\t" +
                            "Supporting Faulty Reads, Percent" +
                            "\n")
    for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) +
                                "\t" +
                                str(Flame_Main_COUNT[1]) +
                                "\t" +
                                str(Flame_Main_COUNT[2]) +
                                "\n")
    Flame_Main_OUTPUT.close
elif Flame_Main_Counter3_1 == True and Flame_Main_Counter3_2 == False :
    Flame_Main_OUTPUT.write("Gene Position" +
                            "\t" +
                            "Supporting Faulty Reads, Absolute" +
                            "\t" +
                            "Supporting Faulty Reads, Percent" +
                            "\t" +
                            "Adjacent Splice Signal" +
                            "\n")
    for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) +
                                "\t" +
                                str(Flame_Main_COUNT[1]) +
                                "\t" +
                                str(Flame_Main_COUNT[2]) +
                                "\t" +
                                str(Flame_Main_COUNT[3]) +
                                "\n")
    Flame_Main_OUTPUT.close
elif Flame_Main_Counter3_1 == False and Flame_Main_Counter3_2 == True:
    Flame_Main_OUTPUT.write("Gene Position" +
                            "\t" +
                            "Supporting Faulty Reads, Absolute" +
                            "\t" +
                            "Supporting Faulty Reads, Percent" +
                            "\t" +
                            "Short Read Support" +
                            "\n")
    for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) +
                                "\t" +
                                str(Flame_Main_COUNT[1]) +
                                "\t" +
                                str(Flame_Main_COUNT[2]) +
                                "\t" +
                                str(Flame_Main_COUNT[3]) +
                                "\n")
    Flame_Main_OUTPUT.close
elif Flame_Main_Counter3_1 == True and Flame_Main_Counter3_2 == True:
    Flame_Main_OUTPUT.write("Gene Position" +
                            "\t" +
                            "Supporting Faulty Reads, Absolute" +
                            "\t" +
                            "Supporting Faulty Reads, Percent" +
                            "\t" +
                            "Adjacent Splice Signal" +
                            "\t" +
                            "Short Read Support" +
                            "\n")
    for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) +
                                "\t" +
                                str(Flame_Main_COUNT[1]) +
                                "\t" +
                                str(Flame_Main_COUNT[2]) +
                                "\t" +
                                str(Flame_Main_COUNT[3]) +
                                "\t" +
                                str(Flame_Main_COUNT[4]) +
                                "\n")
    Flame_Main_OUTPUT.close
else:
    print("Error, Type: Counter3")

print("-----------\tInitiate Quantification of each Splice Permutaiton, Correct\t-----------")
Flame_Main_QUANTIFYCORRECT = FF.QUANTIFYFUNC(Flame_Main_TRANSLATECORREC)
Flame_Main_READTOTALLENGTH = []
##Print out quantification:
if Flame_Main_QUANTIFYCORRECT != {}:
    Flame_Main_OUTPUT = open("%s.QuantCorrect.tsv" %args.OUTPUT, "w+")
    for k, v in Flame_Main_QUANTIFYCORRECT.items():
        #-----------The Function for Summing the Total Length of the Splice Variant-----------#
        Flame_Main_READTOTALLENGTH = []
        for Flame_Main_CORRECTEXON in k.split("-"):
            for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
                if Flame_Main_CORRECTEXON == Flame_Main_COUNT[0]:
                    Flame_Main_READTOTALLENGTH.append(int(Flame_Main_COUNT[2]))
                elif Flame_Main_CORRECTEXON != Flame_Main_COUNT[0]:
                    pass
        #-----------The Function for Summing the Total Length of the Splice Variant-----------#
        Flame_Main_OUTPUT.write(str(v) +
                                "\t" +
                                str(sum(Flame_Main_READTOTALLENGTH)) +
                                "\t" +
                                str(len(Flame_Main_READTOTALLENGTH)) +
                                "\t" +
                                str(k) +
                                "\n")
    Flame_Main_OUTPUT.close()
else:
    pass

print("-----------\tInitiate Quantification of each Splice Permutaiton, Faulty\t-----------")
Flame_Main_QUANTIFYFAULTY = FF.QUANTIFYFUNC(Flame_Main_TRANSLATEFAULTY)
##Print out FAULTYquantification:
if Flame_Main_QUANTIFYFAULTY != {}:
    Flame_Main_OUTPUT = open("%s.QuantFaulty.tsv" %args.OUTPUT, "w+")
    for k, v in Flame_Main_QUANTIFYFAULTY.items():        
        Flame_Main_OUTPUT.write(str(v) +
                                "\t" +
                                str(k) +
                                "\n")
    Flame_Main_OUTPUT.close()
else:
    pass


#Create an option choice of one wants to produce a Adjacency Matrix:
print("-----------\tInitiate Creation of Empty Adjacency Matrix\t\t\t-----------")
Flame_Main_ADJMTX2 = FF.EMPTYADJMTXFUNC(Flame_Main_SPLICECANDIDATES)
print("-----------\tInitiate Creation of Adjacency Matrix, Faulty\t\t\t-----------")
Flame_Main_ADJMTX2 = FF.FAULTYADJMTXFUNC(Flame_Main_POTENTIALS,
                                         Flame_Main_SPLICECANDIDATES,
                                         Flame_Main_ADJMTX2,
                                         Flame_Main_RANGESIZE)
##Print out Adjacency Matrix, Faulty:
###Print the Column- and Rowheaders.
Flame_Main_OUTPUT = open("%s.AdjacencyFaulty.tsv" %args.OUTPUT, "w+")
for Flame_Main_COUNT in Flame_Main_ADJMTX2:
    if Flame_Main_Counter4 <= (len(Flame_Main_SPLICECANDIDATES)-1):
        Flame_Main_OUTPUT.write("\t" + str(Flame_Main_SPLICECANDIDATES[Flame_Main_Counter4][0]))
        Flame_Main_Counter4 += 1
    elif Flame_Main_Counter4 > (len(Flame_Main_SPLICECANDIDATES)-1):
        Flame_Main_OUTPUT.write("\t" + "Below 1%")
        Flame_Main_Counter4 += 1
Flame_Main_OUTPUT.write("\n")
###Filling the Adjacency Matrix itself.
for Flame_Main_COUNT1 in Flame_Main_ADJMTX2:
    if Flame_Main_Counter5 <= (len(Flame_Main_SPLICECANDIDATES)-1):
        Flame_Main_OUTPUT.write(str(Flame_Main_SPLICECANDIDATES[Flame_Main_Counter5][0]) + "\t")
        for Flame_Main_COUNT2 in Flame_Main_COUNT1:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT2) + "\t")
        Flame_Main_Counter5 += 1
        Flame_Main_OUTPUT.write("\n")
    elif Flame_Main_Counter5 == (len(Flame_Main_SPLICECANDIDATES)):
        Flame_Main_OUTPUT.write("Below 1%" + "\t")
        for Flame_Main_COUNT2 in Flame_Main_COUNT1:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT2) + "\t")
        Flame_Main_OUTPUT.write("\n")                        
Flame_Main_OUTPUT.close()

##The Verbose Output Option:
###Print out the raw potential ranges:
if args.VERBOSE:
    Flame_Main_OUTPUT = open("%s.RawRanges.txt" %args.OUTPUT, "w+")
    for Flame_Main_COUNT in Flame_Main_POTENTIALS:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
    Flame_Main_OUTPUT.close()


### Print out the Reference
if args.VERBOSE:
    Flame_Main_OUTPUT = open("%s.Reference.txt" %args.OUTPUT, "w+")
    Flame_Main_OUTPUT.write("Exon Name" +
                            "\t" +
                            "Exon Start Site" +
                            "\t" +
                            "Exon Length" +
                            "\t" +
                            "Exon Stop Site" +
                            "\n")
    for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) +
                                "\t" +
                                str(Flame_Main_COUNT[1]) +
                                "\t" +
                                str(Flame_Main_COUNT[2]) +
                                "\t" +
                                str(Flame_Main_COUNT[3]) +
                                "\n")
    Flame_Main_OUTPUT.close()
