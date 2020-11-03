#!/usr/bin/python3
#Import Packages:
import re #time, sys, os, re
import argparse
import pysam

def CREATEREF(CRINPUT, CRLIST):
    for EXON1 in CRINPUT.read().split("\n"):
        if "\texon\t" in EXON1 and len(EXON1) >= 1:
            REFSPLICESTART1 = (int(EXON1.split("\t")[3]) - STARTPOSITION) #HAVE A FUNCTION TO SPECIFY START SITE 
            REFSPLICESTOP1 = (int(EXON1.split("\t")[4]) - STARTPOSITION)
            REFSPLICENAME1 = re.split("-|\"", EXON1)[6]
            REFSPLICELEN1 = (REFSPLICESTOP1 - REFSPLICESTART1)
            REFSPLICECOMB1 = [REFSPLICENAME1,
                           REFSPLICESTART1,
                           REFSPLICELEN1,
                           REFSPLICESTOP1]
            CRLIST.append(REFSPLICECOMB1)

#Create a new function from the whole function. Make it smoother
#def

#Input command:
parser = argparse.ArgumentParser(description = "FLAME: Full Length Adjecency Matrix Enumeration") #CREATE A DEFAULT?!
parser.add_argument("-I", dest = "INPUT", help = "Input file")
parser.add_argument("-R", dest = "REF", help = "Reference File in Fasta format")
parser.add_argument("--range", dest = "RANGE", help = "Variance Range", default = 20)
parser.add_argument("-G", dest = "GTF", help = "Reference File in GTF format")
parser.add_argument("-B", dest = "SAM", help = "Shortread Sequence")
parser.add_argument("-O", dest = "OUTPUT1", help = "Correct Output file", default = "FLAME")
#parser.add_argument("-O2", dest = "OUTPUT2", help = "Faulty Output file")
args = parser.parse_args()
print("\n-----------------------------------------------------------------------------------")
print("\n-----------\tFLAME: Full Length Adjecency Matrix Enumeration\t\t-----------\n")

print("-----------------------------------------------------------------------------------")
print("\n-----------\tInitiating FLAME\t\t\t\t\t-----------")
print("Input:\t\t{}\n\
Shortread:\t{}\n\
Reference:\t{}\n\
GTF:\t\t{}\n\
Range:\t\t{}\n\
Output:\t\t{}-Suffix\n".format(
    args.INPUT,
    args.SAM,
    args.REF,
    args.GTF,
    args.RANGE,
    args.OUTPUT1
))
print("-----------------------------------------------------------------------------------\n")

#Variables:
##Input Files and References:
INPUT = open(args.INPUT, "r")
REFFILE = open(args.REF, "r")
REF = ""
GTFFILE = open(args.GTF, "r")
RANGESIZE = int(args.RANGE)
SHORTREAD = pysam.AlignmentFile("%s" %args.SAM, "rb")
NSTDLIST = [] #CHANGE NAME, Storage of the reference as a nested List.
WINDOWSIZE = 2 #FIX SO THAT THIS WINDOWSIZE IS FLEXIBLE
STARTPOSITION = 138373 #FIX SO THAT THIS STARTPOSITION IS FLEXIBLE
##FLAME: Filter Function:
FAULTYREADS = [] #Storage of the Faulty Reads, raw.
CORRECREADS = [] #Storage of the Correct Reads, raw.
SPLICESTARTALL1 = []
SPLICELENALL1 = []
SPLICECOUNT1 = 0
EXONCOUNTER1 = 0
Counter1 = 0
FAULTCOUNTER1 = 0
SPLICESTART1 = 0
SPLICELEN1 = 0
SPLICESTOP1 = 0
SPLICECOMB1 = []
RANGE1 = []
##FLAME: Translate Function, Correct:
SPLICESTARTALL2 = []
SPLICELENALL2 = []
SPLICECOUNT2 = 0
Counter2 = 0
TMPSTRING1 = "" #Temporary String to be filled with Translated Read Only.
SPLICESTART2 = 0
SPLICELEN2 = 0
SPLICESTOP2 = 0
SPLICECOMB2 = []
RANGE2 = []
TRANSLATECORREC = [] #Storage of the Correct Reads, Translated.
##FLAME: Translate Function, Faulty:
SPLICESTARTALL3 = []
SPLICELENALL3 = []
SPLICECOUNT3 = 0
EXONCOUNTER2 = 0
TMPSTRING2 = "" #Temporary String to be filled with Translated Read + Faulty Ranges.
Counter3 = 0
FAULTCOUNTER2 = 0
SPLICESTART3 = 0
SPLICELEN3 = 0
SPLICESTOP3= 0
SPLICECOMB3 = []
RANGE3 = []
TRANSLATEFAULTY = [] #Storage of the Faulty Reads, Translated
TRUENESS = False #Checkpoint to see if all the exons have been processed.
##FLAME: Quantification, Correct:
QUANTIFY = {} #Storage of the existing Exon Permutations (Key) and their quantification (Value), Correct.
##FLAME: Quantification, Faulty:
FQUANTIFY = {} #Storage of the existing Exon Permutations (Key) and their quantification (Value), Faulty.
##FLAME: Creation of the Adjecency Matrix, Correct:
LENNSTD = 0
ADJMTX1 = [] #Storage of the Adjecency Matrix, Correct, and their quantification. Basically a storage of Exon connections.
NUMBER = 0
Counter4 = 0
Counter5 = 0
TMPSPLICECOMB1 = "" #Temporary String to be filled with ExonN (The "Start Exon")
TMPSPLICECOMB2 = "" #Temporary String to be filled with ExonN+1 (The "End Exon")
##FLAME: Detection of Novel Splice/Exons based on Frequency, Part 1: Splice each potential new exons into Sequence Ranges:
EXONNAME = [] #Creating a "Blacklist" to filter out the Translated "Correct" reads from the Faulty Reads leaving only the Faulty Ranges.
SPLICEALL = [] #Breaking the Faulty reads into individual "Exons"
POTNEWEXON = [] #Storage for all the "Exonic" regions that does not exist in the reference.
##FLAME: Detection of Novel Splice Exons based on Frequency, Part 2: Quantification and Frequency Analysis of the Start+Stop:
GENEREFERENCE = [] #Create a list containing equal number of elements as nucleotides according to the reference.
COUNT9RANGE = [] #Temporary storage for the "Exonic" range.
STARTPOS = [] #Storage of the Variance/Range function, Start of the Exon.
ENDPOS = [] #Storage of the Variance/Range function, End of the Exon.
##FLAME: Detection of Novel Splice Exons based on Frequency, Part 3: Analysis of possible adjecent Splice Sequence (GU+AG):
SPLICEGU = False #Flag to detect if there exist a "GU" signal adjecently.
SPLICEAG = False #Flag to detect if there exist a "AG" signal adjecently.
CHECKPOINT = [] 
RANGE4 = []
##FLAME: Detection of Novel Splice Exons based on Frequency, Part 3.5: Storing of Candidates in a TSV as a Nested list.
Counter6 = 0
CANDIDATESPLICE = [] #Storage of the potential novel exon splice sites.
THRESHOLD = float(0) #Setting the Threshold as 1% of all the Faulty Reads.
##FLAME: Detection of Novel Splice Exons based on Frequency, Part 4: Using Shortread to confirm the existence of splice sites.

##FLAME: Detection of Novel Splice Exons based on Frequency, Part 5: Adjecency Matrix Analysis of Faulty Reads?
Counter7 = 0 #Empty
ADJMTX2 = [] #Empty
Counter8 = 0 #Empty
Counter9 = 0 #Empty
Counter10 = 0 #Empty
##FLAME: Printing Files:
Counter11 = 0
Counter12 = 0
Counter13 = 0
Counter14 = 0

#Preapre the Fasta reference:
for REFLINE in REFFILE:
    if ">" in REFLINE:
        pass
    else:
        REF += REFLINE.rstrip()

#Prepare GTF reference:
CREATEREF(GTFFILE, NSTDLIST)


#Filter Function:
print("-----------\tInitiate Filter Function\t\t\t\t-----------")
for WHOLEREAD1 in INPUT.read().split("\n"):
    if len(WHOLEREAD1) >= 1:
        SPLICESTARTALL1 = WHOLEREAD1.split("\t")[11]
        SPLICELENALL1 = WHOLEREAD1.split("\t")[10]
        SPLICECOUNT1 = int(WHOLEREAD1.split("\t")[9])
        EXONCOUNTER1 = 0
        for COUNT1 in range(SPLICECOUNT1):
            Counter1 = 0
            FAULTCOUNTER1 = 0
            SPLICESTART1 = int(SPLICESTARTALL1.split(",")[COUNT1])
            SPLICELEN1 = int(SPLICELENALL1.split(",")[COUNT1])
            SPLICESTOP1 = (SPLICESTART1 + SPLICELEN1)
            SPLICECOMB1 = [SPLICESTART1,
                          SPLICELEN1,
                          SPLICESTOP1]
            while FAULTCOUNTER1 != len(NSTDLIST): #CHANGE SO THAT THIS IS DEPENDENT ON YOU REFERENCE SIZE
                RANGE1 = [list(range(NSTDLIST[Counter1][1]-RANGESIZE,
                                     NSTDLIST[Counter1][1]+RANGESIZE)),
                          list(range(NSTDLIST[Counter1][2]-RANGESIZE,
                                     NSTDLIST[Counter1][2]+RANGESIZE)),
                          list(range(NSTDLIST[Counter1][3]-RANGESIZE,
                                     NSTDLIST[Counter1][3]+RANGESIZE))]
                if ((SPLICECOMB1[0] in RANGE1[0]) and
                    (SPLICECOMB1[1] in RANGE1[1]) and
                    (SPLICECOMB1[2] in RANGE1[2])):
                    EXONCOUNTER1 += 1
                    Counter1 += 1
                    break
                elif ((SPLICECOMB1[0] not in RANGE1[0]) or
                      (SPLICECOMB1[1] not in RANGE1[1]) or
                      (SPLICECOMB1[2] not in RANGE1[2])):
                    FAULTCOUNTER1 += 1
                    Counter1 += 1
                    continue
        if EXONCOUNTER1 == SPLICECOUNT1:
            CORRECREADS.append(WHOLEREAD1)
        elif EXONCOUNTER1 != SPLICECOUNT1:
            FAULTYREADS.append(WHOLEREAD1)


#Translate Function, Correct:
print("-----------\tInitiate Translate Function, Corrected\t\t\t-----------")
for WHOLEREAD2 in CORRECREADS:
    SPLICESTARTALL2 = WHOLEREAD2.split("\t")[11]
    SPLICELENALL2 = WHOLEREAD2.split("\t")[10]
    SPLICECOUNT2 = int(WHOLEREAD2.split("\t")[9])
    Counter2 = 0
    COUNT2 = 0
    TMPSTRING1 = ""
    while COUNT2 < SPLICECOUNT2:
        SPLICESTART2 = int(SPLICESTARTALL2.split(",")[COUNT2])
        SPLICELEN2 = int(SPLICELENALL2.split(",")[COUNT2])
        SPLICESTOP2 = (SPLICESTART2 + SPLICELEN2)
        SPLICECOMB2 = [SPLICESTART2,
                       SPLICELEN2,
                       SPLICESTOP2]
        RANGE2 = [list(range(NSTDLIST[Counter2][1]-RANGESIZE,
                             NSTDLIST[Counter2][1]+RANGESIZE)),
                  list(range(NSTDLIST[Counter2][2]-RANGESIZE,
                             NSTDLIST[Counter2][2]+RANGESIZE)),
                  list(range(NSTDLIST[Counter2][3]-RANGESIZE,
                             NSTDLIST[Counter2][3]+RANGESIZE))]
        if ((SPLICECOMB2[0] in RANGE2[0]) and
            (SPLICECOMB2[1] in RANGE2[1]) and
            (SPLICECOMB2[2] in RANGE2[2])):
            TMPSTRING1 += str(NSTDLIST[Counter2][0] + "-")
            COUNT2 += 1
            continue 
        elif ((SPLICECOMB2[0] not in RANGE2[0]) or
              (SPLICECOMB2[1] not in RANGE2[1]) or
              (SPLICECOMB2[2] not in RANGE2[2])):
            Counter2 += 1
            continue
    TRANSLATECORREC.append(TMPSTRING1)

    
#Translate Function: Faulty:
print("-----------\tInitiate Translate Function, Faulty\t\t\t-----------")
TRUENESSTCOUNT = 0
TRUENESSFCOUNT = 0
for WHOLEREAD3 in FAULTYREADS:
    SPLICESTARTALL3 = WHOLEREAD3.split("\t")[11]
    SPLICELENALL3 = WHOLEREAD3.split("\t")[10]
    SPLICECOUNT3 = int(WHOLEREAD3.split("\t")[9])
    EXONCOUNTER2 = 0
    TMPSTRING2 = ""
    for COUNT3 in range(SPLICECOUNT3):
        #TRUENESS = False
        Counter3 = 0
        FAULTCOUNTER2 = 0
        SPLICESTART3 = int(SPLICESTARTALL3.split(",")[COUNT3])
        SPLICELEN3 = int(SPLICELENALL3.split(",")[COUNT3])
        SPLICESTOP3 = (SPLICESTART3 + SPLICELEN3)
        SPLICECOMB3 = [SPLICESTART3,
                      SPLICELEN3,
                      SPLICESTOP3]
        while FAULTCOUNTER2 != len(NSTDLIST): #CHANGE SO THAT IT IS DEPENDENT ON REFERENCE SIZE
            RANGE3 = [list(range(NSTDLIST[Counter3][1]-RANGESIZE,
                                 NSTDLIST[Counter3][1]+RANGESIZE)),
                      list(range(NSTDLIST[Counter3][2]-RANGESIZE,
                                 NSTDLIST[Counter3][2]+RANGESIZE)),
                      list(range(NSTDLIST[Counter3][3]-RANGESIZE,
                                 NSTDLIST[Counter3][3]+RANGESIZE))]
            if ((SPLICECOMB3[0] in RANGE3[0]) and
                (SPLICECOMB3[1] in RANGE3[1]) and
                (SPLICECOMB3[2] in RANGE3[2])):
                EXONCOUNTER2 += 1
                TMPSTRING2 += str(NSTDLIST[Counter3][0] + ",")
                Counter3 += 1
                TRUENESS = False
                break
            elif ((SPLICECOMB3[0] not in RANGE3[0]) or
                  (SPLICECOMB3[1] not in RANGE3[1]) or
                  (SPLICECOMB3[2] not in RANGE3[2])):
                Counter3 += 1
                FAULTCOUNTER2 += 1
                TRUENESS = True
                continue
        #print(TRUENESS)
        if TRUENESS:
            TMPSTRING2 += str(str(SPLICECOMB3[0]) + "-" + str(SPLICECOMB3[2]) + ",")
            TRUENESSTCOUNT += 1
        elif TRUENESS == False:
            TRUENESSFCOUNT += 1
            pass
    TRANSLATEFAULTY.append(TMPSTRING2)

print(TRUENESSTCOUNT, TRUENESSFCOUNT)

#Quantification, Correct Reads:
print("-----------\tInitiate Quantification, Correct\t\t\t-----------")
for COUNT4 in TRANSLATECORREC:
    if COUNT4 in QUANTIFY.keys():
        QUANTIFY[COUNT4] += 1
    if COUNT4 not in QUANTIFY.keys():
        QUANTIFY[COUNT4] = 1


#Quantification, Faulty Reads
print("-----------\tInitiate FQuantification, Faulty\t\t\t-----------")
for COUNT4 in TRANSLATEFAULTY:
    if COUNT4 in FQUANTIFY.keys():
        FQUANTIFY[COUNT4] += 1
    if COUNT4 not in FQUANTIFY.keys():
        FQUANTIFY[COUNT4] = 1


#Create Adjecency Matrix, Correct:
print("-----------\tInitiate Creation of Empty Adjecency Matrix\t\t-----------")
LENNSTD = (len(NSTDLIST)+1)
for COUNT5 in range(LENNSTD): #Create an empty matrix function.
    ADJMTX1.append([0]*LENNSTD)
print("-----------\tCreation of Adjecency Matrix, Correct\t\t\t-----------")
for WHOLEREAD3 in TRANSLATECORREC:
    #print(WHOLEREAD3)
    NUMBER = (len(WHOLEREAD3.split("-"))-1)
    Counter4 = 0
    Counter5 = 0
    #The counter for adjecency matrix:
    for COUNT62 in range(NUMBER-1): #This for loop is to create the individual combinations of Start, Length, Stop.
        TMPSPLICECOMB1 = str(WHOLEREAD3.split("-")[COUNT62])
        TMPSPLICECOMB2 = str(WHOLEREAD3.split("-")[COUNT62+1])
        Counter4 = 0
        Counter5 = 0
        for COUNT63 in NSTDLIST: #This for loop is to find out which adjecency matrix position the Splice end site is located.
            if (TMPSPLICECOMB1 == COUNT63[0]):
                break
            elif (TMPSPLICECOMB1 != COUNT63[0]):
                Counter4 += 1
        for COUNT64 in NSTDLIST: #This for loop is to find out which adjecency matrix position the Splice start site is located,
            if (TMPSPLICECOMB2 == COUNT64[0]):
                break
            elif (TMPSPLICECOMB2 != COUNT64[0]):
                Counter5 += 1
        ADJMTX1[Counter5][Counter4] += 1


#Detection of Novel Splice/Exons based on Frequency:
print("-----------\tFrequency Analysis of Faulty Regions\t\t\t-----------")
##Creating a database for all the Exonic names:
for EXONNAME1 in NSTDLIST:
    EXONNAME.append(EXONNAME1[0])
##The quantifier and frequency analysis:
for WHOLEREAD4 in TRANSLATEFAULTY:
    SPLICEALL1 = WHOLEREAD4.split(",")[:-1]
    for COUNT7 in SPLICEALL1:
        if COUNT7 in EXONNAME:
            pass
        elif COUNT7 not in EXONNAME:
            POTNEWEXON.append(COUNT7)
print(POTNEWEXON)
'''
##Option 1: Sliding Window.
GENERANGE = [NSTDLIST[0][1],
             NSTDLIST[-1][3]]

print("-----------\tCreation of Adjecency Matrix, Faulty\t\t-----------")
#Create Adjecency Matrix, Faulty:
for COUNT8 in range(GENERANGE[0], GENERANGE[1], WINDOWSIZE): #CHANGE THE WINDOW SIZE SO IT IS CUSTOMIZEABLE
    ADJMTX2.append([0]*(len(range(GENERANGE[0], GENERANGE[1], WINDOWSIZE))))
#GENEWINDOW = (range(GENERANGE[0], GENERANGE[1], 100)) #CHANGE THE WINDOW SIZE SO IT IS CUSTOMIZEABLE
#print(GENEWINDOW1)
#print(len(GENEWINDOW1))

#for i in GENEWINDOW1:
#    print(len(i), "AAAAAAAAAAAAA")
#print(len(GENEWINDOW1), "BBBBBBBBBB")

#Fill the Adjecency Matrix, Faulty:
for COUNT9 in POTNEWEXON:
    COUNT9RANGE = COUNT9.split("-")
    Counter8 = 0
    Counter9 = 0
    for GENEWINDOW in range(len(ADJMTX2)-1):
        #print("A",  GENEWINDOW2, "B", GENE)
        if (GENEWINDOW*WINDOWSIZE) <= int(COUNT9RANGE[0]) < ((GENEWINDOW+ 1)*WINDOWSIZE): #CHANGE THE WINDOW SIZE SO IT IS CUSTOMIZEABLE
            #print(Counter9, "BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
            break
        else:
            Counter8 += 1
    for GENEWINDOW in range(len(ADJMTX2)-1):
        if (GENEWINDOW*WINDOWSIZE) <= int(COUNT9RANGE[1]) < ((GENEWINDOW + 1)*WINDOWSIZE): #CHANGE THE WINDOW SIZE SO IT IS CUSTOMIZEABLE
            #print(Counter9, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            break
        else:
            Counter9 += 1
    #print(Counter8, Counter9)
    ADJMTX2[Counter9][Counter8] += 1
#for COUNT10 in COUNT9RANGE:
#print(ADJMTX2)
#print(len(ADJMTX2))
'''

#Option 2: Range + Frequency:
#Novel Splice Detection:
print("-----------\tNovel Exon Detection: Frequency & Range\t\t\t-----------")
for COUNT8 in range(int(NSTDLIST[-1][3])+RANGESIZE):
    GENEREFERENCE.append(0)
for COUNT9 in POTNEWEXON:
    COUNT9RANGE = COUNT9.split("-") #This splits the two numbers (Start and End). This makes it so one cannot see the connections.
    if int(COUNT9RANGE[1]) < int(NSTDLIST[-1][3])+RANGESIZE:
        STARTPOS = list(range(int(COUNT9RANGE[0])-WINDOWSIZE,
                              int(COUNT9RANGE[0])+WINDOWSIZE))
        ENDPOS = list(range(int(COUNT9RANGE[1])-WINDOWSIZE,
                            int(COUNT9RANGE[1])+WINDOWSIZE))
        for COUNT10 in STARTPOS:
            GENEREFERENCE[COUNT10] += 1
        for COUNT10 in ENDPOS:
            GENEREFERENCE[COUNT10] += 1
    else:
        pass
        #print("Discarded Read")

##Create a TSV with: 1. Splice Position, 2. Supporting Reads (Faulty), 3. Supporting Percentage (Faulty)
for COUNT11 in GENEREFERENCE:
    THRESHOLD = float(float(0.01)*len(FAULTYREADS)) #Make this so one can specify the Threshold?
    if COUNT11 > THRESHOLD:
        CANDIDATESPLICE.append([Counter6, COUNT11, round((float(COUNT11)/len(FAULTYREADS))*100, 2)])
        Counter6 += 1
    else:
        Counter6 += 1


print("-----------\tConfirmation of Novel Splice Junctions: Splice Site\t-----------")
##Splice Confirmation1: Detection of Splice sites. They are also added to the previous [list].
for COUNT12 in CANDIDATESPLICE:
    SPLICEGU = False
    SPLICEAG = False
    CHECKPOINT1 = REF[((int(COUNT12[0])+STARTPOSITION)-3):((int(COUNT12[0])+STARTPOSITION)+4)]
    #print(COUNT12[0], REF[int(COUNT12[0])+STARTPOSITION], CHECKPOINT1, CHECKPOINT1[0:4], CHECKPOINT1[3:7]) #-------->FIX!<--------
    if "GT" in CHECKPOINT1[3:7]:
        SPLICEGU = True
    if "AG" in CHECKPOINT1[0:4]:
        SPLICEAG = True
    if SPLICEGU == False and SPLICEAG == False:
        COUNT12.append("None")
    elif SPLICEGU == True and SPLICEAG == False:
        COUNT12.append("GU")
    elif SPLICEGU == False and SPLICEAG == True:
        COUNT12.append("AG")
    elif SPLICEGU == True and SPLICEAG == True:
        COUNT12.append("Both")


print("-----------\tConfirmation of Novel Splice Junctions: Shortread\t-----------")
#Confirmation using short read sequencing, Part 1 - Splice Junction Detection:
SHORTREADSPLICE = {}
for COUNT13 in SHORTREAD:
    SHORTSTART = int(str(COUNT13).split("\t")[3])
    SHORTCIGAR = str(COUNT13).split("\t")[5]
    if "N" in SHORTCIGAR:
        SHORTLEN = sum(list(map(int, re.findall(r'\d+', SHORTCIGAR))))
        SHORT = [SHORTSTART,
                 SHORTCIGAR,
                 SHORTSTART+SHORTLEN]
        if (NSTDLIST[0][1]+STARTPOSITION) < SHORT[0] < (NSTDLIST[-1][3]+STARTPOSITION):
            MATCH1 = None
            MATCH2 = None
            MATCH3 = None
            MATCH4 = None
            if re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                        SHORT[1],
                        re.I):
                MATCH1 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                  SHORT[1],
                                  re.I)
            elif re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                          SHORT[1],
                          re.I):
                MATCH2 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                  SHORT[1],
                                  re.I)
            elif re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                          SHORT[1],
                          re.I):
                MATCH3 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                  SHORT[1],
                                  re.I)
            elif re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                          SHORT[1],
                          re.I):
                MATCH4 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                  SHORT[1],
                                  re.I)
            if MATCH1:
                ITEMS = list(MATCH1.groups())
            elif MATCH2:
                ITEMS = list(MATCH2.groups())
            elif MATCH3:
                ITEMS = list(MATCH3.groups())
            elif MATCH4:
                ITEMS = list(MATCH4.groups())
            else:
                print("EEEEEEEEEEEEEEEEEEEEEEEEEEEE", SHORT[1])
            if ITEMS[1] == "S":
                del ITEMS[:2]
            if ITEMS[-1] == "S":
                del ITEMS[-2:]            
            TMPVALUE1 = SHORTSTART
            for COUNT14 in range(len(ITEMS)):
                CIGARCOUNT = 0
                CIGAROPERATOR = 0
                if COUNT14 % 2 != 0 and COUNT14 > 0:
                    CIGARCOUNT = int(ITEMS[COUNT14-1])
                    CIGAROPERATOR = ITEMS[COUNT14]
                    if CIGAROPERATOR != "N":
                        TMPVALUE1 = TMPVALUE1 + CIGARCOUNT
                        pass
                    elif CIGAROPERATOR == "N":
                        TMPVALUE2 = TMPVALUE1-STARTPOSITION
                        TMPVALUE1 = TMPVALUE1 + CIGARCOUNT
                        if TMPVALUE2 in SHORTREADSPLICE:
                            SHORTREADSPLICE[TMPVALUE2] += 1
                        elif TMPVALUE2 not in SHORTREADSPLICE:
                            SHORTREADSPLICE[TMPVALUE2] = 1
                        TMPVALUE2 = TMPVALUE2 + CIGARCOUNT
                        if TMPVALUE2 in SHORTREADSPLICE:
                            SHORTREADSPLICE[TMPVALUE2] += 1
                        elif TMPVALUE2 not in SHORTREADSPLICE:
                            SHORTREADSPLICE[TMPVALUE2] = 1
                else:
                    pass
                
#Confirmation using short read sequencing, Part 2 - Crossreferencing with Longread:
print("-----------\tFaulty Adjecency Matrix?\t\t\t\t-----------")
for COUNT15 in CANDIDATESPLICE:
    Counter7 = 0
    for k, v in SHORTREADSPLICE.items():
        if COUNT15[0] == k:
            COUNT15.extend(["Yes", v])
            break
        elif COUNT15[0] != k:
            Counter7 += 1
    if Counter7 == len(SHORTREADSPLICE):
        COUNT15.extend(["No", "N/A"])

'''
#Use the new table to generate a adjecency matrix?
print(CANDIDATESPLICE[-1])
LENPOTNEWEXON = (len(CANDIDATESPLICE)+1)
for COUNT16 in range(LENPOTNEWEXON): #Create an empty matrix function.
    ADJMTX2.append([0]*LENPOTNEWEXON)
#print(CANDIDATESPLICE)
for COUNT17 in POTNEWEXON:
    Counter8 = 0
    Counter9 = 0
    TMPPOTEXONRANGE1 = list(range(int(COUNT17.split("-")[0])-RANGESIZE,
                                  int(COUNT17.split("-")[0])+RANGESIZE))
    TMPPOTEXONRANGE2 = list(range(int(COUNT17.split("-")[1])-RANGESIZE,
                                  int(COUNT17.split("-")[1])+RANGESIZE))
    for COUNT181 in CANDIDATESPLICE:
        #print(COUNT181[0], POTEXON)
        if COUNT181[0] in TMPPOTEXONRANGE1:
            #print("AAAAAAAAAAAAAAAA")
            break
        elif COUNT181[0] not in TMPPOTEXONRANGE1:
            Counter8 += 1
            #print("BBBBBBBBBBBBBBB")
    for COUNT182 in CANDIDATESPLICE:
        if COUNT182[0] in TMPPOTEXONRANGE2:
            #print("CCCCCCCCCCCCCCC")
            break
        elif COUNT182[0] not in TMPPOTEXONRANGE2:
            Counter9 += 1
    ADJMTX2[Counter9][Counter8] += 1
'''
            
#Print Functions:
print("-----------\tPrinting Files\t\t\t\t\t\t-----------")

##Print out Correct:
OUTPUT = open("%s.Correct.bed" %args.OUTPUT1, "w+")
for i in CORRECREADS:
    OUTPUT.write(str(i) + "\n")
OUTPUT.close()

##Print out Faulty:
OUTPUT = open("%s.Faulty.bed" %args.OUTPUT1, "w+")
for i in FAULTYREADS:
    OUTPUT.write(str(i) + "\n")
OUTPUT.close()

##Print out quantification:
OUTPUT = open("%s.Quantification.txt" %args.OUTPUT1, "w+")
for k, v in QUANTIFY.items():
    OUTPUT.write(str(v) +
                 "\t" +
                 str(k) +
                 "\n")
OUTPUT.close()

##Print out FAULTYquantification:
OUTPUT = open("%s.QuantificationF.txt" %args.OUTPUT1, "w+") #Change Name of the Output file. ",Faulty"?
for k, v in FQUANTIFY.items():
    OUTPUT.write(str(v) +
                 "\t" +
                 str(k) +
                 "\n")
OUTPUT.close()

##Print out Adjecency Matrix, Correct:
###Print the Column- and Rowheaders.
OUTPUT = open("%s.AdjecencyCorrect.tsv" %args.OUTPUT1, "w+")
for i in ADJMTX1:
    if Counter11 <= (len(NSTDLIST)-1):
        #print(Counter6)
        OUTPUT.write("\t" + str(NSTDLIST[Counter11][0]))
        Counter11 += 1
    elif Counter11 > (len(NSTDLIST)-1):
        OUTPUT.write("\t" + "END" + "\t")
OUTPUT.write("\n")
###Filling the Adjecency Matrix itself.
for i in ADJMTX1:
    if Counter12 <= (len(NSTDLIST)-1):
        OUTPUT.write(str(NSTDLIST[Counter12][0]) + "\t")
        for j in i:
            OUTPUT.write(str(j) + "\t")
        OUTPUT.write("\n")
        Counter12 += 1
    elif Counter12 > (len(NSTDLIST)-1):
        pass
OUTPUT.close()

##Print out Faulty:
OUTPUT = open("%s.FaultyTrans.txt" %args.OUTPUT1, "w+")
for i in TRANSLATEFAULTY:
    OUTPUT.write(str(i) + "\n")
OUTPUT.close

##Print out Potential Splice Sites:
OUTPUT = open("%s.PotentialSplice.tsv" %args.OUTPUT1, "w+")
OUTPUT.write("Gene Position" +
             "\t" +
             "Supporting Faulty Reads, Absolute" +
             "\t" +
             "Supporting Faulty Reads, Percent" +
             "\t" +
             "Adjecent Splice Signal" +
             "\t" +
             "Short Read Support" +
             "\t" +
             "Number of Supporting Short Reads" +
             "\n")
for i in CANDIDATESPLICE:
    OUTPUT.write(str(i[0]) +
                 "\t" +
                 str(i[1]) +
                 "\t" +
                 str(i[2]) +
                 "\t" +
                 str(i[3]) +
                 "\t" +
                 str(i[4]) +
                 "\t" +
                 str(i[5]) +
                 "\n")
OUTPUT.close

##Print out quantification:
OUTPUT = open("%s.Shortread.txt" %args.OUTPUT1, "w+")
OUTPUT.write("Number of Supporting Shortreads" +
             "\t" +
             "Gene Position" +
             "\n")
for k, v in SHORTREADSPLICE.items():
    OUTPUT.write(str(v) +
                 "\t" +
                 str(k) +
                 "\n")
OUTPUT.close()

##Print out Adjecency Matrix, Faulty:
###Print the Column- and Rowheaders.
OUTPUT = open("%s.AdjecencyFaulty.tsv" %args.OUTPUT1, "w+")
for i in ADJMTX2:
    if Counter13 <= (len(CANDIDATESPLICE)-1):
        #print(Counter6)
        OUTPUT.write("\t" + str(CANDIDATESPLICE[Counter13][0]))
        Counter13 += 1
    elif Counter13 > (len(CANDIDATESPLICE)-1):
        OUTPUT.write("\t" + "END" + "\t")
OUTPUT.write("\n")
###Filling the Adjecency Matrix itself.
for i in ADJMTX2:
    if Counter14 <= (len(CANDIDATESPLICE)-1):
        OUTPUT.write(str(CANDIDATESPLICE[Counter14][0]) + "\t")
        for j in i:
            OUTPUT.write(str(j) + "\t")
        OUTPUT.write("\n")
        Counter14 += 1
    elif Counter14 > (len(CANDIDATESPLICE)-1):
        pass
OUTPUT.close()

OUTPUT = open("%s.Asdf" %args.OUTPUT1, "w+")
for i in POTNEWEXON:
    OUTPUT.write(str(i) + "\n")
OUTPUT.close()
