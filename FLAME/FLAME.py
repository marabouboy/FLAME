#!/usr/bin/python3
#FLAME 0.1.4
#Import Packages:
import argparse, io
import pysam
import FLAME_FUNC.FLAME_FUNC as FF #Break this down into each part based on the command.

#####################
# Need to be fixed List:
# Fix the comments within FLAME.py. It is a bit wonkey with different explainations for different functions. Be consistent.
#
#####################

#Input command:
parser = argparse.ArgumentParser(description = "FLAME: Full Length Adjacency Matrix Enumeration")
##Obligatory Inputs:
parser.add_argument("-I", dest = "INPUT", help = "Input file [Required]")
parser.add_argument("-GTF", dest = "GTF", help = "Reference File in GTF format")
parser.add_argument("-G", dest = "GENE", help = "Targeted Gene [Required]", default = "Transcriptome Mode")
parser.add_argument("--range", dest = "RANGE", help = "Variance Range", default = 20)
parser.add_argument("--min", dest = "MINIMUM", help = "Minimum Read Coverage", default = 10)
parser.add_argument("--ratio", dest = "RATIO", help = "Minimum Annotation Ratio", default = float(0.25))
parser.add_argument("-O", dest = "OUTPUT", help = "Output Prefix", default = "Flame")
##Optional Inputs:
parser.add_argument("-R", dest = "REF", help = "Reference File in Fasta format", default = 0)
#parser.add_argument("-THRESHOLD?", dest = "THRESHOLD", help = "Threshold for the Frequency Analysis, the lower threshold the more comprehensive but the longer processing power required", default = 0.01)?
parser.add_argument("-B", dest = "SAM", help = "Shortread Sequence", default = 0)
parser.add_argument("--verbose", dest = "VERBOSE", help = "Verbose output", action="store_true")

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
    args.RANGE, #ADD ANOTHER FILED SPECIFYING THE MINIMUM RATIO FOR FLAME-TRANSCRIPTOMEWIDE?
    args.OUTPUT
))


#Variables:
##Input Files and References:
Flame_Main_INPUT1 = []
for Flame_Main_COUNT in open(args.INPUT, "r"):
    Flame_Main_INPUT1.append(Flame_Main_COUNT)
try:
    if args.GTF != None:
        Flame_Main_GTFFILE = open(args.GTF, "r") #Storage[List] of the GTF file.
    else:
        Flame_Main_GTFFILE = io.StringIO("")
except:
    print("GTF-File Error")

Flame_Main_RANGESIZE = int(args.RANGE) #Storage(Integer) of the Windowsize.
Flame_Main_MINIMUMREADS = int(args.MINIMUM) #Storage(Integer) of the minimum number of reads.
Flame_Main_RATIOTHRESH = float(args.RATIO) #Storage(Float) of the ratio required for a gene to pass through and be Translated and Quantified.
Flame_Main_GENENAME = args.GENE #Storage("String") of the Specific Gene.
Flame_Main_REFERENCELIST = [] #Storage[Nested List] of the reference.
Flame_Main_SIGNIFGENES = []
Flame_Main_FREQUENCYWINDOWSIZE = 2 #FIX: Make this interactive so one can input the window size for the Frequency Analysis.
Flame_Main_FREQYENCYTHRESHOLD = float(0.01) #FIX: Make this interactive so one can input the window size for the Frequency Threshold.

#-------------------------------------------------------------------------------------------------------------------------------------#
#FLAME: Transcriptome Mode:
if args.GENE == "Transcriptome Mode": #<----------------- Change here
    Flame_Main_GENELIST1 = []
    Flame_Main_GENELIST2 = []
    Flame_Main_TMPGENE1 = ""
    Flame_Main_TMPGENE2 = ""
    Flame_Main_INPUT2 = []
    Flame_Main_PROGRESSMAX = len(Flame_Main_INPUT1)
    Flame_Main_Counter1 = True
    Flame_Main_Counter2 = True
    Flame_Main_Counter3 = True
    
    #FLAME-WT: 1. Segment the number of genes into a list:
    for Flame_Main_COUNT1 in Flame_Main_GTFFILE.read().split("\n"):
        try:
            if (any(Flame_Main_COUNT1.split("\t")[2]) and
                ("exon" in Flame_Main_COUNT1.split("\t")[2])):
                try:
                    Flame_Main_TMPGENE1 = Flame_Main_COUNT1.split(";")[2].split(" ")[2][1:-1]
                except: #Need this "extra" except clause?
                    pass 
                if Flame_Main_TMPGENE1 == Flame_Main_TMPGENE2:
                    pass
                else:
                    Flame_Main_TMPGENE1 = Flame_Main_TMPGENE1
                    Flame_Main_GENELIST1.append(Flame_Main_TMPGENE1)
            else:
                pass
        except:
            pass
    Flame_Main_GENELIST2 = list(set(Flame_Main_GENELIST1))
    Flame_Main_GENELIST2.sort()

    #FLAME-WT: 2. Loop through this list in order to create a tmp-ref-file:
    for Flame_Main_COUNT1 in Flame_Main_GENELIST2:
        Flame_Main_GTFFILE.seek(0, 0)
        print("\n%s:" %Flame_Main_COUNT1)
        print("-----------\tInitiate Creation of Reference\t\t\t\t\t-----------")
        Flame_Main_REFERENCELIST = FF.CREATEREFFUNC(Flame_Main_GTFFILE,
                                                    Flame_Main_COUNT1)
        print("-----------\tInitiate Extraction of Relevant Reads\t\t\t\t-----------")
        Flame_Main_INPUT2 = FF.SEGMENTFUNC(Flame_Main_INPUT1,
                                           Flame_Main_REFERENCELIST)
        
        #FLAME-WT: 3. Run the FILTERFUNC on the seperated reads:
        print("-----------\tInitiate Filter Function\t\t\t\t\t-----------")
        Flame_Main_ANNOTATEDREADS, Flame_Main_INCONGRUENTREADS = FF.FILTERFUNC(Flame_Main_INPUT2,
                                                                               Flame_Main_REFERENCELIST,
                                                                               Flame_Main_RANGESIZE)
        
        #FLAME-WT: 4. Acquire the Ratio of Annoated Reads vs Incongruent Reads:
        try:
            Flame_Main_RATIO = float(len(Flame_Main_ANNOTATEDREADS)/
                                 (len(Flame_Main_ANNOTATEDREADS)+len(Flame_Main_INCONGRUENTREADS)))
        except (ZeroDivisionError):
            print("%s: Zero Division, No Read-Coverage" % Flame_Main_COUNT1)
            continue
        
        #5a. If AnnotatedReads/TotalReads >= 0.25 (25%); simply run it through TRANSLATEFUNC & QUANTFUNC.
        #print("Number of Reads:", len(Flame_Main_INPUT2))<------
        if len(Flame_Main_INPUT2) < Flame_Main_MINIMUMREADS:
            print("%s: Number of reads under gene-region below minimum-read threshold." %Flame_Main_COUNT1)
            continue
        if Flame_Main_RATIO >= Flame_Main_RATIOTHRESH:
            Flame_Main_TRANSLATEANNOTATED, Flame_Main_TRANSLATEINCONGRUENT = FF.TRANSLATEFUNC(Flame_Main_REFERENCELIST,
                                                                                              Flame_Main_RANGESIZE,
                                                                                              Flame_Main_ANNOTATEDREADS)
            print("-----------\tInitiate Quantification of each Splice Permutaiton, Annotated\t-----------")
            Flame_Main_QUANTIFYANNOTATED = FF.QUANTIFYFUNC(Flame_Main_TRANSLATEANNOTATED)
            Flame_Main_READTOTALLENGTH = []
            ##Print out quantification:
            if Flame_Main_Counter1:
                Flame_Main_OUTPUT = open("%s.QuantAnnotated.tsv" %args.OUTPUT, "w+")
                Flame_Main_OUTPUT.write("Gene" + "\t" +
                                        "Count" + "\t" +
                                        "Length of Isoform" + "\t" +
                                        "Number of Exons" + "\t" +
                                        "Isoform Permutation" + "\n")
                for k, v in Flame_Main_QUANTIFYANNOTATED.items():
                    #-----------The Function for Summing the Total Length of the Splice Variant-----------#
                    Flame_Main_READTOTALLENGTH = []
                    for Flame_Main_ANNOTATEDEXON in k.split(","):
                        for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
                            if Flame_Main_ANNOTATEDEXON == Flame_Main_COUNT[0]:
                                Flame_Main_READTOTALLENGTH.append(int(Flame_Main_COUNT[3]))
                            elif Flame_Main_ANNOTATEDEXON != Flame_Main_COUNT[0]:
                                pass
                    #-----------The Function for Summing the Total Length of the Splice Variant-----------#
                    Flame_Main_OUTPUT.write(str(Flame_Main_COUNT1) + "\t" +
                                            str(v) + "\t" +
                                            str(sum(Flame_Main_READTOTALLENGTH)) + "\t" +
                                            str(len(Flame_Main_READTOTALLENGTH)) + "\t" +
                                            str(k) + "\n")
                Flame_Main_OUTPUT.close()
                Flame_Main_Counter1 = False
            else:
                if Flame_Main_QUANTIFYANNOTATED != {}:
                    Flame_Main_OUTPUT = open("%s.QuantAnnotated.tsv" %args.OUTPUT, "a")
                    for k, v in Flame_Main_QUANTIFYANNOTATED.items():
                        #-----------The Function for Summing the Total Length of the Splice Variant-----------#
                        Flame_Main_READTOTALLENGTH = []
                        for Flame_Main_ANNOTATEDEXON in k.split(","):
                            for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
                                if Flame_Main_ANNOTATEDEXON == Flame_Main_COUNT[0]:
                                    Flame_Main_READTOTALLENGTH.append(int(Flame_Main_COUNT[3]))
                                elif Flame_Main_ANNOTATEDEXON != Flame_Main_COUNT[0]:
                                    pass
                        #-----------The Function for Summing the Total Length of the Splice Variant-----------#
                        Flame_Main_OUTPUT.write(str(Flame_Main_COUNT1) + "\t" +
                                                str(v) + "\t" +
                                                str(sum(Flame_Main_READTOTALLENGTH)) + "\t" +
                                                str(len(Flame_Main_READTOTALLENGTH)) + "\t" +
                                                str(k) + "\n")
                    Flame_Main_OUTPUT.close()
                else:
                    pass
            if Flame_Main_Counter2:
                Flame_Main_OUTPUT = open("%s.Reference.txt" %args.OUTPUT, "w+")
                Flame_Main_OUTPUT.write("Gene" + "\t" +
                                        "Exon Name" + "\t" +
                                        "Chromosome"+ "\t" +
                                        "Exon Start Site" + "\t" +
                                        "Exon Length" + "\t" +
                                        "Exon Stop Site" + "\n")
                for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
                    Flame_Main_OUTPUT.write(str(Flame_Main_COUNT1)+ "\t" +
                                            str(Flame_Main_COUNT[0]) + "\t" +
                                            str(Flame_Main_COUNT[1]) + "\t" +
                                            str(Flame_Main_COUNT[2]) + "\t" +
                                            str(Flame_Main_COUNT[3]) + "\t" +
                                            str(Flame_Main_COUNT[4]) + "\n")
                Flame_Main_OUTPUT.close()
                Flame_Main_Counter2 = False
            else:
                Flame_Main_OUTPUT = open("%s.Reference.txt" %args.OUTPUT, "a")
                for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
                    Flame_Main_OUTPUT.write(str(Flame_Main_COUNT1)+ "\t" +
                                            str(Flame_Main_COUNT[0]) + "\t" +
                                            str(Flame_Main_COUNT[1]) + "\t" +
                                            str(Flame_Main_COUNT[2]) + "\t" +
                                            str(Flame_Main_COUNT[3]) + "\t" +
                                            str(Flame_Main_COUNT[4]) + "\n")
                Flame_Main_OUTPUT.close()
        else:
            #5b. If AnnotatedReads/TotalReads <= 0.25; flag the gene for having low consensus with the reference and recommend further inspection.
            Flame_Main_SIGNIFGENES.append(Flame_Main_COUNT1)
            if Flame_Main_Counter3:
                Flame_Main_OUTPUT = open("%s.SignificantGenes.txt" %args.OUTPUT, "w+")
                Flame_Main_OUTPUT.write("Gene" +  "\t" +
                                        "Number of Reads" + "\t"+
                                        "Ratio Annotated/Incongruent" + "\n")
                Flame_Main_OUTPUT.write(str(Flame_Main_COUNT1) + "\t" +
                                        str(len(Flame_Main_INPUT2)) + "\t" +
                                        str(round(Flame_Main_RATIO, 2)) + "\n")
                Flame_Main_OUTPUT.close()
                Flame_Main_Counter3 = False
            else:
                Flame_Main_OUTPUT = open("%s.SignificantGenes.txt" %args.OUTPUT, "a")
                Flame_Main_OUTPUT.write(str(Flame_Main_COUNT1) + "\t" +
                                        str(len(Flame_Main_INPUT2)) + "\t" +
                                        str(round(Flame_Main_RATIO, 2)) + "\n")
                Flame_Main_OUTPUT.close()

#-------------------------------------------------------------------------------------------------------------------------------------#
#FLAME: Single-Gene Mode:
else: #<----------------- Change here
    
    ##The optional choices:
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
        #If-statement that checks whether the Sequence Mapping/Alignment file is in SAM format or the compressed BAM format.
        if str(args.SAM).endswith(".sam") or str(args.SAM).endswith(".SAM"):
            Flame_Main_SHORTREAD1 = pysam.AlignmentFile("%s" %args.SAM, "r")
        elif str(args.SAM).endswith(".bam") or str(args.SAM).endswith(".BAM"):
            Flame_Main_SHORTREAD1 = pysam.AlignmentFile("%s" %args.SAM, "rb")
        else:
            print("Shortread input file format not recognized ([File].bam/[File].BAM/[File].sam/[File].SAM)")
            
    if args.VERBOSE:
        print("Verbose:\tOn\n")
    print("-------------------------------------------------------------------------------------------\n")

    ###The Optional Removal of both.
    if args.REF == 0 and args.SAM == 0:
        print("\n")

    ##Variables used for the printing functions:
    Flame_Main_Counter1 = 0
    Flame_Main_Counter2 = 0
    Flame_Main_Counter3_1 = False
    Flame_Main_Counter3_2 = False
    Flame_Main_Counter4 = 0
    Flame_Main_Counter5 = 0
    
    #Function to transform the GTF-file into an efficient python nested list: Exon[Name, Chromosome, Start, Length, Stop].
    Flame_Main_REFERENCELIST = FF.CREATEREFFUNC(Flame_Main_GTFFILE,
                                                args.GENE) #((The GTF Input file), (The name of the Gene))

    
    print("-----------\tInitiate Filter Function\t\t\t\t\t-----------")
    #Function to classify the reads as either "Annotated" or "Incongruent" and store them.
    Flame_Main_ANNOTATEDREADS, Flame_Main_INCONGRUENTREADS = FF.FILTERFUNC(Flame_Main_INPUT1,
                                                                           Flame_Main_REFERENCELIST,
                                                                           Flame_Main_RANGESIZE) #Input, Reference, Rangesize
    
    ##Print out Annotated:
    if Flame_Main_ANNOTATEDREADS != []:
        Flame_Main_OUTPUT = open("%s.Annotated.bed" %args.OUTPUT, "w+")
        for Flame_Main_COUNT in Flame_Main_ANNOTATEDREADS:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
        Flame_Main_OUTPUT.close()
    else:
        pass

    ##Print out Incongruent:
    if Flame_Main_INCONGRUENTREADS != []:
        Flame_Main_OUTPUT = open("%s.Incongruent.bed" %args.OUTPUT, "w+")
        for Flame_Main_COUNT in Flame_Main_INCONGRUENTREADS:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
        Flame_Main_OUTPUT.close()
    else:
        pass   

    #Translate Function. They have built in Print function. Also make a function that will make either Annotated/Incongruent reads = 0 to avoid downstream computation, if one chooses.
    Flame_Main_TRANSLATEANNOTATED, Flame_Main_TRANSLATEINCONGRUENT = FF.TRANSLATEFUNC(Flame_Main_REFERENCELIST,
                                                                                      Flame_Main_RANGESIZE,
                                                                                      Flame_Main_ANNOTATEDREADS,
                                                                                      Flame_Main_INCONGRUENTREADS)

    if args.VERBOSE:
        ##Print out Annotated:
        if Flame_Main_TRANSLATEANNOTATED != []:
            Flame_Main_OUTPUT = open("%s.AnnotatedTrans.txt" %args.OUTPUT, "w+")
            for Flame_Main_COUNT in Flame_Main_TRANSLATEANNOTATED:
                Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
            Flame_Main_OUTPUT.close()
        else:
            pass
    
    if args.VERBOSE:
        ##Print out Incongruent:
        if Flame_Main_TRANSLATEINCONGRUENT != []:
            Flame_Main_OUTPUT = open("%s.IncongruentTrans.txt" %args.OUTPUT, "w+")
            for Flame_Main_COUNT in Flame_Main_TRANSLATEINCONGRUENT:
                Flame_Main_OUTPUT.write(str(Flame_Main_COUNT) + "\n")
            Flame_Main_OUTPUT.close()
        else:
            pass

    #Annotated Adjacency Matrix Function:
    print("-----------\tInitiate Creation of Empty Adjacency Matrix\t\t\t-----------")
    #Command itself:
    Flame_Main_ADJMTX1 = FF.EMPTYADJMTXFUNC(Flame_Main_REFERENCELIST)
    print("-----------\tInitiate Creation of Adjacency Matrix, Annotated\t\t-----------")
    #Command itself:
    Flame_Main_ADJMTX1 = FF.ANNOTATEDADJMTXFUNC(Flame_Main_TRANSLATEANNOTATED,
                                                Flame_Main_REFERENCELIST,
                                                Flame_Main_ADJMTX1)
    ##Print out Adjacency Matrix, Annotated:
    ###Print the Column- and Rowheaders.
    Flame_Main_OUTPUT = open("%s.AdjacencyAnnotated.tsv" %args.OUTPUT, "w+")
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


    print("-----------\tInitiate Singling of Incongruent Exons\t\t\t\t-----------")
    #Command itself:
    Flame_Main_POTENTIALS = FF.INCONGRUENTSEPERATORFUNC(Flame_Main_TRANSLATEINCONGRUENT,
                                                        Flame_Main_REFERENCELIST)
    
    print("-----------\tInitiate Novel Splice Site Detection, Part1: Frequency\t\t-----------")
    #Command itself:
    Flame_Main_GENEREFERENCE = FF.FREQUENCYSITEFUNC(Flame_Main_POTENTIALS,
                                                    Flame_Main_REFERENCELIST,
                                                    Flame_Main_RANGESIZE,
                                                    Flame_Main_FREQUENCYWINDOWSIZE)
    
    print("-----------\tInitiate Novel Splice Site Detection, Part2: Threshold\t\t-----------")
    #Prepare GTF reference:
    Flame_Main_SPLICECANDIDATES = []
    #Command itself:
    Flame_Main_SPLICECANDIDATES = FF.FREQUENCYTHRESHFUNC(Flame_Main_GENEREFERENCE,
                                                         Flame_Main_FREQYENCYTHRESHOLD,
                                                         Flame_Main_INCONGRUENTREADS)
    
    if args.REF != 0:
        print("-----------\tInitiate Novel Splice Site Detection, Part3: Splice Signal\t-----------")
        #Command itself:
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
        Flame_Main_SPLICECANDIDATES, Flame_Main_SPLICESITECOUNT = FF.SHORTREADFUNC(Flame_Main_SHORTREAD1,
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
    ##################################################
    ##--->ADD THE CHROMOSOME WITHIN THE OUTPUT?!<---##
    ##################################################
    Flame_Main_OUTPUT = open("%s.PotentialSplice.tsv" %args.OUTPUT, "w+")
    if Flame_Main_Counter3_1 == False and Flame_Main_Counter3_2 == False:
        Flame_Main_OUTPUT.write("Gene Position" + "\t" +
                                "Supporting Incongruent Reads, Absolute" + "\t" +
                                "Supporting Incongruent Reads, Percent" + "\n")
        for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) + "\t" +
                                    str(Flame_Main_COUNT[1]) + "\t" +
                                    str(Flame_Main_COUNT[2]) + "\n")
        Flame_Main_OUTPUT.close
    elif Flame_Main_Counter3_1 == True and Flame_Main_Counter3_2 == False :
        Flame_Main_OUTPUT.write("Gene Position" + "\t" +
                                "Supporting Incongruent Reads, Absolute" + "\t" +
                                "Supporting Incongruent Reads, Percent" + "\t" +
                                "Adjacent Splice Signal" + "\n")
        for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) + "\t" +
                                    str(Flame_Main_COUNT[1]) + "\t" +
                                    str(Flame_Main_COUNT[2]) + "\t" +
                                    str(Flame_Main_COUNT[3]) + "\n")
        Flame_Main_OUTPUT.close
    elif Flame_Main_Counter3_1 == False and Flame_Main_Counter3_2 == True:
        Flame_Main_OUTPUT.write("Gene Position" + "\t" +
                                "Supporting Incongruent Reads, Absolute" + "\t" +
                                "Supporting Incongruent Reads, Percent" + "\t" +
                                "Short Read Support" + "\n")
        for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) + "\t" +
                                    str(Flame_Main_COUNT[1]) + "\t" +
                                    str(Flame_Main_COUNT[2]) + "\t" +
                                    str(Flame_Main_COUNT[3]) + "\n")
        Flame_Main_OUTPUT.close
    elif Flame_Main_Counter3_1 == True and Flame_Main_Counter3_2 == True:
        Flame_Main_OUTPUT.write("Gene Position" + "\t" +
                                "Supporting Incongruent Reads, Absolute" + "\t" +
                                "Supporting Incongruent Reads, Percent" + "\t" +
                                "Adjacent Splice Signal" + "\t" +
                                "Short Read Support" + "\n")
        for Flame_Main_COUNT in Flame_Main_SPLICECANDIDATES:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) + "\t" +
                                    str(Flame_Main_COUNT[1]) + "\t" +
                                    str(Flame_Main_COUNT[2]) + "\t" +
                                    str(Flame_Main_COUNT[3]) + "\t" +
                                    str(Flame_Main_COUNT[4]) + "\n")
        Flame_Main_OUTPUT.close
    else:
        print("Error, Type: Counter3")

    print("-----------\tInitiate Quantification of each Splice Permutaiton, Annotated\t-----------")
    #Command itself:
    Flame_Main_QUANTIFYANNOTATED = FF.QUANTIFYFUNC(Flame_Main_TRANSLATEANNOTATED)
    Flame_Main_READTOTALLENGTH = []
    ##Print out quantification:
    if Flame_Main_QUANTIFYANNOTATED != {}:
        Flame_Main_OUTPUT = open("%s.QuantAnnotated.tsv" %args.OUTPUT, "w+")
        Flame_Main_OUTPUT.write("Count" + "\t" +
                                "Length of Isoform" + "\t" +
                                "Number of Exons" + "\t" +
                                "Isoform Permutation" + "\n")
        for k, v in Flame_Main_QUANTIFYANNOTATED.items():
            #-----------The Function for Summing the Total Length of the Splice Variant-----------#
            Flame_Main_READTOTALLENGTH = []
            for Flame_Main_ANNOTATEDEXON in k.split(","):
                for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
                    if Flame_Main_ANNOTATEDEXON == Flame_Main_COUNT[0]:
                        Flame_Main_READTOTALLENGTH.append(int(Flame_Main_COUNT[3]))
                    elif Flame_Main_ANNOTATEDEXON != Flame_Main_COUNT[0]:
                        pass
            #-----------The Function for Summing the Total Length of the Splice Variant-----------#
            Flame_Main_OUTPUT.write(str(v) + "\t" +
                                    str(sum(Flame_Main_READTOTALLENGTH)) + "\t" +
                                    str(len(Flame_Main_READTOTALLENGTH)) + "\t" +
                                    str(k) + "\n")
        Flame_Main_OUTPUT.close()
    else:
        pass

    print("-----------\tInitiate Quantification of each Splice Permutaiton, Incongruent\t-----------")
    #Command itself:
    Flame_Main_QUANTIFYINCONGRUENT = FF.QUANTIFYFUNC(Flame_Main_TRANSLATEINCONGRUENT)
    ##Print out INCONGRUENT quantification:
    if Flame_Main_QUANTIFYINCONGRUENT != {}:
        Flame_Main_OUTPUT = open("%s.QuantIncongruent.tsv" %args.OUTPUT, "w+")
        Flame_Main_OUTPUT.write("Count" + "\t" +
                                "Isoform" + "\n")
        for k, v in Flame_Main_QUANTIFYINCONGRUENT.items():        
            Flame_Main_OUTPUT.write(str(v) + "\t" +
                                    str(k) + "\n")
        Flame_Main_OUTPUT.close()
    else:
        pass


    #Create an option choice of one wants to produce a Adjacency Matrix:
    print("-----------\tInitiate Creation of Empty Adjacency Matrix\t\t\t-----------")
    Flame_Main_ADJMTX2 = FF.EMPTYADJMTXFUNC(Flame_Main_SPLICECANDIDATES)
    print("-----------\tInitiate Creation of Adjacency Matrix, Incongruent\t\t-----------")
    Flame_Main_ADJMTX2 = FF.INCONGRUENTADJMTXFUNC(Flame_Main_POTENTIALS,
                                                  Flame_Main_SPLICECANDIDATES,
                                                  Flame_Main_ADJMTX2,
                                                  Flame_Main_RANGESIZE)
    ##Print out Adjacency Matrix, Incongruent:
    ###Print the Column- and Rowheaders.
    Flame_Main_OUTPUT = open("%s.AdjacencyIncongruent.tsv" %args.OUTPUT, "w+")
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
        Flame_Main_OUTPUT.write("Exon Name" + "\t" +
                                "Chromosome" + "\t" +
                                "Exon Start Site" + "\t" +
                                "Exon Length" + "\t" +
                                "Exon Stop Site" + "\n")
        for Flame_Main_COUNT in Flame_Main_REFERENCELIST:
            Flame_Main_OUTPUT.write(str(Flame_Main_COUNT[0]) + "\t" +
                                    str(Flame_Main_COUNT[1]) + "\t" +
                                    str(Flame_Main_COUNT[2]) + "\t" +
                                    str(Flame_Main_COUNT[3]) + "\t" +
                                    str(Flame_Main_COUNT[4]) + "\n")
        Flame_Main_OUTPUT.close()
