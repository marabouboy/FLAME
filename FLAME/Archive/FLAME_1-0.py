#!/usr/bin/python3
#Import Packages:
import re #time, sys, os, re
import argparse
import pysam
#CHANGE THIS SO IT FOLLOWS THE REST OF THE PROGRAM, HAVE A FUNCTION THAT ONE CAN USE THE WHOLE GENE, AS WELL AS HAVING ONE THAT FIGURES THE STARTPOSITIOIN BY ITSELF. FILTER AND EXTRACT ONLY THE "exon"'s WE'RE INTERESTED IN. ALSO BUILT IN NAMING SYSTEM?: 1ST LAYER: I-II-III-IV-V-VI..., 2ND LAYER: a-b-c-d-e-f-g..., 3RD LAYER: 1-2-3-4-5-6...?
def CREATEREF(CRINPUT, CRLIST, STARTPOSITION):
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


#Central Filter function that will sort the reads into two different lists (FILTERCORRECTOUTPUT, FILTERFAULTYOUTPUT) depending if they have only contain perfectly matching exons to the reference or if they have a "faulty" exon.
def FILTERFUNC(FILTERINPUT, REF, FILTERRANGESIZE):
    #-----------Outside Counters and variables-----------#
    Flame_Filter_EXONCOUNTER = 0 
    Flame_Filter_Counter = 0
    Flame_Filter_FAULTYCOUNTER = 0
    FILTERCORRECTOUTPUT = []
    FILTERFAULTYOUTPUT = []
    Flame_Filter_SINGLE_READ_STRT_ALL = []
    Flame_Filter_SINGLE_READ_LEN_ALL = []
    Flame_Filter_SINGLE_READ_COUNT = 0
    Flame_Filter_SPLICESTART = 0
    Flame_Filter_SPLICELEN = 0
    Flame_Filter_SPLICESTOP = 0
    Flame_Filter_SPLICECOMB = []
    Flame_Filter_RANGE = []
    #-----------The Function Itself-----------#
    for Flame_Filter_SINGLE_READ in FILTERINPUT.read().split("\n"): #Loop: Going through each read. #Change it so that the input read (As an outside read file) is read outside?
        if len(Flame_Filter_SINGLE_READ) >= 1: #If statement that filters out empty reads.
            #For ease of use (Below): Singling out all the relevant information from the bed12 read:
            Flame_Filter_SINGLE_READ_STRT_ALL = Flame_Filter_SINGLE_READ.split("\t")[11]
            Flame_Filter_SINGLE_READ_LEN_ALL =  Flame_Filter_SINGLE_READ.split("\t")[10]
            Flame_Filter_SINGLE_READ_COUNT = int(Flame_Filter_SINGLE_READ.split("\t")[9])
            Flame_Filter_EXONCOUNTER = 0
            for Flame_Filter_COUNT1 in range(Flame_Filter_SINGLE_READ_COUNT): #Loop: Going through each exon.
                Flame_Filter_Counter = 0
                Flame_Filter_FAULTYCOUNTER = 0
                #For ease of use (Below): Specifying each characteristic (Start, Length, Stop) for each exon in each read.
                Flame_Filter_SPLICESTART = int(Flame_Filter_SINGLE_READ_STRT_ALL.split(",")[Flame_Filter_COUNT1])
                Flame_Filter_SPLICELEN = int(Flame_Filter_SINGLE_READ_LEN_ALL.split(",")[Flame_Filter_COUNT1])
                Flame_Filter_SPLICESTOP = (Flame_Filter_SPLICESTART + Flame_Filter_SPLICELEN)
                #For ease of use (Below): Collapsing the characteristics into a single list.
                Flame_Filter_SPLICECOMB = [Flame_Filter_SPLICESTART,
                                           Flame_Filter_SPLICELEN,
                                           Flame_Filter_SPLICESTOP]
                while Flame_Filter_FAULTYCOUNTER != len(REF): #Loop: Going through each reference exon trying to match with the exon.
                    Flame_Filter_RANGE = [list(range(REF[Flame_Filter_Counter][1]-FILTERRANGESIZE,
                                                     REF[Flame_Filter_Counter][1]+FILTERRANGESIZE)),
                                          list(range(REF[Flame_Filter_Counter][2]-FILTERRANGESIZE,
                                                     REF[Flame_Filter_Counter][2]+FILTERRANGESIZE)),
                                          list(range(REF[Flame_Filter_Counter][3]-FILTERRANGESIZE,
                                                     REF[Flame_Filter_Counter][3]+FILTERRANGESIZE))] #Variance Function (Meant to account for experimental, sequencing, base-calling, and/or alignment errors.)
                    if ((Flame_Filter_SPLICECOMB[0] in Flame_Filter_RANGE[0]) and
                        (Flame_Filter_SPLICECOMB[1] in Flame_Filter_RANGE[1]) and
                        (Flame_Filter_SPLICECOMB[2] in Flame_Filter_RANGE[2])): #If statement that will be flagged if the Exon matches the reference, otherwise, keep looping until you reach the end of the reference file.
                        Flame_Filter_EXONCOUNTER += 1
                        Flame_Filter_Counter += 1
                        break
                    elif ((Flame_Filter_SPLICECOMB[0] not in Flame_Filter_RANGE[0]) or
                          (Flame_Filter_SPLICECOMB[1] not in Flame_Filter_RANGE[1]) or
                          (Flame_Filter_SPLICECOMB[2] not in Flame_Filter_RANGE[2])):
                        Flame_Filter_FAULTYCOUNTER += 1
                        Flame_Filter_Counter += 1
                        continue
            if Flame_Filter_EXONCOUNTER == Flame_Filter_SINGLE_READ_COUNT: #If statement that will determine if the read contains an exon that is not annotated by the reference.
                FILTERCORRECTOUTPUT.append(Flame_Filter_SINGLE_READ)
            elif Flame_Filter_EXONCOUNTER != Flame_Filter_SINGLE_READ_COUNT:
                FILTERFAULTYOUTPUT.append(Flame_Filter_SINGLE_READ)
    return FILTERCORRECTOUTPUT, FILTERFAULTYOUTPUT #Output Object Type: [List, List]


#Central Translate function that contain the translate function for both the Correct Reads as well as the Faulty Reads.
def TRANSLATEFUNC(REF, TRANSLATERANGESIZE, TRANSLATECORRECTINPUT = 0, TRANSLATEFAULTYINPUT = 0):
    if TRANSLATECORRECTINPUT != []: #If statement that will activate depending if the input (Input[2]) is not empty. 
        print("-----------\tInitiate Translate Function, Corrected\t\t\t-----------")
        #-----------Outside Counters and variables-----------#
        Flame_Translate_Counter1 = 0
        Flame_Translate_SINGLE_READ_STRT_ALL = []
        Flame_Translate_SINGLE_READ_LEN_ALL = []
        Flame_Translate_SINGLE_READ_COUNT = 0
        Flame_Translate_Counter2 = 0
        Flame_Translate_TEMPORARY_STRING = "" #Variable housing the temporary translated read as a string due to the string append function (X += str(Y + "-"))
        Flame_Translate_SPLICESTART = 0
        Flame_Translate_SPLICELEN = 0
        Flame_Translate_SPLICESTOP = 0
        Flame_Translate_SPLICECOMB = []
        Flame_Translate_RANGE = []
        Flame_Translate_CORRECT = []
        #-----------The Function Itself-----------#
        for Flame_Translate_SINGLE_READ in TRANSLATECORRECTINPUT: #Loop: Going through each correct read. This is reliant on that either you have sorted the reads between "Correct" and "Faulty" from the previous function (FILTERFUNC) or that you have a reliable dataset to use.
            #For ease of use (Below): Singling out all the relevant information for the bed12 read:
            Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
            Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
            Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
            Flame_Translate_Counter1 = 0
            Flame_Translate_Counter2 = 0
            Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
            while Flame_Translate_Counter1 < Flame_Translate_SINGLE_READ_COUNT: #Loop: Iterate through every single exon (Flame_Translate_Counter1)
                Flame_Translate_SPLICESTART = int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_Counter1])
                Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_Counter1])
                Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                #For ease of use (Below): Collapsing the characteristics into a single list.
                Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                              Flame_Translate_SPLICELEN,
                                              Flame_Translate_SPLICESTOP]
                Flame_Translate_RANGE = [list(range(REF[Flame_Translate_Counter2][1]-TRANSLATERANGESIZE,
                                                    REF[Flame_Translate_Counter2][1]+TRANSLATERANGESIZE)),
                                         list(range(REF[Flame_Translate_Counter2][2]-TRANSLATERANGESIZE,
                                                    REF[Flame_Translate_Counter2][2]+TRANSLATERANGESIZE)),
                                         list(range(REF[Flame_Translate_Counter2][3]-TRANSLATERANGESIZE,
                                                    REF[Flame_Translate_Counter2][3]+TRANSLATERANGESIZE))] #Variance Function (Meant to account for experimental, sequencing, base-calling, and/or alignment errors.)
                if ((Flame_Translate_SPLICECOMB[0] in Flame_Translate_RANGE[0]) and
                    (Flame_Translate_SPLICECOMB[1] in Flame_Translate_RANGE[1]) and
                    (Flame_Translate_SPLICECOMB[2] in Flame_Translate_RANGE[2])): #If statement that will translate the read and append it to the temporary string to be outputed into a list (Flame_Translate_CORRECT)
                    Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter2][0] + "-") #FIX: Could insert a part that removes the last dash to help downstream but does not matter in the grander scheme.
                    Flame_Translate_Counter1 += 1
                    continue 
                elif ((Flame_Translate_SPLICECOMB[0] not in Flame_Translate_RANGE[0]) or
                      (Flame_Translate_SPLICECOMB[1] not in Flame_Translate_RANGE[1]) or
                      (Flame_Translate_SPLICECOMB[2] not in Flame_Translate_RANGE[2])): #If statement that allows for iteration for every reference.
                    Flame_Translate_Counter2 += 1
                    continue
            Flame_Translate_CORRECT.append(Flame_Translate_TEMPORARY_STRING) #Function to output the translated read into a list (Flame_Translate_CORRECT) and recycle to next read in the central for-loop.
    elif TRANSLATECORRECTINPUT == []: #If statement that will activate depending if the input (Input[2]) is empty. This is to just create any output for the function.
        Flame_Translate_CORRECT = []
        
    if TRANSLATEFAULTYINPUT != []: #If statement that will activate depending if the input (Input[3]) is not empty.
        print("-----------\tInitiate Translate Function, Faulty\t\t\t-----------")
        #-----------Outside Counters and variables-----------#
        Flame_Translate_Counter2 = 0
        Flame_Translate_SINGLE_READ_STRT_ALL = []
        Flame_Translate_SINGLE_READ_LEN_ALL = []
        Flame_Translate_SINGLE_READ_COUNT = 0
        Flame_Translate_Counter3 = 0
        Flame_Translate_Counter4 = 0
        Flame_Translate_TEMPORARY_STRING = "" #Variable housing the temporary translated read as a string due to the string append function (X += str(Y + "-"))
        Flame_Translate_SPLICESTART = 0
        Flame_Translate_SPLICELEN = 0
        Flame_Translate_SPLICESTOP = 0
        Flame_Translate_SPLICECOMB = []
        Flame_Translate_RANGE = []
        Flame_Translate_BOOLEAN = False
        Flame_Translate_FAULTY = []
        #-----------The Function Itself-----------#
        for Flame_Translate_SINGLE_READ in TRANSLATEFAULTYINPUT: #Loop: Going through each faulty read. This is reliant on that either you have sorted the reads between "Correct" and "Faulty" from the previous function (FILTERFUNC) or that you have a reliable dataset to use.
            #For ease of use (Below): Singling out all the relevant information for the bed12 read:
            Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
            Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
            Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
            Flame_Translate_Counter3 = 0
            Flame_Translate_Counter4 = 0
            Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
            for Flame_Translate_COUNT1 in range(Flame_Translate_SINGLE_READ_COUNT): #Loop: Iterate through every single exon (Flame_Translate_COUNT1)
                Flame_Translate_Counter3 = 0
                Flame_Translate_Counter4 = 0
                Flame_Translate_SPLICESTART = int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_COUNT1])
                Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_COUNT1])
                Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                #For ease of use (Below): Collapsing the characteristics into a single list.
                Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                              Flame_Translate_SPLICELEN,
                                              Flame_Translate_SPLICESTOP]
                while Flame_Translate_Counter4 != len(REF): #Loop: Forcing to run through all the Exons specified in the reference unless match in which it breaks the loop in order to optimize the running time.
                    Flame_Translate_RANGE = [list(range(REF[Flame_Translate_Counter3][1]-TRANSLATERANGESIZE,
                                                        REF[Flame_Translate_Counter3][1]+TRANSLATERANGESIZE)),
                                             list(range(REF[Flame_Translate_Counter3][2]-TRANSLATERANGESIZE,
                                                        REF[Flame_Translate_Counter3][2]+TRANSLATERANGESIZE)),
                                             list(range(REF[Flame_Translate_Counter3][3]-TRANSLATERANGESIZE,
                                                        REF[Flame_Translate_Counter3][3]+TRANSLATERANGESIZE))] #Variance Function (Meant to account for experimental, sequencing, base-calling, and/or alingment errors.)
                    if ((Flame_Translate_SPLICECOMB[0] in Flame_Translate_RANGE[0]) and
                        (Flame_Translate_SPLICECOMB[1] in Flame_Translate_RANGE[1]) and
                        (Flame_Translate_SPLICECOMB[2] in Flame_Translate_RANGE[2])): #If statement that will translate the read, if it matches, and append it to the temporary string and then break from the while loop to save computational power.
                        Flame_Translate_Counter2 += 1
                        Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter3][0] + ",")
                        Flame_Translate_Counter3 += 1
                        Flame_Translate_BOOLEAN = False 
                        break
                    elif ((Flame_Translate_SPLICECOMB[0] not in Flame_Translate_RANGE[0]) or
                          (Flame_Translate_SPLICECOMB[1] not in Flame_Translate_RANGE[1]) or
                          (Flame_Translate_SPLICECOMB[2] not in Flame_Translate_RANGE[2])): #If statement that will iterate through the reference and once it has passed through the entire reference, it will have the Boolean variable (Flame_Translate_BOOLEAN) as True for downstream.
                        Flame_Translate_Counter3 += 1
                        Flame_Translate_Counter4 += 1
                        Flame_Translate_BOOLEAN = True
                        continue
                if Flame_Translate_BOOLEAN: #Boolean if statement that allows for the insertion of the unreferenced exon range as a "unknown" exon.
                    Flame_Translate_TEMPORARY_STRING += str(str(Flame_Translate_SPLICECOMB[0]) +
                                                            "-" +
                                                            str(Flame_Translate_SPLICECOMB[2]) +
                                                            ",")
                elif Flame_Translate_BOOLEAN == False:
                    pass
            Flame_Translate_FAULTY.append(Flame_Translate_TEMPORARY_STRING) #Function to output the semi-translated read into a list (Flame_Translate_FAULTY) and recycle to next read in the central for-loop.
    elif TRANSLATEFAULTYINPUT == []: #If statement that will activate depending if the input (Input[2]) is empty. This is to just create any output for the function.
        Flame_Translate_FAULTY = []
    return Flame_Translate_CORRECT, Flame_Translate_FAULTY #Output Object Type: [List, List]


#Central Quantificiation function that will quantify the combinations of exons.
def QUANTIFYFUNC(QUANTIFYINPUT):
    Flame_Quantify_DICT = {}
    for Flame_Quantify_COUNT1 in QUANTIFYINPUT:
        if Flame_Quantify_COUNT1 in Flame_Quantify_DICT.keys():
            Flame_Quantify_DICT[Flame_Quantify_COUNT1] += 1
        if Flame_Quantify_COUNT1 not in Flame_Quantify_DICT.keys():
            Flame_Quantify_DICT[Flame_Quantify_COUNT1] = 1
    return Flame_Quantify_DICT #Output Object Type: Dictionary


#Central Creation of empty Adjecency Matrix function that will translate the connection between exons.
def EMPTYADJMTXFUNC(REF):
    #-----------Outside Counters and variables-----------#
    Flame_Emptyadjmtx_MATRIX = []
    #-----------The Function Itself-----------#
    Flame_Emptyadjmtx_REF1 = (len(REF)+1)
    for Flame_Emptyadjmtx_COUNT1 in range(Flame_Emptyadjmtx_REF1):
        Flame_Emptyadjmtx_MATRIX.append([0]*Flame_Emptyadjmtx_REF1)
    return Flame_Emptyadjmtx_MATRIX #Output Object Type: Nested List.

#Central Adjencency Matrix filling function that will count the connections.
def CORRECTADJMTXFUNC(CORRECTADJMTXINPUT, REF, EMPTYADJMTX):
    Flame_Rightadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Rightadjmtx_EXONNR = 0
    Flame_Rightadjmtx_Counter1 = 0
    Flame_Rightadjmtx_Counter2 = 0
    Flame_Rightadjmtx_SPLICECOMB1 = ""
    Flame_Rightadjmtx_SPLICECOMB2 = ""
    for Flame_Rightadjmtx_COUNT1 in CORRECTADJMTXINPUT:
        Flame_Rightadjmtx_EXONNR = (len(Flame_Rightadjmtx_COUNT1.split("-")[:-1]))
        for Flame_Rightadjmtx_COUNT2 in range(Flame_Rightadjmtx_EXONNR):
            Flame_Rightadjmtx_SPLICECOMB1 = str(Flame_Rightadjmtx_COUNT1.split("-")[Flame_Rightadjmtx_COUNT2])
            Flame_Rightadjmtx_SPLICECOMB2 = str(Flame_Rightadjmtx_COUNT1.split("-")[Flame_Rightadjmtx_COUNT2+1])
            Flame_Rightadjmtx_Counter1 = 0
            Flame_Rightadjmtx_Counter2 = 0
            for Flame_Rightadjmtx_COUNT2 in REF:
                if Flame_Rightadjmtx_SPLICECOMB1 == Flame_Rightadjmtx_COUNT2[0]:
                    break
                elif Flame_Rightadjmtx_SPLICECOMB1 != Flame_Rightadjmtx_COUNT2[0]:
                    Flame_Rightadjmtx_Counter1 += 1
            for Flame_Rightadjmtx_COUNT2 in REF:
                if Flame_Rightadjmtx_SPLICECOMB2 == Flame_Rightadjmtx_COUNT2[0]:
                    break
                elif Flame_Rightadjmtx_SPLICECOMB2 != Flame_Rightadjmtx_COUNT2[0]:
                    Flame_Rightadjmtx_Counter2 += 1
            Flame_Rightadjmtx_EMPTYADJMTX[Flame_Rightadjmtx_Counter2][Flame_Rightadjmtx_Counter1] += 1
    return Flame_Rightadjmtx_EMPTYADJMTX #Output Object Type: Nested List. 

def FAULTYADJMTXFUNC(FAULTYADJMTXINPUT, REF, EMPTYADJMTX, ADJMTXRANGESIZE):
    Flame_Wrongadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Wrongadjmtx_Counter1 = 0
    Flame_Wrongadjmtx_Counter2 = 0
    Flame_Wrongadjmtx_RANGE1 = []
    Flame_Wrongadjmtx_RANGE2 = []
    for Flame_Wrongadjmtx_COUNT1 in FAULTYADJMTXINPUT:
        Flame_Wrongadjmtx_Counter1 = 0
        Flame_Wrongadjmtx_Counter2 = 0
        Flame_Wrongadjmtx_RANGE1 = list(range(int(Flame_Wrongadjmtx_COUNT1.spllit("-")[0])-ADJMTXRANGESIZE,
                                              int(Flame_Wrongadjmtx_COUNT1.spllit("-")[0])+ADJMTXRANGESIZE))
        Flame_Wrongadjmtx_RANGE2 = list(range(int(Flame_Wrongadjmtx_COUNT1.spllit("-")[1])-ADJMTXRANGESIZE,
                                              int(Flame_Wrongadjmtx_COUNT1.spllit("-")[1])+ADJMTXRANGESIZE))
        for Flame_Wrongadjmtx_COUNT2 in REF:
            if Flame_Wrongadjmtx_COUNT2[0] in Flame_Wrongadjmtx_RANGE1:
                break
            elif Flame_Wrongadjmtx_COUNT2[0] not in Flame_Wrongadjmtx_RANGE1:
                Flame_Wrongadjmtx_Counter1 += 1
        for Flame_Wrongadjmtx_COUNT2 in REF:
            if Flame_Wrongadjmtx_COUNT2[0] in Flame_Wrongadjmtx_RANGE2:
                break
            elif Flame_Wrongadjmtx_COUNT2[0] not in Flame_Wrongadjmtx_RANGE2:
                Flame_Wrongadjmtx_Counter2 += 1
        Flame_Wrongadjmtx_EMPTYADJMTX[Flame_Wrongadjmtx_Counter2][Flame_Wrongadjmtx_Counter1] += 1
    return Flame_Wrongadjmtx_EMPTYADJMTX #Output Object Type: Nested List.

def FAULTYSEPERATORFUNC(FAULTYSEPINPUT, REF):
    Flame_Faultysep_NAMESEP = []
    Flame_Faultysep_STRT_ALL = []
    Flame_Faultysep_POTENTIAL_EXON = []
    #Create an reference and database
    for Flame_Faultysep_COUNT1 in REF:
        Flame_Faultysep_NAMESEP.append(Flame_Faultysep_COUNT1[0])
    #Create a list containing only the Non-reference exon ranges:
    for Flame_Faultysep_COUNT2 in FAULTYSEPINPUT:
        Flame_Faultysep_STRT_ALL = Flame_Faultysep_COUNT2.split(",")[:-1]
        for Flame_Faultysep_COUNT3 in Flame_Faultysep_STRT_ALL:
            if Flame_Faultysep_COUNT3 in Flame_Faultysep_NAMESEP:
                pass
            elif Flame_Faultysep_COUNT3 not in Flame_Faultysep_NAMESEP:
                Flame_Faultysep_POTENTIAL_EXON.append(Flame_Faultysep_COUNT3)
    return Flame_Faultysep_POTENTIAL_EXON #Output Object Type: List

def FREQUENCYSITEFUNC(FREQUENCYSITEINPUT, REF, FREQUENCYSITERANGESIZE, FREQUENCYSITEWINDOWSIZE):
    Flame_FrequencySite_GENEEMULATION = []
    Flame_FrequencySite_COUNT1RANGE = []
    Flame_FrequencySite_START = []
    Flame_FrequencySite_END = []
    for Flame_FrequencySite_COUNT1 in range(int(REF[-1][3]) +
                                        FREQUENCYSITERANGESIZE):
        Flame_FrequencySite_GENEEMULATION.append(0) #Make X.append([0]*Y) function instead?
    for Flame_FrequencySite_COUNT1 in FREQUENCYSITEINPUT:
        Flame_FrequencySite_COUNT1RANGE = Flame_FrequencySite_COUNT1.split("-") #This splits the two numbers (Start and End). This makes it so one cannot see the connections.
        if int(Flame_FrequencySite_COUNT1RANGE[1]) < (int(REF[-1][3]) +
                                                  FREQUENCYSITERANGESIZE):
            Flame_FrequencySite_START = list(range(int(Flame_FrequencySite_COUNT1RANGE[0])-FREQUENCYSITEWINDOWSIZE,
                                               int(Flame_FrequencySite_COUNT1RANGE[0])+FREQUENCYSITEWINDOWSIZE))
            Flame_FrequencySite_END = list(range(int(Flame_FrequencySite_COUNT1RANGE[1])-FREQUENCYSITEWINDOWSIZE,
                                             int(Flame_FrequencySite_COUNT1RANGE[1])+FREQUENCYSITEWINDOWSIZE))
            for Flame_FrequencySite_COUNT2 in Flame_FrequencySite_START:
                Flame_FrequencySite_GENEEMULATION[Flame_FrequencySite_COUNT2] += 1
            for Flame_FrequencySite_COUNT2 in Flame_FrequencySite_END:
                Flame_FrequencySite_GENEEMULATION[Flame_FrequencySite_COUNT2] += 1
        else:
            pass
    return Flame_FrequencySite_GENEEMULATION #Output Object Type: List

def FREQUENCYTHRESHFUNC(FREQUENCYTRESHINPUT, PERCENTTHRESH, REF): 
    Flame_FrequencyThresh_THRESHOLD = float()
    Flame_FrequencyThresh_SPLICECANDIDATES = []
    Flame_FrequencyThresh_Counter1 = 0
    for Flame_FrequencyThresh_COUNT1 in FREQUENCYTRESHINPUT:
        Flame_FrequencyThresh_THRESHOLD = float(float(PERCENTTHRESH)*len(REF))
        if Flame_FrequencyThresh_COUNT1 > Flame_FrequencyThresh_THRESHOLD:
            Flame_FrequencyThresh_SPLICECANDIDATES.append([Flame_FrequencyThresh_Counter1,
                                                           Flame_FrequencyThresh_COUNT1,
                                                           round((float(Flame_FrequencyThresh_COUNT1)/len(REF))*100, 2)])
            Flame_FrequencyThresh_Counter1 += 1
        else:
            Flame_FrequencyThresh_Counter1 += 1
    return Flame_FrequencyThresh_SPLICECANDIDATES #Output Object Type: List
            
def SPLICESIGNALFUNC(SPLICESIGNALINPUT, REF, SPLICESIGNALGENESTART):
    Flame_SpliceSignal_SPLICESIGNALINPUT = SPLICESIGNALINPUT
    Flame_SpliceSignal_WINDOW = []
    Flame_SpliceSignal_GU = False
    Flame_SpliceSignal_AG = False
    for Flame_SpliceSignal_COUNT1 in Flame_SpliceSignal_SPLICESIGNALINPUT:
        Flame_SpliceSignal_GU = False
        Flame_SpliceSignal_AG = False
        Flame_SpliceSignal_WINDOW = REF[((int(Flame_SpliceSignal_COUNT1[0]) + SPLICESIGNALGENESTART)-3):
                                        ((int(Flame_SpliceSignal_COUNT1[0]) + SPLICESIGNALGENESTART)+4)]
        if "GT" in Flame_SpliceSignal_WINDOW[3:7]:
            Flame_SpliceSignal_GU = True
        if "AG" in Flame_SpliceSignal_WINDOW[0:4]:
            Flame_SpliceSignal_AG = True
        if Flame_SpliceSignal_GU == False and Flame_SpliceSignal_AG == False:
            Flame_SpliceSignal_COUNT1.append("None")
        elif Flame_SpliceSignal_GU == True and Flame_SpliceSignal_AG == False:
            Flame_SpliceSignal_COUNT1.append("GU")
        elif Flame_SpliceSignal_GU == False and Flame_SpliceSignal_AG == True:
            Flame_SpliceSignal_COUNT1.append("AG")
        elif Flame_SpliceSignal_GU == True and Flame_SpliceSignal_AG == True:
            Flame_SpliceSignal_COUNT1.append("Both")
    return Flame_SpliceSignal_SPLICESIGNALINPUT #Output Object Type: List

def SHORTREADFUNC(SHORTREADINPUT, SHORTREADCANDIDATES, REF, SHORTREADGENESTART):
    Flame_ShortRead_SPLICESITECOUNT = {}
    Flame_ShortRead_CANDIDATES = SHORTREADCANDIDATES
    Flame_ShortRead_START = 0
    Flame_ShortRead_CIGAR = ""
    Flame_ShortRead_MATCH1 = None
    Flame_ShortRead_MATCH2 = None
    Flame_ShortRead_MATCH3 = None
    Flame_ShortRead_MATCH4 = None
    Flame_ShortRead_ITEMS = []
    Flame_ShortRead_TMPSTART1 = 0
    Flame_ShortRead_TMPSTART2 = 0
    Flame_ShortRead_CIGARCOUNT = 0
    Flame_ShortRead_CIGAROPERATOR = 0
    Flame_ShortRead_Counter1 = 0
    for Flame_ShortRead_COUNT1 in SHORTREADINPUT:
        Flame_ShortRead_START = int(str(Flame_ShortRead_COUNT1).split("\t")[3])
        Flame_ShortRead_CIGAR = str(Flame_ShortRead_COUNT1).split("\t")[5]
        if "N" in Flame_ShortRead_CIGAR:
            Fĺame_ShortRead_LENGTH = sum(list(map(int, re.findall(r'\d+',
                                                                  Flame_ShortRead_CIGAR))))
            Flame_ShortRead_COMB = [Flame_ShortRead_START,
                                    Flame_ShortRead_CIGAR,
                                    (Flame_ShortRead_START+Fĺame_ShortRead_LENGTH)]
            if (REF[0][1]+SHORTREADGENESTART) < Flame_ShortRead_COMB[0] < (REF[-1][3]+SHORTREADGENESTART):
                Flame_ShortRead_MATCH1 = None
                Flame_ShortRead_MATCH2 = None
                Flame_ShortRead_MATCH3 = None
                Flame_ShortRead_MATCH4 = None
                Flame_ShortRead_ITEMS = []
                if re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                            Flame_ShortRead_COMB[1],
                            re.I):
                    Flame_ShortRead_MATCH1 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                                      Flame_ShortRead_COMB[1],
                                                      re.I)
                elif re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                              Flame_ShortRead_COMB[1],
                              re.I):
                    Flame_ShortRead_MATCH2 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                                      Flame_ShortRead_COMB[1],
                                                      re.I)
                elif re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                              Flame_ShortRead_COMB[1],
                              re.I):
                    Flame_ShortRead_MATCH3 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                                      Flame_ShortRead_COMB[1],
                                                      re.I)
                elif re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                              Flame_ShortRead_COMB[1],
                              re.I):
                    Flame_ShortRead_MATCH4 = re.match(r"([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)([0-9]+)([A-Z]+)",
                                                      Flame_ShortRead_COMB[1],
                                                      re.I)
                if Flame_ShortRead_MATCH1:
                    Flame_ShortRead_ITEMS = list(Flame_ShortRead_MATCH1.groups())
                elif Flame_ShortRead_MATCH2:
                    Flame_ShortRead_ITEMS = list(Flame_ShortRead_MATCH2.groups())
                elif Flame_ShortRead_MATCH3:
                    Flame_ShortRead_ITEMS = list(Flame_ShortRead_MATCH3.groups())
                elif Flame_ShortRead_MATCH4:
                    Flame_ShortRead_ITEMS = list(Flame_ShortRead_MATCH4.groups())
                else:
                    print("Error, CIGAR-String:", Flame_ShortRead_COMB)
                if Flame_ShortRead_ITEMS[1] == "S":
                    del Flame_ShortRead_ITEMS[:2]
                if Flame_ShortRead_ITEMS[-1] == "S":
                    del Flame_ShortRead_ITEMS[-2:]            
                Flame_ShortRead_TMPSTART1 = Flame_ShortRead_START
                #print(Flame_ShortRead_TMPSTART1)
                for Flame_ShortRead_COUNT1 in range(len(Flame_ShortRead_ITEMS)):
                    Flame_ShortRead_CIGARCOUNT = 0
                    Flame_ShortRead_CIGAROPERATOR = 0
                    if Flame_ShortRead_COUNT1 % 2 != 0 and Flame_ShortRead_COUNT1 > 0:
                        #print("A")
                        Flame_ShortRead_CIGARCOUNT = int(Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1-1])
                        Flame_ShortRead_CIGAROPERATOR = Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1]
                        if Flame_ShortRead_CIGAROPERATOR != "N":
                            print("B", Flame_ShortRead_CIGAROPERATOR)
                            Flame_ShortRead_TMPSTART1 = Flame_ShortRead_TMPSTART1 + Flame_ShortRead_CIGARCOUNT
                            pass
                        elif Flame_ShortRead_CIGAROPERATOR == "N":
                            print("C")
                            Flame_ShortRead_TMPSTART2 = Flame_ShortRead_TMPSTART1-STARTPOSITION
                            Flame_ShortRead_TMPSTART1 = Flame_ShortRead_TMPSTART1 + Flame_ShortRead_CIGARCOUNT
                            if Flame_ShortRead_TMPSTART2 in Flame_ShortRead_SPLICESITECOUNT:
                                print("D")
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] += 1
                            elif Flame_ShortRead_TMPSTART2 not in Flame_ShortRead_SPLICESITECOUNT:
                                print("E")
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] = 1
                            Flame_ShortRead_TMPSTART2 = Flame_ShortRead_TMPSTART2 + Flame_ShortRead_CIGARCOUNT
                            if Flame_ShortRead_TMPSTART2 in Flame_ShortRead_SPLICESITECOUNT:
                                print("F")
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] += 1
                            elif Flame_ShortRead_TMPSTART2 not in Flame_ShortRead_SPLICESITECOUNT:
                                print("G")
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] = 1
                    else:
                        pass
    print("BBBBBB", Flame_ShortRead_SPLICESITECOUNT)

    for Flame_ShortRead_COUNT2 in Flame_ShortRead_CANDIDATES:
        Flame_ShortRead_Counter1 = 0
        for k, v in Flame_ShortRead_SPLICESITECOUNT.items():
            if Flame_ShortRead_COUNT2[0] == k:
                Flame_ShortRead_COUNT2.extend(["Yes", v])
                break
            elif Flame_ShortRead_COUNT2[0] != k:
                Flame_ShortRead_Counter1 += 1
        if Flame_ShortRead_Counter1 == len(Flame_ShortRead_SPLICESITECOUNT):
            Flame_ShortRead_COUNT2.extend(["No", "N/A"])
        
    return Flame_ShortRead_CANDIDATES




#Input command:
parser = argparse.ArgumentParser(description = "FLAME: Full Length Adjecency Matrix Enumeration") #CREATE A DEFAULT?!
parser.add_argument("-I", dest = "INPUT", help = "Input file")
parser.add_argument("-R", dest = "REF", help = "Reference File in Fasta format")
parser.add_argument("--range", dest = "RANGE", help = "Variance Range", default = 20)
parser.add_argument("-G", dest = "GTF", help = "Reference File in GTF format")
parser.add_argument("-B", dest = "SAM", help = "Shortread Sequence")
parser.add_argument("-O", dest = "OUTPUT1", help = "Output Prefix", default = "flame")
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
Flame_SPLICECANDIDATES = []

CREATEREF(GTFFILE,
          NSTDLIST,
          STARTPOSITION) #(GTF Filter, Empty List)

print("-----------\tInitiate Filter Function\t\t\t\t-----------")
CORRECREADS, FAULTYREADS = FILTERFUNC(INPUT,
                                      NSTDLIST,
                                      RANGESIZE) #Input, Reference, Rangesize, CorrectOutput, 

#CORRECREADS = []
#FAULTYREADS = []

TRANSLATECORREC, TRANSLATEFAULTY = TRANSLATEFUNC(NSTDLIST,
                                                 RANGESIZE,
                                                 CORRECREADS,
                                                 FAULTYREADS)
ADJMTX1 = EMPTYADJMTXFUNC(NSTDLIST)

ADJMTX1 = CORRECTADJMTXFUNC(TRANSLATECORREC,
                            NSTDLIST,
                            ADJMTX1)

POTNEWEXON = FAULTYSEPERATORFUNC(TRANSLATEFAULTY,
                                 NSTDLIST)

GENEREFERENCE = FREQUENCYSITEFUNC(POTNEWEXON,
                                  NSTDLIST,
                                  RANGESIZE,
                                  WINDOWSIZE)

Flame_SPLICECANDIDATES = FREQUENCYTHRESHFUNC(GENEREFERENCE,
                                             0.01,
                                             FAULTYREADS)

Flame_SPLICECANDIDATES = SPLICESIGNALFUNC(Flame_SPLICECANDIDATES,
                                          REF,
                                          STARTPOSITION)

Flame_SPLICECANDIDATES = SHORTREADFUNC(SHORTREAD,
                                       Flame_SPLICECANDIDATES,
                                       NSTDLIST,
                                       STARTPOSITION)
#print(Flame_SPLICECANDIDATES)
#print()
QUANTIFY = {}
FQUANTIFY = {}


'''
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

for COUNT15 in Flame_ShortRead_CANDIDATES:
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
#Confirmation using short read sequencing, Part 2 - Crossreferencing with Longread:
print("-----------\tFaulty Adjecency Matrix?\t\t\t\t-----------")

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
if CORRECREADS != []:
    OUTPUT = open("%s.Correct.bed" %args.OUTPUT1, "w+")
    for i in CORRECREADS:
        OUTPUT.write(str(i) + "\n")
    OUTPUT.close()
else:
    pass
        
##Print out Faulty:
if FAULTYREADS != []:
    OUTPUT = open("%s.Faulty.bed" %args.OUTPUT1, "w+")
    for i in FAULTYREADS:
        OUTPUT.write(str(i) + "\n")
    OUTPUT.close()
else:
    pass   

##Print out quantification:
if QUANTIFY != {}:
    OUTPUT = open("%s.Quantification.txt" %args.OUTPUT1, "w+")
    for k, v in QUANTIFY.items():
        OUTPUT.write(str(v) +
                     "\t" +
                     str(k) +
                     "\n")
    OUTPUT.close()
else:
    pass

##Print out FAULTYquantification:
if FQUANTIFY != {}:
    OUTPUT = open("%s.QuantificationF.txt" %args.OUTPUT1, "w+") #Change Name of the Output file. ",Faulty"?
    for k, v in FQUANTIFY.items():
        OUTPUT.write(str(v) +
                     "\t" +
                     str(k) +
                     "\n")
    OUTPUT.close()
else:
    pass


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

##Print out Correct:
if TRANSLATECORREC != []:
    OUTPUT = open("%s.CorrectTrans.txt" %args.OUTPUT1, "w+")
    for i in TRANSLATECORREC:
        OUTPUT.write(str(i) + "\n")
    OUTPUT.close
else:
    pass

##Print out Faulty:
if TRANSLATEFAULTY != []:
    OUTPUT = open("%s.FaultyTrans.txt" %args.OUTPUT1, "w+")
    for i in TRANSLATEFAULTY:
        OUTPUT.write(str(i) + "\n")
    OUTPUT.close
else:
    pass


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
for i in Flame_SPLICECANDIDATES:
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
'''
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
'''
OUTPUT = open("%s.Asdf" %args.OUTPUT1, "w+")
for i in POTNEWEXON:
    OUTPUT.write(str(i) + "\n")
OUTPUT.close()

