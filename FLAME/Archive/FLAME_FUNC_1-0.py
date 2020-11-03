#!/usr/bin/python3
import re #time, sys, os, re
import argparse
import pysam
def CREATEREFFUNC(CREATEREFRINPUT, CREATEREFGENESTART, CREATEREFNAME): #CHANGE IT SO IT DOES NOT MATTER IF IT IS GTF OR GFF?
    #-----------Outside counters and variables-----------#
    Flame_CreateRef_REFERENCE = []
    Flame_CreateRef_BACKLOGREPEATS = []
    Flame_CreateRef_NAMELIST = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
                                "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",
                                "XXI", "XXII", "XXIII", "XXIV", "XXV", "XXVI","XXVII", "XXVIII", "XXIX", "XXX",
                                "XXXI", "XXXII", "XXXIII", "XXXIV", "XXXV", "XXXVI", "XXXVII", "XXXVIII", "XXXIX", "XL",
                                "XLI", "XLII", "XLIII", "XLIV", "XLV", "XLVI", "XLVII", "XLVIII", "XLIX", "L",
                                "LI", "LII","LIII", "LIV", "LV", "LVI", "LVII", "LVIII", "LIX", "LX",
                                "LXI", "LXII", "LXIII", "LXIV", "LXV","LXVI", "LXVII", "LXVIII", "LXIX", "LXX",
                                "LXXI", "LXXII", "LXXIII", "LXXIV", "LXXV", "LXXVI", "LXXVII", "LXXVIII","LXXIX", "LXXX",
                                "LXXXI", "LXXXII", "LXXXIII", "LXXXIV", "LXXXV", "LXXXVI", "LXXXVII", "LXXXVIII", "LXXXIX", "XC",
                                "XCI", "XCII", "XCIII", "XCIV", "XCV", "XCVI", "XCVII", "XCVIII", "XCIX", "C"]
    Flame_CreateRef_NAMECOUNT = 0
    Flame_CreateRef_START = 0
    Flame_CreateRef_STOP = 0
    Flame_CreateRef_LEN = 0
    Flame_CreateRef_COMB = []
    #-----------The Function Itself-----------#
    for Flame_CreateRef_COUNT1 in CREATEREFRINPUT.read().split("\n"):
        if (("\texon\t" in Flame_CreateRef_COUNT1) and
            (str(CREATEREFNAME) in Flame_CreateRef_COUNT1) and
            ((str(CREATEREFNAME) + "-") not in Flame_CreateRef_COUNT1.split(";")[0]) and
            len(Flame_CreateRef_COUNT1) >= 1):
            Flame_CreateRef_START = int(Flame_CreateRef_COUNT1.split("\t")[3]) - CREATEREFGENESTART
            Flame_CreateRef_STOP = int(Flame_CreateRef_COUNT1.split("\t")[4]) - CREATEREFGENESTART
            Flame_CreateRef_LEN = Flame_CreateRef_STOP - Flame_CreateRef_START
            Flame_CreateRef_COMB = [Flame_CreateRef_NAMELIST[Flame_CreateRef_NAMECOUNT],
                                    Flame_CreateRef_START,
                                    Flame_CreateRef_LEN,
                                    Flame_CreateRef_STOP]
            if [Flame_CreateRef_START,
                Flame_CreateRef_LEN,
                Flame_CreateRef_STOP] in Flame_CreateRef_BACKLOGREPEATS:
                pass
            elif [Flame_CreateRef_START,
                  Flame_CreateRef_LEN,
                  Flame_CreateRef_STOP] not in Flame_CreateRef_BACKLOGREPEATS:
                Flame_CreateRef_BACKLOGREPEATS.append([Flame_CreateRef_START,
                                                       Flame_CreateRef_LEN,
                                                       Flame_CreateRef_STOP])
                Flame_CreateRef_NAMECOUNT += 1
                Flame_CreateRef_REFERENCE.append(Flame_CreateRef_COMB)
    return Flame_CreateRef_REFERENCE


#Central Filter function that will sort the reads into two different lists (FILTERCORRECTOUTPUT, FILTERFAULTYOUTPUT) depending if they have only contain perfectly matching exons to the reference or if they have a "faulty" exon.
def FILTERFUNC(FILTERINPUT, REF, FILTERRANGESIZE):
    #-----------Outside counters and variables-----------#
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
        print("-----------\tInitiate Translate Function, Corrected\t\t\t\t-----------")
        #-----------Outside counters and variables-----------#
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
                    Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter2][0] + "-") #FIX: Could insert a part that removes the last dash to help downstream but does not matter in the grander scheme. Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element.
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
        print("-----------\tInitiate Translate Function, Faulty\t\t\t\t-----------")
        #-----------Outside counters and variables-----------#
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
                        Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter3][0] + ",") #Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element. list1 = ['1','2','3','4'], s = "-", s = s.join(list1) 
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
                                                            ",") #Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element.
                elif Flame_Translate_BOOLEAN == False:
                    pass
            Flame_Translate_FAULTY.append(Flame_Translate_TEMPORARY_STRING) #Function to output the semi-translated read into a list (Flame_Translate_FAULTY) and recycle to next read in the central for-loop.
    elif TRANSLATEFAULTYINPUT == []: #If statement that will activate depending if the input (Input[2]) is empty. This is to just create any output for the function.
        Flame_Translate_FAULTY = []
    return Flame_Translate_CORRECT, Flame_Translate_FAULTY #Output Object Type: [List, List]


#Central Quantificiation function that will quantify the combinations of exons.
def QUANTIFYFUNC(QUANTIFYINPUT):
    #-----------Outside counters and variables-----------#
    Flame_Quantify_DICT = {}
    for Flame_Quantify_COUNT1 in QUANTIFYINPUT:
        if Flame_Quantify_COUNT1 in Flame_Quantify_DICT.keys():
            Flame_Quantify_DICT[Flame_Quantify_COUNT1] += 1
        if Flame_Quantify_COUNT1 not in Flame_Quantify_DICT.keys():
            Flame_Quantify_DICT[Flame_Quantify_COUNT1] = 1
    return Flame_Quantify_DICT #Output Object Type: Dictionary


#Central Creation of empty Adjecency Matrix function that will translate the connection between exons.
def EMPTYADJMTXFUNC(REF):
    #-----------Outside counters and variables-----------#
    Flame_Emptyadjmtx_MATRIX = []
    #-----------The Function Itself-----------#
    Flame_Emptyadjmtx_REF1 = (len(REF)+1)
    for Flame_Emptyadjmtx_COUNT1 in range(Flame_Emptyadjmtx_REF1):
        Flame_Emptyadjmtx_MATRIX.append([0]*Flame_Emptyadjmtx_REF1)
    return Flame_Emptyadjmtx_MATRIX #Output Object Type: Nested List.

#ADD COMMENTS!
#Central Adjencency Matrix filling function that will count the connections.
def CORRECTADJMTXFUNC(CORRECTADJMTXINPUT, REF, EMPTYADJMTX):
    #-----------Outside counters and variables-----------#
    Flame_Rightadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Rightadjmtx_EXONNR = 0
    Flame_Rightadjmtx_Counter1 = 0
    Flame_Rightadjmtx_Counter2 = 0
    Flame_Rightadjmtx_SPLICECOMB1 = ""
    Flame_Rightadjmtx_SPLICECOMB2 = ""
    #-----------The Function Itself-----------#
    for Flame_Rightadjmtx_COUNT1 in CORRECTADJMTXINPUT:
        Flame_Rightadjmtx_EXONNR = (len(Flame_Rightadjmtx_COUNT1.split("-")[:-1])) #If one changes the "TMPSTRING" function, then one can avoid the use of this "Remove last element function"
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

#ADD COMMENTS!
def FAULTYADJMTXFUNC(FAULTYADJMTXINPUT, REF, EMPTYADJMTX, ADJMTXRANGESIZE):
    #-----------Outside counters and variables-----------#
    Flame_Wrongadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Wrongadjmtx_Counter1 = 0
    Flame_Wrongadjmtx_Counter2 = 0
    Flame_Wrongadjmtx_RANGE1 = []
    Flame_Wrongadjmtx_RANGE2 = []
    #-----------The Function Itself-----------#
    for Flame_Wrongadjmtx_COUNT1 in FAULTYADJMTXINPUT:
        Flame_Wrongadjmtx_Counter1 = 0
        Flame_Wrongadjmtx_Counter2 = 0
        Flame_Wrongadjmtx_RANGE1 = list(range(int(Flame_Wrongadjmtx_COUNT1.split("-")[0])-ADJMTXRANGESIZE,
                                              int(Flame_Wrongadjmtx_COUNT1.split("-")[0])+ADJMTXRANGESIZE))
        Flame_Wrongadjmtx_RANGE2 = list(range(int(Flame_Wrongadjmtx_COUNT1.split("-")[1])-ADJMTXRANGESIZE,
                                              int(Flame_Wrongadjmtx_COUNT1.split("-")[1])+ADJMTXRANGESIZE))
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

#ADD COMMENTS!
def FAULTYSEPERATORFUNC(FAULTYSEPINPUT, REF):
    #-----------Outside counters and variables-----------#
    Flame_Faultysep_NAMESEP = []
    Flame_Faultysep_STRT_ALL = []
    Flame_Faultysep_POTENTIAL_EXON = []
    #-----------The Function Itself-----------#
    #Create an reference and database
    for Flame_Faultysep_COUNT1 in REF:
        Flame_Faultysep_NAMESEP.append(Flame_Faultysep_COUNT1[0])
    #Create a list containing only the Non-reference exon ranges:
    for Flame_Faultysep_COUNT2 in FAULTYSEPINPUT:
        Flame_Faultysep_STRT_ALL = Flame_Faultysep_COUNT2.split(",")[:-1] #If one changes the "TMPSTRING" function, then one can avoid the use of this Remove last element function"
        for Flame_Faultysep_COUNT3 in Flame_Faultysep_STRT_ALL:
            if Flame_Faultysep_COUNT3 in Flame_Faultysep_NAMESEP:
                pass
            elif Flame_Faultysep_COUNT3 not in Flame_Faultysep_NAMESEP:
                Flame_Faultysep_POTENTIAL_EXON.append(Flame_Faultysep_COUNT3)
    return Flame_Faultysep_POTENTIAL_EXON #Output Object Type: List

#ADD COMMENTS!
def FREQUENCYSITEFUNC(FREQUENCYSITEINPUT, REF, FREQUENCYSITERANGESIZE, FREQUENCYSITEWINDOWSIZE):
    #-----------Outside counters and variables-----------#
    Flame_FrequencySite_GENEEMULATION = []
    Flame_FrequencySite_COUNT1RANGE = []
    Flame_FrequencySite_START = []
    Flame_FrequencySite_END = []
    #-----------The Function Itself-----------#
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

#ADD COMMENTS!
def FREQUENCYTHRESHFUNC(FREQUENCYTRESHINPUT, PERCENTTHRESH, REF):
    #-----------Outside counters and variables-----------#
    Flame_FrequencyThresh_THRESHOLD = float()
    Flame_FrequencyThresh_SPLICECANDIDATES = []
    Flame_FrequencyThresh_Counter1 = 0
    #-----------The Function Itself-----------#
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

#ADD COMMENTS!
def SPLICESIGNALFUNC(SPLICESIGNALINPUT, REF, SPLICESIGNALGENESTART):
    #-----------Outside counters and variables-----------#
    Flame_SpliceSignal_SPLICESIGNALINPUT = SPLICESIGNALINPUT
    Flame_SpliceSignal_WINDOW = []
    Flame_SpliceSignal_GU = False
    Flame_SpliceSignal_AG = False
    #-----------The Function Itself-----------#
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

#ADD COMMENTS!
def SHORTREADFUNC(SHORTREADINPUT, SHORTREADCANDIDATES, REF, SHORTREADGENESTART):
    #-----------Outside counters and variables-----------#
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
    #-----------The Function Itself-----------#
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
                for Flame_ShortRead_COUNT1 in range(len(Flame_ShortRead_ITEMS)):
                    Flame_ShortRead_CIGARCOUNT = 0
                    Flame_ShortRead_CIGAROPERATOR = 0
                    if Flame_ShortRead_COUNT1 % 2 != 0 and Flame_ShortRead_COUNT1 > 0:
                        Flame_ShortRead_CIGARCOUNT = int(Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1-1])
                        Flame_ShortRead_CIGAROPERATOR = Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1]
                        if Flame_ShortRead_CIGAROPERATOR != "N":
                            Flame_ShortRead_TMPSTART1 = Flame_ShortRead_TMPSTART1 + Flame_ShortRead_CIGARCOUNT
                            pass
                        elif Flame_ShortRead_CIGAROPERATOR == "N":
                            Flame_ShortRead_TMPSTART2 = Flame_ShortRead_TMPSTART1-SHORTREADGENESTART
                            Flame_ShortRead_TMPSTART1 = Flame_ShortRead_TMPSTART1 + Flame_ShortRead_CIGARCOUNT
                            if Flame_ShortRead_TMPSTART2 in Flame_ShortRead_SPLICESITECOUNT:
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] += 1
                            elif Flame_ShortRead_TMPSTART2 not in Flame_ShortRead_SPLICESITECOUNT:
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] = 1
                            Flame_ShortRead_TMPSTART2 = Flame_ShortRead_TMPSTART2 + Flame_ShortRead_CIGARCOUNT
                            if Flame_ShortRead_TMPSTART2 in Flame_ShortRead_SPLICESITECOUNT:
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] += 1
                            elif Flame_ShortRead_TMPSTART2 not in Flame_ShortRead_SPLICESITECOUNT:
                                Flame_ShortRead_SPLICESITECOUNT[Flame_ShortRead_TMPSTART2] = 1
                    else:
                        pass
    #ADD COMMENTS!
    for Flame_ShortRead_COUNT2 in Flame_ShortRead_CANDIDATES:
        #-----------Outside counters and variables-----------#
        Flame_ShortRead_Counter1 = 0
        #-----------The Function Itself-----------#
        for k, v in Flame_ShortRead_SPLICESITECOUNT.items():
            if Flame_ShortRead_COUNT2[0] == k:
                Flame_ShortRead_COUNT2.extend(["Yes", v])
                break
            elif Flame_ShortRead_COUNT2[0] != k:
                Flame_ShortRead_Counter1 += 1
        if Flame_ShortRead_Counter1 == len(Flame_ShortRead_SPLICESITECOUNT):
            Flame_ShortRead_COUNT2.extend(["No", "N/A"])
        
    return Flame_ShortRead_CANDIDATES



#def main() #FIX: Create a Main function if you have everything ready and just want to run everything in one go.
#if __name__ == "__main__":
# main()