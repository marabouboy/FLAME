#!/usr/bin/python3
#Version 0.1.5
import re, time, sys
import argparse
import pysam

#Central function to produce the Progress bar that is viewed in the terminal. Should be self-contained.
def PROGRESSBAR(PROGRESS):
    Flame_Progress_BARLENGTH = 50 #Change this to modify the size of the progressbar
    Flame_Progress_STATUS = ""
    Flame_Progress_BLOCK = 0
    if isinstance(PROGRESS, int):
        PROGRESS = float(PROGRESS)
    if not isinstance(PROGRESS, float):
        PROGRESS = 0
        Flame_Progress_STATUS = "Error: Progress Variable must be a Float\r\n"
    if PROGRESS < 0:
        PROGRESS = 0
        Flame_Progress_STATUS = "Halt...\r\n"
    if PROGRESS >= 1:
        PROGRESS = 1
        Flame_Progress_STATUS = "Done...\r\n"
    Flame_Progress_BLOCK = int(round(Flame_Progress_BARLENGTH*PROGRESS))
    #Central function to format the progress bar in the terminal
    Flame_Progress_TEXT = "\rStatus: [{0}] {1}% {2}".format( "#"*Flame_Progress_BLOCK +
                                                              "-"*(Flame_Progress_BARLENGTH-Flame_Progress_BLOCK),
                                                              int(round(PROGRESS,2)*100),
                                                              Flame_Progress_STATUS)
    sys.stdout.write(Flame_Progress_TEXT)
    sys.stdout.flush()


#Central function to create an local-variable containing all of the exons in a nested-list-format.
def CREATEREFFUNC(CREATEREFRINPUT, CREATEREFNAME): #CHANGE IT SO IT DOES NOT MATTER IF IT IS GTF OR GFF?
    #-----------Outside counters and variables-----------#
    Flame_CreateRef_REFERENCE = []
    Flame_CreateRef_BACKLOGREPEATS = []
    #The function to name all of the exons. Can be manually changed into whatever fits the user.
    # x equals value if condition else othervalue
    Flame_CreateRef_NAMELIST = [str(x) for x in range(1, 801, 1)]
    Flame_CreateRef_NAMECOUNT = 0
    Flame_CreateRef_START = 0
    Flame_CreateRef_STOP = 0
    Flame_CreateRef_LEN = 0
    Flame_CreateRef_COMB = []
    #-----------The Function Itself-----------#
    for Flame_CreateRef_COUNT1 in CREATEREFRINPUT:
        if (("\texon\t" in Flame_CreateRef_COUNT1) and
            (str("\"" + CREATEREFNAME + "\"") in Flame_CreateRef_COUNT1) and
            ((str(CREATEREFNAME) + "-") not in Flame_CreateRef_COUNT1.split(";")[0]) and
            ((str(CREATEREFNAME) + "-") not in Flame_CreateRef_COUNT1.split(";")[1]) and
            len(Flame_CreateRef_COUNT1) >= 1):
            Flame_CreateRef_CHROM = Flame_CreateRef_COUNT1.split("\t")[0]
            Flame_CreateRef_START = int(Flame_CreateRef_COUNT1.split("\t")[3]) 
            Flame_CreateRef_STOP = int(Flame_CreateRef_COUNT1.split("\t")[4]) 
            Flame_CreateRef_LEN = Flame_CreateRef_STOP - Flame_CreateRef_START + 1
            Flame_CreateRef_COMB = ["",
                                    Flame_CreateRef_CHROM,
                                    Flame_CreateRef_START,
                                    Flame_CreateRef_LEN,
                                    Flame_CreateRef_STOP]
            if Flame_CreateRef_COMB in Flame_CreateRef_BACKLOGREPEATS:
                pass
            elif Flame_CreateRef_COMB not in Flame_CreateRef_BACKLOGREPEATS:
                Flame_CreateRef_BACKLOGREPEATS.append(Flame_CreateRef_COMB)
    Flame_CreateRef_REFERENCE = sorted(Flame_CreateRef_BACKLOGREPEATS,
                                       key = lambda x: (x[2], x[4]))
    for Flame_CreateRef_COUNT1 in Flame_CreateRef_REFERENCE:
        Flame_CreateRef_COUNT1[0] = Flame_CreateRef_NAMELIST[Flame_CreateRef_NAMECOUNT]
        Flame_CreateRef_NAMECOUNT += 1
    return Flame_CreateRef_REFERENCE

#Central function to create an local-variable containing all of the exons in a nested-list-format.
def SEGMENTFUNC(SEGMENTINPUT1, REF):
    #-----------Outside counters and variables-----------#
    Flame_Segment_PROGRESSMAX = len(SEGMENTINPUT1)
    Flame_Segment_PROGRESSCOUNT1 = 0
    Flame_Segment_PROGRESSCOUNT2 = 0
    Flame_Segment_Counter1 = 0
    Flame_Segment_Counter2 = 0
    Flame_Segment_OUTPUT = []
    #-----------The Function Itself-----------#
    if len(set(list(Chrom[1] for Chrom in REF))) == 1:
        while Flame_Segment_Counter1 < len(SEGMENTINPUT1):
            if str(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[0]) != REF[0][1]:
                #The reason for the +100 is for a faster cycling of all reads. If one wants to be more accurate, then one needs to lessen this number, opposite being also possible.
                Flame_Segment_Counter1 += 100
                Flame_Segment_PROGRESSCOUNT1 += 100
                Flame_Segment_PROGRESSCOUNT2 += 100
            else:
                if ((str(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[0]) == str(REF[Flame_Segment_Counter2][1])) and
                    (int(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[1]) >= min(list(Start[2] for Start in REF))) and
                    (int(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[2]) <= max(list(End[4] for End in REF)))):
                    Flame_Segment_OUTPUT.append(SEGMENTINPUT1[Flame_Segment_Counter1])
                    Flame_Segment_Counter1 += 1
                    Flame_Segment_PROGRESSCOUNT1 += 1
                    Flame_Segment_PROGRESSCOUNT2 += 1
                else:
                    Flame_Segment_Counter1 += 1
                    Flame_Segment_PROGRESSCOUNT1 += 1
                    Flame_Segment_PROGRESSCOUNT2 += 1
    else:
        while Flame_Segment_Counter1 < len(SEGMENTINPUT1):
            for Flame_Segment_COUNT in REF:
                if ((str(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[0]) == str(Flame_Segment_COUNT[1])) and
                    (int(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[1]) >= int(Flame_Segment_COUNT[2])) and
                    (int(SEGMENTINPUT1[Flame_Segment_Counter1].split("\t")[2]) <= int(Flame_Segment_COUNT[4]))):
                    Flame_Segment_OUTPUT.append(SEGMENTINPUT1[Flame_Segment_Counter1])
                else:
                    pass
            Flame_Segment_Counter1 += 1
            Flame_Segment_PROGRESSCOUNT1 += 1
            Flame_Segment_PROGRESSCOUNT2 += 1
    #-----------The Progress Bar Start-----------#
    if Flame_Segment_PROGRESSCOUNT2 != Flame_Segment_PROGRESSMAX:
        if Flame_Segment_PROGRESSCOUNT1 >= int(round(Flame_Segment_PROGRESSMAX*0.01, 2)):
            PROGRESSBAR(Flame_Segment_PROGRESSCOUNT2 / Flame_Segment_PROGRESSMAX)
            Flame_Segment_PROGRESSCOUNT1 = 0
            time.sleep(0.001)
    elif Flame_Segment_PROGRESSCOUNT2 == Flame_Segment_PROGRESSMAX:
        PROGRESSBAR(Flame_Segment_PROGRESSCOUNT2/Flame_Segment_PROGRESSMAX)
        time.sleep(0.001)
    #-----------The Progress Bar End-----------#    
    return(Flame_Segment_OUTPUT)


#Central Filter function that will sort the reads into two different lists (FILTERANNOTATEDOUTPUT, FILTERINCONGRUENTOUTPUT) depending if they have only contain perfectly matching exons to the reference or if they have a "faulty" exon.
def FILTERFUNC(FILTERINPUT1, REF, FILTERRANGESIZE):
    #-----------Outside counters and variables-----------#
    Flame_Filter_PROGRESSMAX = sum(1 for Flame_Filter_LINE in FILTERINPUT1)
    Flame_Filter_PROGRESSCOUNT1 = 0
    Flame_Filter_PROGRESSCOUNT2 = 0
    Flame_Filter_EXONCOUNTER = 0 
    Flame_Filter_Counter = 0
    Flame_Filter_INCONGRUENTCOUNTER = 0
    FILTERANNOTATEDOUTPUT = []
    FILTERINCONGRUENTOUTPUT = []
    Flame_Filter_SINGLE_READ_STRT_ALL = []
    Flame_Filter_SINGLE_READ_LEN_ALL = []
    Flame_Filter_SINGLE_READ_COUNT = 0
    Flame_Filter_SPLICESTART = 0
    Flame_Filter_SPLICELEN = 0
    Flame_Filter_SPLICESTOP = 0
    Flame_Filter_SPLICECOMB = []
    #-----------The Function Itself-----------#
    #Going through each longread sequencing read.
    for Flame_Filter_SINGLE_READ in FILTERINPUT1:  #FIX: Change it so that the input read (As an outside read file) is read outside?
        #If statement that filters out empty reads.
        if len(Flame_Filter_SINGLE_READ) >= 1: 
            #For ease of use: Singling out all the relevant information from the bed12 read:
            Flame_Filter_SINGLE_READ_STRTPNT = int(Flame_Filter_SINGLE_READ.split("\t")[1])
            Flame_Filter_SINGLE_READ_STRT_ALL = Flame_Filter_SINGLE_READ.split("\t")[11]
            Flame_Filter_SINGLE_READ_LEN_ALL =  Flame_Filter_SINGLE_READ.split("\t")[10]
            Flame_Filter_SINGLE_READ_COUNT = int(Flame_Filter_SINGLE_READ.split("\t")[9])
            Flame_Filter_EXONCOUNTER = 0
            #Going through each exon.
            for Flame_Filter_COUNT1 in range(Flame_Filter_SINGLE_READ_COUNT): 
                Flame_Filter_Counter = 0
                Flame_Filter_INCONGRUENTCOUNTER = 0
                #For ease of use: Specifying each characteristic (Start, Length, Stop) for each exon in each read.
                Flame_Filter_SPLICESTART = (int(Flame_Filter_SINGLE_READ_STRT_ALL.split(",")[Flame_Filter_COUNT1]) +
                                            Flame_Filter_SINGLE_READ_STRTPNT)
                Flame_Filter_SPLICELEN = int(Flame_Filter_SINGLE_READ_LEN_ALL.split(",")[Flame_Filter_COUNT1])
                Flame_Filter_SPLICESTOP = (Flame_Filter_SPLICESTART + Flame_Filter_SPLICELEN)
                #For ease of use: Collapsing the characteristics into a single list.
                Flame_Filter_SPLICECOMB = [Flame_Filter_SPLICESTART,
                                           Flame_Filter_SPLICELEN,
                                           Flame_Filter_SPLICESTOP]
                #Going through each reference exon trying to match with the exon.
                while Flame_Filter_INCONGRUENTCOUNTER != len(REF):
                    #If statement that will be flagged if the Exon matches the reference.
                    if ((REF[Flame_Filter_Counter][2]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[0] <= REF[Flame_Filter_Counter][2]+FILTERRANGESIZE) and
                        (REF[Flame_Filter_Counter][3]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[1] <= REF[Flame_Filter_Counter][3]+FILTERRANGESIZE) and
                        (REF[Flame_Filter_Counter][4]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[2] <= REF[Flame_Filter_Counter][4]+FILTERRANGESIZE)):
                        Flame_Filter_EXONCOUNTER += 1
                        Flame_Filter_Counter += 1                        
                        break
                    #If otherwise, keep looping until you reach the end of the reference file.
                    elif ((not (REF[Flame_Filter_Counter][2]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[0] <= REF[Flame_Filter_Counter][2]+FILTERRANGESIZE)) or
                          (not (REF[Flame_Filter_Counter][3]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[1] <= REF[Flame_Filter_Counter][3]+FILTERRANGESIZE)) or
                          (not (REF[Flame_Filter_Counter][4]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[2] <= REF[Flame_Filter_Counter][4]+FILTERRANGESIZE))):
                        Flame_Filter_INCONGRUENTCOUNTER += 1
                        Flame_Filter_Counter += 1
                        continue
            #If statement that will determine if the read contains an exon that is not annotated by the reference.
            if Flame_Filter_EXONCOUNTER == Flame_Filter_SINGLE_READ_COUNT: 
                FILTERANNOTATEDOUTPUT.append(Flame_Filter_SINGLE_READ)
            elif Flame_Filter_EXONCOUNTER != Flame_Filter_SINGLE_READ_COUNT:
                FILTERINCONGRUENTOUTPUT.append(Flame_Filter_SINGLE_READ)
        #-----------The Progress Bar Start-----------#
        Flame_Filter_PROGRESSCOUNT1 += 1
        Flame_Filter_PROGRESSCOUNT2 += 1
        if Flame_Filter_PROGRESSCOUNT2 != Flame_Filter_PROGRESSMAX:
            if Flame_Filter_PROGRESSCOUNT1 >= int(round(Flame_Filter_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_Filter_PROGRESSCOUNT2 / Flame_Filter_PROGRESSMAX)
                Flame_Filter_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_Filter_PROGRESSCOUNT2 == Flame_Filter_PROGRESSMAX:
            PROGRESSBAR(Flame_Filter_PROGRESSCOUNT2/Flame_Filter_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return FILTERANNOTATEDOUTPUT, FILTERINCONGRUENTOUTPUT #Output Object Type: [List, List]


#Central Translate function that contain the translate function for both the Annotated Reads as well as the Incongruent Reads.
def TRANSLATEFUNC(REF, TRANSLATERANGESIZE, TRANSLATEANNOTATEDINPUT = [], TRANSLATEINCONGRUENTINPUT = []):
    #If statement that will activate depending if the input (Input[2]) is not empty. Linked with Line 230.
    if TRANSLATEANNOTATEDINPUT != []: 
        print("-----------\tInitiate Translate Function, Annotated\t\t\t\t-----------")
        #-----------Outside counters and variables-----------#
        Flame_Translate_PROGRESSMAX = sum(1 for Flame_Translate_LINE in TRANSLATEANNOTATEDINPUT)
        Flame_Translate_PROGRESSCOUNT1 = 0
        Flame_Translate_PROGRESSCOUNT2 = 0
        Flame_Translate_Counter1 = 0
        Flame_Translate_SINGLE_READ_STRT_ALL = []
        Flame_Translate_SINGLE_READ_LEN_ALL = []
        Flame_Translate_SINGLE_READ_COUNT = 0
        Flame_Translate_RefCounter = 0
        Flame_Translate_Counter2 = 0
        Flame_Translate_TEMPORARY_STRING = "" #Variable housing the temporary translated read as a string due to the string append function (X += str(Y + "-"))
        Flame_Translate_SPLICESTART = 0
        Flame_Translate_SPLICELEN = 0
        Flame_Translate_SPLICESTOP = 0
        Flame_Translate_REFSTART = []
        Flame_Translate_SPLICECOMB = []
        Flame_Translate_ANNOTATED = []
        #-----------The Function Itself-----------#
        #Going through each correct read. This is reliant on that either you have sorted the reads between "Annotated" and "Incongruent" from the previous function (FILTERFUNC) or that you have a reliable dataset to use.
        for Flame_Translate_SINGLE_READ in TRANSLATEANNOTATEDINPUT:
            #For ease of use: Singling out all the relevant information for the bed12 read:
            Flame_Translate_STRTPNT = int(Flame_Translate_SINGLE_READ.split("\t")[1])
            Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
            Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
            Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
            Flame_Translate_REFSTARTS = [[item[2]-100 for item in REF], [item[2]+100 for item in REF]]
            Flame_Translate_RefCounter = 0
            #print(Flame_Translate_REFSTARTS)
            Flame_Translate_Counter1 = 0
            Flame_Translate_Counter2 = 0
            Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
            #Iterate through every single exon and due to TRANSLATEANNOTATEDINPUT read count for each exon being known, a while-loop is implemented.
            while Flame_Translate_Counter1 < Flame_Translate_SINGLE_READ_COUNT:
                Flame_Translate_SPLICESTART = (int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_Counter1]) +
                                               Flame_Translate_STRTPNT)
                Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_Counter1])
                Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                #For ease of use: Collapsing the characteristics into a single list.
                Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                              Flame_Translate_SPLICELEN,
                                              Flame_Translate_SPLICESTOP]                            
                #If statement that will be flagged if the Exon matches the reference +- variance-variable (TRANSLATERANGESIZE). Linked with line 210.
                if ((REF[Flame_Translate_Counter2][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter2][2]+TRANSLATERANGESIZE) and
                    (REF[Flame_Translate_Counter2][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter2][3]+TRANSLATERANGESIZE) and
                    (REF[Flame_Translate_Counter2][4]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter2][4]+TRANSLATERANGESIZE)):
                    Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter2][0] + ",") #FIX: Could insert a part that removes the last element to help downstream but does not matter in the grander scheme. Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element.
                    Flame_Translate_Counter1 += 1
                    Flame_Translate_Counter2 -= 2 #Temporary Band-aid for when the Variable-function window is too large and covers other entire exons. Fix this with a more sophisticated solution!
                    continue
                #If otherwise, keep looping until you reach the end of the reference file. Linked with line 203.
                elif ((not (REF[Flame_Translate_Counter2][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter2][2]+TRANSLATERANGESIZE)) or
                      (not (REF[Flame_Translate_Counter2][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter2][3]+TRANSLATERANGESIZE)) or
                      (not (REF[Flame_Translate_Counter2][4]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter2][4]+TRANSLATERANGESIZE))):
                    Flame_Translate_Counter2 += 1
            Flame_Translate_ANNOTATED.append(Flame_Translate_TEMPORARY_STRING)
            #-----------The Progress Bar Start-----------#
            Flame_Translate_PROGRESSCOUNT1 += 1
            Flame_Translate_PROGRESSCOUNT2 += 1
            if Flame_Translate_PROGRESSCOUNT2 != Flame_Translate_PROGRESSMAX:
                if Flame_Translate_PROGRESSCOUNT1 >= int(round(Flame_Translate_PROGRESSMAX*0.01, 2)):
                    PROGRESSBAR(Flame_Translate_PROGRESSCOUNT2 / Flame_Translate_PROGRESSMAX)
                    Flame_Translate_PROGRESSCOUNT1 = 0
                    time.sleep(0.001)
            elif Flame_Translate_PROGRESSCOUNT2 == Flame_Translate_PROGRESSMAX:
                PROGRESSBAR(Flame_Translate_PROGRESSCOUNT2/Flame_Translate_PROGRESSMAX)
                time.sleep(0.001)
            #-----------The Progress Bar End-----------#
    #If statement that will activate depending if the input (Input[2]) is empty. Linked with Line 164.
    else: 
        #print("-----------\tInitiate Translate Function, Annotated\t\t\t\t-----------")
        #print("Status: [] Error... Empty Input")
        Flame_Translate_ANNOTATED = []
        pass
        
    if TRANSLATEINCONGRUENTINPUT != []: #If statement that will activate depending if the input (Input[3]) is not empty. Linked with Line 312.
        print("-----------\tInitiate Translate Function, Incongruent\t\t\t-----------")
        #-----------Outside counters and variables-----------#
        Flame_Translate_PROGRESSMAX = sum(1 for Flame_Translate_LINE in TRANSLATEINCONGRUENTINPUT)
        Flame_Translate_PROGRESSCOUNT1 = 0
        Flame_Translate_PROGRESSCOUNT2 = 0
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
        Flame_Translate_BOOLEAN = False
        Flame_Translate_INCONGRUENT = []
        #-----------The Function Itself-----------#
        #If-statement that checks whether the reference is empty to account for an empty reference scenario:
        if len(REF) > 0:
            #Going through each faulty read. This is reliant on that either you have sorted the reads between "Annotated" and "Incongruent" from the previous function (FILTERFUNC) or that you have a reliable dataset to use.
            for Flame_Translate_SINGLE_READ in TRANSLATEINCONGRUENTINPUT:
                #For ease of use: Singling out all the relevant information for the bed12 read:
                Flame_Translate_STRTPNT = int(Flame_Translate_SINGLE_READ.split("\t")[1])
                Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
                Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
                Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
                Flame_Translate_Counter3 = 0
                Flame_Translate_Counter4 = 0
                Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
                #Iterate through every single exon (Flame_Translate_COUNT1)
                for Flame_Translate_COUNT1 in range(Flame_Translate_SINGLE_READ_COUNT): 
                    Flame_Translate_Counter3 = 0
                    Flame_Translate_Counter4 = 0
                    Flame_Translate_SPLICESTART = (int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_COUNT1]) +
                                                   Flame_Translate_STRTPNT)
                    Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_COUNT1])
                    Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                    #For ease of use: Collapsing the characteristics into a single list.
                    Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                                  Flame_Translate_SPLICELEN,
                                                  Flame_Translate_SPLICESTOP]
                    #Forcing to run through all the Exons specified in the reference unless match in which it breaks the loop in order to optimize the running time.                
                    while Flame_Translate_Counter4 != len(REF):
                        #If statement that will be flagged if the Exon matches the reference +- variance-variable (TRANSLATERANGESIZE). Linked with line 289.
                        if ((REF[Flame_Translate_Counter3][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter3][2]+TRANSLATERANGESIZE) and
                            (REF[Flame_Translate_Counter3][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter3][3]+TRANSLATERANGESIZE) and
                            (REF[Flame_Translate_Counter3][4]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter3][4]+TRANSLATERANGESIZE)):
                            Flame_Translate_Counter2 += 1
                            Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter3][0] + ",") #FIX: Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element. list1 = ['1','2','3','4'], s = "-", s = s.join(list1) 
                            Flame_Translate_Counter3 += 1
                            Flame_Translate_BOOLEAN = False
                            break
                        #If otherwise, keep looping until you reach the end of the reference file. Linked with line 281.
                        elif ((not (REF[Flame_Translate_Counter3][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter3][2]+TRANSLATERANGESIZE)) or
                              (not (REF[Flame_Translate_Counter3][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter3][3]+TRANSLATERANGESIZE)) or
                              (not (REF[Flame_Translate_Counter3][4]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter3][4]+TRANSLATERANGESIZE))):
                            Flame_Translate_Counter3 += 1
                            Flame_Translate_Counter4 += 1
                            Flame_Translate_BOOLEAN = True
                            continue
                    #Boolean if statement that allows for the insertion of the unreferenced exon range as a "unknown" exon.
                    if Flame_Translate_BOOLEAN: 
                        Flame_Translate_TEMPORARY_STRING += str(str(Flame_Translate_SPLICECOMB[0]) +
                                                                "-" +
                                                                str(Flame_Translate_SPLICECOMB[2]) +
                                                                ",") #FIX: Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element.
                    elif Flame_Translate_BOOLEAN == False:
                        pass
                Flame_Translate_INCONGRUENT.append(Flame_Translate_TEMPORARY_STRING) #Function to output the semi-translated read into a list (Flame_Translate_INCONGRUENT) and recycle to next read in the central for-loop.
                #-----------The Progress Bar Start-----------#
                Flame_Translate_PROGRESSCOUNT1 += 1
                Flame_Translate_PROGRESSCOUNT2 += 1
                if Flame_Translate_PROGRESSCOUNT2 != Flame_Translate_PROGRESSMAX:
                    if Flame_Translate_PROGRESSCOUNT1 >= int(round(Flame_Translate_PROGRESSMAX*0.01, 2)):
                        PROGRESSBAR(Flame_Translate_PROGRESSCOUNT2 / Flame_Translate_PROGRESSMAX)
                        Flame_Translate_PROGRESSCOUNT1 = 0
                        time.sleep(0.001)
                elif Flame_Translate_PROGRESSCOUNT2 == Flame_Translate_PROGRESSMAX:
                    PROGRESSBAR(Flame_Translate_PROGRESSCOUNT2/Flame_Translate_PROGRESSMAX)
                    time.sleep(0.001)
                #-----------The Progress Bar End-----------#
        elif len(REF) == 0:
            for Flame_Translate_SINGLE_READ in TRANSLATEINCONGRUENTINPUT:
                #For ease of use: Singling out all the relevant information for the bed12 read:
                Flame_Translate_STRTPNT = int(Flame_Translate_SINGLE_READ.split("\t")[1])
                Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
                Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
                Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
                Flame_Translate_Counter3 = 0
                Flame_Translate_Counter4 = 0
                Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
                #Iterate through every single exon (Flame_Translate_COUNT1)
                for Flame_Translate_COUNT1 in range(Flame_Translate_SINGLE_READ_COUNT): 
                    Flame_Translate_Counter3 = 0
                    Flame_Translate_Counter4 = 0
                    Flame_Translate_SPLICESTART = (int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_COUNT1]) +
                                                   Flame_Translate_STRTPNT)
                    Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_COUNT1])
                    Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                    #For ease of use: Collapsing the characteristics into a single list.
                    Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                                  Flame_Translate_SPLICELEN,
                                                  Flame_Translate_SPLICESTOP]
                    Flame_Translate_TEMPORARY_STRING += str(str(Flame_Translate_SPLICECOMB[0]) +
                                                                "-" +
                                                                str(Flame_Translate_SPLICECOMB[2]) +
                                                                ",") #FIX: Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element.
                    
                Flame_Translate_INCONGRUENT.append(Flame_Translate_TEMPORARY_STRING)
    else: #If statement that will activate depending if the input (Input[2]) is empty. Linked with Line 236
        #print("-----------\tInitiate Translate Function, Incongruent\t\t\t-----------")
        #print("Status: [] Error... Empty Input")
        Flame_Translate_INCONGRUENT = []
        pass
    return Flame_Translate_ANNOTATED, Flame_Translate_INCONGRUENT #Output Object Type: [List, List]

#Central Quantificiation function that will quantify the combinations of exons.
def QUANTIFYFUNC(QUANTIFYINPUT):
    #-----------Outside counters and variables-----------#
    Flame_Quantify_PROGRESSMAX = sum(1 for Flame_Quantify_LINE in QUANTIFYINPUT)
    Flame_Quantify_PROGRESSCOUNT1 = 0
    Flame_Quantify_PROGRESSCOUNT2 = 0
    Flame_Quantify_DICT = {}
    #-----------The Function Itself-----------#
    if QUANTIFYINPUT != []:
        #Quantifying the Splice Permutaitons by going through a dictionary and looking for new keys.
        for Flame_Quantify_COUNT1 in QUANTIFYINPUT:
            #If the are matching keys, just add to it's count.
            if Flame_Quantify_COUNT1 in Flame_Quantify_DICT.keys():
                Flame_Quantify_DICT[Flame_Quantify_COUNT1] += 1
            else:
                #If there is no matching keys, create a new splice permutation combination.
                Flame_Quantify_DICT[Flame_Quantify_COUNT1] = 1
        #else:
            #if Flame_Quantify_COUNT1 not in Flame_Quantify_DICT.keys():
            #    Flame_Quantify_DICT[Flame_Quantify_COUNT1] = 1
            #-----------The Progress Bar Start-----------#
            Flame_Quantify_PROGRESSCOUNT1 += 1
            Flame_Quantify_PROGRESSCOUNT2 += 1
            if Flame_Quantify_PROGRESSCOUNT2 != Flame_Quantify_PROGRESSMAX:
                if Flame_Quantify_PROGRESSCOUNT1 >= int(round(Flame_Quantify_PROGRESSMAX*0.01, 2)):
                    PROGRESSBAR(Flame_Quantify_PROGRESSCOUNT2 / Flame_Quantify_PROGRESSMAX)
                    Flame_Quantify_PROGRESSCOUNT1 = 0
                    time.sleep(0.001)
            elif Flame_Quantify_PROGRESSCOUNT2 == Flame_Quantify_PROGRESSMAX:
                PROGRESSBAR(Flame_Quantify_PROGRESSCOUNT2/Flame_Quantify_PROGRESSMAX)
                time.sleep(0.001)
            #-----------The Progress Bar End-----------#
    elif QUANTIFYINPUT == []:
        print("Status: [] Error... Empty Input")
    return Flame_Quantify_DICT #Output Object Type: Dictionary


#Central Creation of empty Adjacency Matrix function that will translate the connection between exons.
def EMPTYADJMTXFUNC(REF):
    #-----------Outside counters and variables-----------#
    Flame_Emptyadjmtx_PROGRESSMAX = (sum(1 for Flame_Emptyadjmtx_LINE in REF)+1)
    Flame_Emptyadjmtx_PROGRESSCOUNT1 = 0
    Flame_Emptyadjmtx_PROGRESSCOUNT2 = 0
    Flame_Emptyadjmtx_MATRIX = []
    #-----------The Function Itself-----------#
    Flame_Emptyadjmtx_REF1 = (len(REF)+1)
    #Create an empty adjacency matrix based on the number of elements as the input. This part is responsible for the correct number of columns.
    for Flame_Emptyadjmtx_COUNT1 in range(Flame_Emptyadjmtx_REF1):
        Flame_Emptyadjmtx_MATRIX.append([0]*Flame_Emptyadjmtx_REF1) #This part is responsible for the correct number of rows.
        #-----------The Progress Bar Start-----------#
        Flame_Emptyadjmtx_PROGRESSCOUNT1 += 1
        Flame_Emptyadjmtx_PROGRESSCOUNT2 += 1
        if Flame_Emptyadjmtx_PROGRESSCOUNT2 != Flame_Emptyadjmtx_PROGRESSMAX:
            if Flame_Emptyadjmtx_PROGRESSCOUNT1 >= int(round(Flame_Emptyadjmtx_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_Emptyadjmtx_PROGRESSCOUNT2 / Flame_Emptyadjmtx_PROGRESSMAX)
                Flame_Emptyadjmtx_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_Emptyadjmtx_PROGRESSCOUNT2 == Flame_Emptyadjmtx_PROGRESSMAX:
            PROGRESSBAR(Flame_Emptyadjmtx_PROGRESSCOUNT2/Flame_Emptyadjmtx_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_Emptyadjmtx_MATRIX #Output Object Type: Nested List.


#Central Adjancency Matrix filling function that will count the connections within the Annotated Dataset.
def ANNOTATEDADJMTXFUNC(ANNOTATEDADJMTXINPUT, REF, EMPTYADJMTX):
    #-----------Outside counters and variables-----------#
    Flame_Rightadjmtx_PROGRESSMAX = sum(1 for Flame_Rightadjmtx_LINE in ANNOTATEDADJMTXINPUT)
    Flame_Rightadjmtx_PROGRESSCOUNT1 = 0
    Flame_Rightadjmtx_PROGRESSCOUNT2 = 0
    Flame_Rightadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Rightadjmtx_EXONNR = 0
    Flame_Rightadjmtx_Counter1 = 0
    Flame_Rightadjmtx_Counter2 = 0 
    Flame_Rightadjmtx_SPLICECOMB1 = ""
    Flame_Rightadjmtx_SPLICECOMB2 = ""
    #-----------The Function Itself-----------#
    #If statement that will activate depending if the input is empty. Linked with Line 431
    if ANNOTATEDADJMTXINPUT != []:
        #Iterates through the input (Translated Annotated Reads).
        for Flame_Rightadjmtx_COUNT1 in ANNOTATEDADJMTXINPUT: 
            Flame_Rightadjmtx_EXONNR = (len(Flame_Rightadjmtx_COUNT1.split(",")[:-1])) #If one changes the "TMPSTRING" function, then one can avoid the use of this "Remove last element function"
            for Flame_Rightadjmtx_COUNT2 in range(Flame_Rightadjmtx_EXONNR):
                Flame_Rightadjmtx_SPLICECOMB1 = str(Flame_Rightadjmtx_COUNT1.split(",")[Flame_Rightadjmtx_COUNT2])
                Flame_Rightadjmtx_SPLICECOMB2 = str(Flame_Rightadjmtx_COUNT1.split(",")[Flame_Rightadjmtx_COUNT2+1])
                Flame_Rightadjmtx_Counter1 = 0
                Flame_Rightadjmtx_Counter2 = 0
                #Crosschecks with the reference. If it matches, breaks out of the loop, otherwise, add to a counter to know which position in the Adjacency matrix it needs to add to.
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
            #-----------The Progress Bar Start-----------#
            Flame_Rightadjmtx_PROGRESSCOUNT1 += 1
            Flame_Rightadjmtx_PROGRESSCOUNT2 += 1
            if Flame_Rightadjmtx_PROGRESSCOUNT2 != Flame_Rightadjmtx_PROGRESSMAX:
                if Flame_Rightadjmtx_PROGRESSCOUNT1 >= int(round(Flame_Rightadjmtx_PROGRESSMAX*0.01, 2)):
                    PROGRESSBAR(Flame_Rightadjmtx_PROGRESSCOUNT2 / Flame_Rightadjmtx_PROGRESSMAX)
                    Flame_Rightadjmtx_PROGRESSCOUNT1 = 0
                    time.sleep(0.001)
            elif Flame_Rightadjmtx_PROGRESSCOUNT2 == Flame_Rightadjmtx_PROGRESSMAX:
                PROGRESSBAR(Flame_Rightadjmtx_PROGRESSCOUNT2/Flame_Rightadjmtx_PROGRESSMAX)
                time.sleep(0.001)
            #-----------The Progress Bar End-----------#
    #If statement that will activate depending if the input is empty. Linked with Line 431
    elif ANNOTATEDADJMTXINPUT == []:
        print("Status: [] Error... Empty Input")
    return Flame_Rightadjmtx_EMPTYADJMTX #Output Object Type: Nested List. 

#Central Adjancency Matrix filling function that will count the connections within the Incongruent Dataset.
def INCONGRUENTADJMTXFUNC(INCONGRUENTADJMTXINPUT, REF, EMPTYADJMTX, ADJMTXRANGESIZE):
    #-----------Outside counters and variables-----------#
    Flame_Wrongadjmtx_PROGRESSMAX = sum(1 for Flame_Wrongadjmtx_LINE in INCONGRUENTADJMTXINPUT)
    Flame_Wrongadjmtx_PROGRESSCOUNT1 = 0
    Flame_Wrongadjmtx_PROGRESSCOUNT2 = 0
    Flame_Wrongadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Wrongadjmtx_Counter1 = 0
    Flame_Wrongadjmtx_Counter2 = 0
    tmpcounter = 0
    Flame_Wrongadjmtx_RANGE1 = []
    Flame_Wrongadjmtx_RANGE2 = []
    #-----------The Function Itself-----------#
    #Iterates through the input (Potential New Exons)
    for Flame_Wrongadjmtx_COUNT1 in INCONGRUENTADJMTXINPUT: 
        Flame_Wrongadjmtx_Counter1 = 0
        Flame_Wrongadjmtx_Counter2 = 0
        #Crosschecks with the Reference (Splice Candidates (X > 1% frequency). This loop is for the Row position.
        for Flame_Wrongadjmtx_COUNT2 in REF: 
            if (int(Flame_Wrongadjmtx_COUNT1.split("-")[0])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[0])+ADJMTXRANGESIZE):
                break
            elif (not (int(Flame_Wrongadjmtx_COUNT1.split("-")[0])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[0])+ADJMTXRANGESIZE)):
                Flame_Wrongadjmtx_Counter1 += 1
        #Crosschecks with the Reference (Splice Candidates (X > 1% frequency). This loop is for the Col position. 
        for Flame_Wrongadjmtx_COUNT2 in REF: 
            if (int(Flame_Wrongadjmtx_COUNT1.split("-")[1])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[1])+ADJMTXRANGESIZE):
                break
            elif (not (int(Flame_Wrongadjmtx_COUNT1.split("-")[1])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[1])+ADJMTXRANGESIZE)):
                Flame_Wrongadjmtx_Counter2 += 1
        #Add to the cordinates within the input-"empty"-adjacency matrix.
        Flame_Wrongadjmtx_EMPTYADJMTX[Flame_Wrongadjmtx_Counter2][Flame_Wrongadjmtx_Counter1] += 1
        #-----------The Progress Bar Start-----------#
        Flame_Wrongadjmtx_PROGRESSCOUNT1 += 1
        Flame_Wrongadjmtx_PROGRESSCOUNT2 += 1
        if Flame_Wrongadjmtx_PROGRESSCOUNT2 != Flame_Wrongadjmtx_PROGRESSMAX:
            if Flame_Wrongadjmtx_PROGRESSCOUNT1 >= int(round(Flame_Wrongadjmtx_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_Wrongadjmtx_PROGRESSCOUNT2 / Flame_Wrongadjmtx_PROGRESSMAX)
                Flame_Wrongadjmtx_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_Wrongadjmtx_PROGRESSCOUNT2 == Flame_Wrongadjmtx_PROGRESSMAX:
            PROGRESSBAR(Flame_Wrongadjmtx_PROGRESSCOUNT2/Flame_Wrongadjmtx_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_Wrongadjmtx_EMPTYADJMTX #Output Object Type: Nested List.

#Central Function to Seperate the potential novel exon region from the Incongruent Dataset.
def INCONGRUENTSEPERATORFUNC(INCONGRUENTSEPINPUT, REF, MODE):
    #-----------Outside counters and variables-----------#
    #Flame_Inconsep_PROGRESSMAX moved to after the filtering step, added patch 0.1.5
    Flame_Inconsep_FILTER_COUNT1 = 0
    Flame_Inconsep_FILTER_DELLIST = []
    Flame_Inconsep_FILTER_FLAG = True
    Flame_Inconsep_PROGRESSCOUNT1 = 0
    Flame_Inconsep_PROGRESSCOUNT2 = 0
    Flame_Inconsep_NAMESEP = []
    Flame_Inconsep_STRT_ALL = []
    Flame_Inconsep_POTENTIAL_EXON = []
    #-----------The Function Itself-----------#
    if not MODE:
        for i in [x.split(",") for x in INCONGRUENTSEPINPUT]:
            Flame_Inconsep_FILTER_FLAG = True
            for j in i:
                try:
                    int(j)
                    Flame_Inconsep_FILTER_FLAG = False
                    break
                except ValueError:
                    pass
            if Flame_Inconsep_FILTER_FLAG:
                Flame_Inconsep_FILTER_DELLIST.append(Flame_Inconsep_FILTER_COUNT1)
            Flame_Inconsep_FILTER_COUNT1 += 1
        for i in list(reversed(Flame_Inconsep_FILTER_DELLIST)):
            del INCONGRUENTSEPINPUT[i]
    else:
        pass
    Flame_Inconsep_PROGRESSMAX = sum(1 for Flame_Inconsep_LINE in INCONGRUENTSEPINPUT)
    #Create an reference and database
    for Flame_Inconsep_COUNT1 in REF:
        Flame_Inconsep_NAMESEP.append(Flame_Inconsep_COUNT1[0])
    #Create a list containing only the Non-reference exon ranges:
    for Flame_Inconsep_COUNT2 in INCONGRUENTSEPINPUT:
        Flame_Inconsep_STRT_ALL = Flame_Inconsep_COUNT2.split(",")[:-1] #FIX: If one changes the "TMPSTRING" function, then one can avoid the use of this Remove last element function"
        for Flame_Inconsep_COUNT3 in Flame_Inconsep_STRT_ALL: #Only if Flame_Inconsep_COUNT3 does not match, append the Exon range to a list that will be returned.
            if Flame_Inconsep_COUNT3 in Flame_Inconsep_NAMESEP:
                pass
            elif Flame_Inconsep_COUNT3 not in Flame_Inconsep_NAMESEP:
                Flame_Inconsep_POTENTIAL_EXON.append(Flame_Inconsep_COUNT3)
        #-----------The Progress Bar Start-----------#
        Flame_Inconsep_PROGRESSCOUNT1 += 1
        Flame_Inconsep_PROGRESSCOUNT2 += 1
        if Flame_Inconsep_PROGRESSCOUNT2 != Flame_Inconsep_PROGRESSMAX:
            if Flame_Inconsep_PROGRESSCOUNT1 >= int(round(Flame_Inconsep_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_Inconsep_PROGRESSCOUNT2 / Flame_Inconsep_PROGRESSMAX)
                Flame_Inconsep_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_Inconsep_PROGRESSCOUNT2 == Flame_Inconsep_PROGRESSMAX:
            PROGRESSBAR(Flame_Inconsep_PROGRESSCOUNT2/Flame_Inconsep_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_Inconsep_POTENTIAL_EXON #Output Object Type: List


#Central Function to measure the frequency of potentially novel splice site signals (S3).
def FREQUENCYSITEFUNC(FREQUENCYSITEINPUT, REF, FREQUENCYSITERANGESIZE, FREQUENCYSITEWINDOWSIZE):
    #-----------Outside counters and variables-----------#
    Flame_FrequencySite_PROGRESSMAX = sum(1 for Flame_FrequencySite_LINE in FREQUENCYSITEINPUT)
    Flame_FrequencySite_PROGRESSCOUNT1 = 0
    Flame_FrequencySite_PROGRESSCOUNT2 = 0
    Flame_FrequencySite_GENEEMULATION = []
    Flame_FrequencySite_COUNT1RANGE = []
    Flame_FrequencySite_START = []
    Flame_FrequencySite_END = []
    Flame_FrequencySite_MAX = 0
    Flame_FrequencySite_Counter1 = 0
    #-----------The Function Itself-----------#
    if len(REF) == 0: #If-statement that checks whether the Reference is empty or not.
        for Flame_FrequencySite_COUNT1 in FREQUENCYSITEINPUT:
            if int(Flame_FrequencySite_COUNT1.split("-")[1]) > Flame_FrequencySite_MAX:
                Flame_FrequencySite_MAX = int(Flame_FrequencySite_COUNT1.split("-")[1])
            elif int(Flame_FrequencySite_COUNT1.split("-")[1]) <= Flame_FrequencySite_MAX:
                pass
    elif len(REF) > 0: #Going through the reference and singeling out the last possible "genomic"-position.
        for Flame_FrequencySite_COUNT1 in REF:
            if REF[Flame_FrequencySite_Counter1][4] > Flame_FrequencySite_MAX:
                Flame_FrequencySite_MAX = REF[Flame_FrequencySite_Counter1][4]
                Flame_FrequencySite_Counter1 += 1
            elif REF[Flame_FrequencySite_Counter1][4] <= Flame_FrequencySite_MAX:
                Flame_FrequencySite_Counter1 += 1
    #Create an empty dataset that starts at "genomic"-position 1 and ends at the last possible "genomic"-position dependent on the GTF file.
    for Flame_FrequencySite_COUNT2 in range((Flame_FrequencySite_MAX) +
                                            FREQUENCYSITERANGESIZE): #Create a smaller list for only the gene range?
        Flame_FrequencySite_GENEEMULATION.append(0) #FIX: Make X.append([0]*Y) function instead?
    #Go through each potential novel exon-range and increase the start&stop position's, and their adjecent nucleotides, value. This to highlight potential novel S3.
    for Flame_FrequencySite_COUNT2 in FREQUENCYSITEINPUT:
        Flame_FrequencySite_COUNT1RANGE = Flame_FrequencySite_COUNT2.split("-") #This splits the two numbers (Start and End). This makes it so one cannot see the connections.
        if int(Flame_FrequencySite_COUNT1RANGE[1]) < ((Flame_FrequencySite_MAX) +
                                                      FREQUENCYSITERANGESIZE):
            Flame_FrequencySite_START = list(range(int(Flame_FrequencySite_COUNT1RANGE[0])-FREQUENCYSITEWINDOWSIZE,
                                                   int(Flame_FrequencySite_COUNT1RANGE[0])+FREQUENCYSITEWINDOWSIZE)) #FIX: Add "+1" as they only have 4 values instead of 5 values (X +- 2)?
            Flame_FrequencySite_END = list(range(int(Flame_FrequencySite_COUNT1RANGE[1])-FREQUENCYSITEWINDOWSIZE,
                                                 int(Flame_FrequencySite_COUNT1RANGE[1])+FREQUENCYSITEWINDOWSIZE))
            for Flame_FrequencySite_COUNT3 in Flame_FrequencySite_START:
                Flame_FrequencySite_GENEEMULATION[Flame_FrequencySite_COUNT3] += 1
            for Flame_FrequencySite_COUNT3 in Flame_FrequencySite_END:
                Flame_FrequencySite_GENEEMULATION[Flame_FrequencySite_COUNT3] += 1
        else:
            pass
        #-----------The Progress Bar Start-----------#
        Flame_FrequencySite_PROGRESSCOUNT1 += 1
        Flame_FrequencySite_PROGRESSCOUNT2 += 1
        if Flame_FrequencySite_PROGRESSCOUNT2 != Flame_FrequencySite_PROGRESSMAX:
            if Flame_FrequencySite_PROGRESSCOUNT1 >= int(round(Flame_FrequencySite_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_FrequencySite_PROGRESSCOUNT2 / Flame_FrequencySite_PROGRESSMAX)
                Flame_FrequencySite_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_FrequencySite_PROGRESSCOUNT2 == Flame_FrequencySite_PROGRESSMAX:
            PROGRESSBAR(Flame_FrequencySite_PROGRESSCOUNT2/Flame_FrequencySite_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_FrequencySite_GENEEMULATION #Output Object Type: List

def FREQUENCYSITEOPTIFUNC(SIMULATEDGENOME):
    #-----------Outside counters and variables-----------#
    Flame_FrequencyOpti_PROGRESSMAX = sum(1 for Flame_FrequencyOpti_LINE in SIMULATEDGENOME)
    Flame_FrequencyOpti_PROGRESSCOUNT1 = 0
    Flame_FrequencyOpti_PROGRESSCOUNT2 = 0
    Flame_FrequencyOpti_NONZEROLIST = []
    Flame_FrequencyOpti_Counter1 = 0 #Poisition index counter
    
    #-----------The Function Itself-----------#
    #Going through the "genome" (Simulated) and finding where the gene of interest starts for optimizations sake, should have a time complexity function = O(N)
    for Flame_FrequencyOpti_COUNTER1 in SIMULATEDGENOME:
        if Flame_FrequencyOpti_COUNTER1 > 0:
            Flame_FrequencyOpti_NONZEROLIST.append((Flame_FrequencyOpti_Counter1,
                                                    Flame_FrequencyOpti_COUNTER1))
        else:
            pass
        Flame_FrequencyOpti_Counter1 += 1
        #-----------The Progress Bar Start-----------#
        Flame_FrequencyOpti_PROGRESSCOUNT1 += 1
        Flame_FrequencyOpti_PROGRESSCOUNT2 += 1
        if Flame_FrequencyOpti_PROGRESSCOUNT2 != Flame_FrequencyOpti_PROGRESSMAX:
            if Flame_FrequencyOpti_PROGRESSCOUNT1 >= int(round(Flame_FrequencyOpti_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_FrequencyOpti_PROGRESSCOUNT2 / Flame_FrequencyOpti_PROGRESSMAX)
                Flame_FrequencyOpti_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_FrequencyOpti_PROGRESSCOUNT2 == Flame_FrequencyOpti_PROGRESSMAX:
            PROGRESSBAR(Flame_FrequencyOpti_PROGRESSCOUNT2/Flame_FrequencyOpti_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_FrequencyOpti_NONZEROLIST #Output Object Type: List of tuples


#Central Function to set a threshold and high-light the S3 that exceed said threshold.
def FREQUENCYTHRESHFUNC(SIMULATEDGENOME, SIMULATEDGENOME_START, PERCENTTHRESH, INPUTLIST):
    #-----------Outside counters and variables-----------#
    Flame_FrequencyThresh_PROGRESSMAX = sum(1 for Flame_FrequencyThresh_LINE in SIMULATEDGENOME)
    Flame_FrequencyThresh_PROGRESSCOUNT1 = SIMULATEDGENOME_START[0][0]
    Flame_FrequencyThresh_PROGRESSCOUNT2 = SIMULATEDGENOME_START[0][0]
    Flame_FrequencyThresh_THRESHOLD = float()
    Flame_FrequencyThresh_SPLICECANDIDATES = []
    Flame_FrequencyThresh_Counter1 = SIMULATEDGENOME_START[0][0]
    #-----------The Function Itself-----------#
    #Going through each "genomic"-nucleotide position to filter out any position that does not exceed the threshold leaving only high-frequency potential novel S3.
    for Flame_FrequencyThresh_COUNT1 in SIMULATEDGENOME[SIMULATEDGENOME_START[0][0]:len(SIMULATEDGENOME)]:
        Flame_FrequencyThresh_THRESHOLD = float(PERCENTTHRESH*min(len(INPUTLIST), 10000)) #Remove zeros?
        if Flame_FrequencyThresh_COUNT1 > Flame_FrequencyThresh_THRESHOLD:
            Flame_FrequencyThresh_SPLICECANDIDATES.append([Flame_FrequencyThresh_Counter1,
                                                           Flame_FrequencyThresh_COUNT1,
                                                           round((float(Flame_FrequencyThresh_COUNT1)/len(INPUTLIST))*100, 2)])
            Flame_FrequencyThresh_Counter1 += 1
        else:
            Flame_FrequencyThresh_Counter1 += 1
        #-----------The Progress Bar Start-----------#
        Flame_FrequencyThresh_PROGRESSCOUNT1 += 1
        Flame_FrequencyThresh_PROGRESSCOUNT2 += 1
        if Flame_FrequencyThresh_PROGRESSCOUNT2 != Flame_FrequencyThresh_PROGRESSMAX:
            if Flame_FrequencyThresh_PROGRESSCOUNT1 >= int(round(Flame_FrequencyThresh_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_FrequencyThresh_PROGRESSCOUNT2 / Flame_FrequencyThresh_PROGRESSMAX)
                Flame_FrequencyThresh_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_FrequencyThresh_PROGRESSCOUNT2 == Flame_FrequencyThresh_PROGRESSMAX:
            PROGRESSBAR(Flame_FrequencyThresh_PROGRESSCOUNT2/Flame_FrequencyThresh_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_FrequencyThresh_SPLICECANDIDATES #Output Object Type: List


#Central Function to survey adjacent nucleotides of the highlighted potential novel S3 for GU-AG signals. 
def SPLICESIGNALFUNC(SPLICESIGNALINPUT, REF):
    #-----------Outside counters and variables-----------#
    Flame_SpliceSignal_PROGRESSMAX = sum(1 for Flame_SpliceSignal_LINE in SPLICESIGNALINPUT)
    Flame_SpliceSignal_PROGRESSCOUNT1 = 0
    Flame_SpliceSignal_PROGRESSCOUNT2 = 0
    Flame_SpliceSignal_SPLICESIGNALINPUT = SPLICESIGNALINPUT
    Flame_SpliceSignal_WINDOW = []
    Flame_SpliceSignal_GU = False
    Flame_SpliceSignal_AG = False
    #-----------The Function Itself-----------#
    #Going through each highlighted S3
    for Flame_SpliceSignal_COUNT1 in Flame_SpliceSignal_SPLICESIGNALINPUT:
        Flame_SpliceSignal_GU = False
        Flame_SpliceSignal_AG = False
        Flame_SpliceSignal_WINDOW = REF[(int(Flame_SpliceSignal_COUNT1[0])-3): #FIX: Make window customizeable?
                                        (int(Flame_SpliceSignal_COUNT1[0])+4)]
        #Check for either GU/GT and/or AG signal.
        if "GT" in Flame_SpliceSignal_WINDOW[3:7]:
            Flame_SpliceSignal_GU = True
        if "AG" in Flame_SpliceSignal_WINDOW[0:4]:
            Flame_SpliceSignal_AG = True
        #List out the 4 different outcomes: None, GU, AG, Both.
        if Flame_SpliceSignal_GU == False and Flame_SpliceSignal_AG == False:
            Flame_SpliceSignal_COUNT1.append("None")
        elif Flame_SpliceSignal_GU == True and Flame_SpliceSignal_AG == False:
            Flame_SpliceSignal_COUNT1.append("GU")
        elif Flame_SpliceSignal_GU == False and Flame_SpliceSignal_AG == True:
            Flame_SpliceSignal_COUNT1.append("AG")
        elif Flame_SpliceSignal_GU == True and Flame_SpliceSignal_AG == True:
            Flame_SpliceSignal_COUNT1.append("Both")
        #-----------The Progress Bar Start-----------#
        Flame_SpliceSignal_PROGRESSCOUNT1 += 1
        Flame_SpliceSignal_PROGRESSCOUNT2 += 1
        if Flame_SpliceSignal_PROGRESSCOUNT2 != Flame_SpliceSignal_PROGRESSMAX:
            if Flame_SpliceSignal_PROGRESSCOUNT1 >= int(round(Flame_SpliceSignal_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_SpliceSignal_PROGRESSCOUNT2 / Flame_SpliceSignal_PROGRESSMAX)
                Flame_SpliceSignal_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_SpliceSignal_PROGRESSCOUNT2 == Flame_SpliceSignal_PROGRESSMAX:
            PROGRESSBAR(Flame_SpliceSignal_PROGRESSCOUNT2/Flame_SpliceSignal_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    return Flame_SpliceSignal_SPLICESIGNALINPUT #Output Object Type: List

#Central Function to involve short-read confirmation using the short-read BAM/SAM file.
def SHORTREADFUNC(SHORTREADINPUT1, SHORTREADCANDIDATES, REF):
    #-----------Outside counters and variables-----------#
    Flame_ShortRead_PROGRESSMAX = sum(1 for Flame_ShortRead_LINE in SHORTREADINPUT1)
    pysam.HTSFile.reset(SHORTREADINPUT1)
    Flame_ShortRead_PROGRESSCOUNT1 = 0
    Flame_ShortRead_PROGRESSCOUNT2 = 0
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
    #Going through each aligned short-read from the BAM/SAM
    for Flame_ShortRead_COUNT1 in SHORTREADINPUT1:
        #For ease of use: Singling out all the relevant information for the BAM/SAM read:
        Flame_ShortRead_START = int(str(Flame_ShortRead_COUNT1).split("\t")[3])
        Flame_ShortRead_CIGAR = str(Flame_ShortRead_COUNT1).split("\t")[5]
        Flame_ShortRead_CHROM = str(Flame_ShortRead_COUNT1).split("\t")[2]
        if "N" in Flame_ShortRead_CIGAR:
            Fĺame_ShortRead_LENGTH = sum(list(map(int, re.findall(r'\d+',
                                                                  Flame_ShortRead_CIGAR))))
            #For ease of use: Collapsing the characteristics into a single list.
            Flame_ShortRead_COMB = [Flame_ShortRead_START,
                                    Flame_ShortRead_CIGAR,
                                    (Flame_ShortRead_START+Fĺame_ShortRead_LENGTH)]
            Flame_ShortRead_MATCH1 = None
            Flame_ShortRead_MATCH2 = None
            Flame_ShortRead_MATCH3 = None
            Flame_ShortRead_MATCH4 = None
            Flame_ShortRead_ITEMS = []
            #Central if-statement to localize how many splice events a single short-read sequence has (4, 3, 2 or 1) and spliting the CIGARSTRING-delimiters between them.
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
            else:
                print("Error, CIGAR-String:", Flame_ShortRead_COMB)
            #For ease of use: Collapsing the four different options of splice event into a single variable.
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
            #To remove Softclipped parts.
            if Flame_ShortRead_ITEMS[1] == "S":
                del Flame_ShortRead_ITEMS[:2]
            if Flame_ShortRead_ITEMS[-1] == "S":
                del Flame_ShortRead_ITEMS[-2:]            
            Flame_ShortRead_TMPSTART1 = Flame_ShortRead_START
            #To iterate through all of the seperated CIGARSTRING-parts.
            for Flame_ShortRead_COUNT1 in range(len(Flame_ShortRead_ITEMS)):
                Flame_ShortRead_CIGARCOUNT = 0
                Flame_ShortRead_CIGAROPERATOR = 0
                #To filter out any reads that has an uneven number of CIGARSTRING-events due to splice-events requiring uneven number of CIGARSTRING elements to exist (MATCH-SPLICE-MATCH, MATCH-INSERTION-MATCH-SPLICE-MATCH, and so on)
                if Flame_ShortRead_COUNT1 % 2 != 0 and Flame_ShortRead_COUNT1 > 0:
                    Flame_ShortRead_CIGARCOUNT = int(Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1-1]) #The Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1-1] is not a mistake.
                    Flame_ShortRead_CIGAROPERATOR = Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1]
                    #To filter out any read that does not have the CIGAR-operator of "N" which denotes splice event.
                    if Flame_ShortRead_CIGAROPERATOR != "N":
                        Flame_ShortRead_TMPSTART1 = Flame_ShortRead_TMPSTART1 + Flame_ShortRead_CIGARCOUNT
                        pass
                    elif Flame_ShortRead_CIGAROPERATOR == "N":
                        Flame_ShortRead_TMPSTART2 = Flame_ShortRead_TMPSTART1
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
        #-----------The Progress Bar Start-----------#
        Flame_ShortRead_PROGRESSCOUNT1 += 1
        Flame_ShortRead_PROGRESSCOUNT2 += 1
        if Flame_ShortRead_PROGRESSCOUNT2 != Flame_ShortRead_PROGRESSMAX:
            if Flame_ShortRead_PROGRESSCOUNT1 >= int(round(Flame_ShortRead_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_ShortRead_PROGRESSCOUNT2 / Flame_ShortRead_PROGRESSMAX)
                Flame_ShortRead_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_ShortRead_PROGRESSCOUNT2 == Flame_ShortRead_PROGRESSMAX:
            PROGRESSBAR(Flame_ShortRead_PROGRESSCOUNT2/Flame_ShortRead_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar End-----------#
    #The function to print out the results from previous function by appending either a value to denote the number of supporting short-reads or N/A as in Not-Available.
    for Flame_ShortRead_COUNT2 in Flame_ShortRead_CANDIDATES:
        #-----------Outside counters and variables-----------#
        Flame_ShortRead_Counter1 = 0
        #-----------The Function Itself-----------#
        for k, v in Flame_ShortRead_SPLICESITECOUNT.items():
            if Flame_ShortRead_COUNT2[0] == k:
                Flame_ShortRead_COUNT2.extend([v])
                break
            elif Flame_ShortRead_COUNT2[0] != k:
                Flame_ShortRead_Counter1 += 1
        if Flame_ShortRead_Counter1 == len(Flame_ShortRead_SPLICESITECOUNT):
            Flame_ShortRead_COUNT2.extend(["N/A"])
    return Flame_ShortRead_CANDIDATES, Flame_ShortRead_SPLICESITECOUNT



#def main() #FIX: Create a Main function if you have everything ready and just want to run everything in one go.
#if __name__ == "__main__":
# main()
