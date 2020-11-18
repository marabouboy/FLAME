#!/usr/bin/python3
import re, time, sys
import argparse
import pysam

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
    Flame_Progress_TEXT = "\rStatus: [{0}] {1}% {2}".format( "#"*Flame_Progress_BLOCK +
                                                              "-"*(Flame_Progress_BARLENGTH-Flame_Progress_BLOCK),
                                                              int(round(PROGRESS,2)*100),
                                                              Flame_Progress_STATUS)
    sys.stdout.write(Flame_Progress_TEXT)
    sys.stdout.flush()

    
def CREATEREFFUNC(CREATEREFRINPUT, CREATEREFNAME): #CHANGE IT SO IT DOES NOT MATTER IF IT IS GTF OR GFF?
    #-----------Outside counters and variables-----------#
    Flame_CreateRef_REFERENCE = []
    Flame_CreateRef_BACKLOGREPEATS = []
    Flame_CreateRef_NAMELIST = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
                                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
                                "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                                "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",
                                "51", "52", "53", "54", "55", "56", "57", "58", "59", "60",
                                "61", "62", "63", "64", "65", "66", "67", "68", "69", "70",
                                "71", "72", "73", "74", "75", "76", "77", "78", "79", "80",
                                "81", "82", "83", "84", "85", "86", "87", "88", "89", "90",
                                "91", "92", "93", "94", "95", "96", "97", "98", "99", "100",
                                "101", "102", "103", "104", "105", "106", "107", "108", "109", "110",
                                "111", "112", "113", "114", "115", "116", "117", "118", "119", "120",
                                "121", "122", "123", "124", "125", "126", "127", "128", "129", "130",
                                "131", "132", "133", "134", "135", "136", "137", "138", "139", "140",
                                "141", "142", "143", "144", "145", "146", "147", "148", "149", "150"]
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
            ((str(CREATEREFNAME) + "-") not in Flame_CreateRef_COUNT1.split(";")[1]) and
            len(Flame_CreateRef_COUNT1) >= 1):
            Flame_CreateRef_START = int(Flame_CreateRef_COUNT1.split("\t")[3]) 
            Flame_CreateRef_STOP = int(Flame_CreateRef_COUNT1.split("\t")[4]) 
            Flame_CreateRef_LEN = Flame_CreateRef_STOP - Flame_CreateRef_START
            Flame_CreateRef_COMB = ["",
                                    Flame_CreateRef_START,
                                    Flame_CreateRef_LEN,
                                    Flame_CreateRef_STOP]
            if Flame_CreateRef_COMB in Flame_CreateRef_BACKLOGREPEATS:
                pass
            elif Flame_CreateRef_COMB not in Flame_CreateRef_BACKLOGREPEATS:
                Flame_CreateRef_BACKLOGREPEATS.append(Flame_CreateRef_COMB)
    Flame_CreateRef_REFERENCE = sorted(Flame_CreateRef_BACKLOGREPEATS,
                                       key = lambda x: x[1])
    for Flame_CreateRef_COUNT1 in Flame_CreateRef_REFERENCE:
        Flame_CreateRef_COUNT1[0] = Flame_CreateRef_NAMELIST[Flame_CreateRef_NAMECOUNT]
        Flame_CreateRef_NAMECOUNT += 1           
    return Flame_CreateRef_REFERENCE


#Central Filter function that will sort the reads into two different lists (FILTERCORRECTOUTPUT, FILTERFAULTYOUTPUT) depending if they have only contain perfectly matching exons to the reference or if they have a "faulty" exon.
def FILTERFUNC(FILTERINPUT1, FILTERINPUT2, REF, FILTERRANGESIZE):
    #-----------Outside counters and variables-----------#
    Flame_Filter_PROGRESSMAX = sum(1 for Flame_Filter_LINE in FILTERINPUT2)
    Flame_Filter_PROGRESSCOUNT1 = 0
    Flame_Filter_PROGRESSCOUNT2 = 0
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
    #-----------The Function Itself-----------#
    for Flame_Filter_SINGLE_READ in FILTERINPUT1.read().split("\n"): #Loop: Going through each read. #Change it so that the input read (As an outside read file) is read outside?
        if len(Flame_Filter_SINGLE_READ) >= 1: #If statement that filters out empty reads.
            #For ease of use (Below): Singling out all the relevant information from the bed12 read:
            Flame_Filter_SINGLE_READ_STRTPNT = int(Flame_Filter_SINGLE_READ.split("\t")[1])
            Flame_Filter_SINGLE_READ_STRT_ALL = Flame_Filter_SINGLE_READ.split("\t")[11]
            Flame_Filter_SINGLE_READ_LEN_ALL =  Flame_Filter_SINGLE_READ.split("\t")[10]
            Flame_Filter_SINGLE_READ_COUNT = int(Flame_Filter_SINGLE_READ.split("\t")[9])
            Flame_Filter_EXONCOUNTER = 0
            for Flame_Filter_COUNT1 in range(Flame_Filter_SINGLE_READ_COUNT): #Loop: Going through each exon.
                Flame_Filter_Counter = 0
                Flame_Filter_FAULTYCOUNTER = 0
                #For ease of use (Below): Specifying each characteristic (Start, Length, Stop) for each exon in each read.
                Flame_Filter_SPLICESTART = (int(Flame_Filter_SINGLE_READ_STRT_ALL.split(",")[Flame_Filter_COUNT1]) +
                                            Flame_Filter_SINGLE_READ_STRTPNT)
                Flame_Filter_SPLICELEN = int(Flame_Filter_SINGLE_READ_LEN_ALL.split(",")[Flame_Filter_COUNT1])
                Flame_Filter_SPLICESTOP = (Flame_Filter_SPLICESTART + Flame_Filter_SPLICELEN)
                #For ease of use (Below): Collapsing the characteristics into a single list.
                Flame_Filter_SPLICECOMB = [Flame_Filter_SPLICESTART,
                                           Flame_Filter_SPLICELEN,
                                           Flame_Filter_SPLICESTOP]
                while Flame_Filter_FAULTYCOUNTER != len(REF): #Loop: Going through each reference exon trying to match with the exon.
                    if ((REF[Flame_Filter_Counter][1]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[0] <= REF[Flame_Filter_Counter][1]+FILTERRANGESIZE) and
                        (REF[Flame_Filter_Counter][2]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[1] <= REF[Flame_Filter_Counter][2]+FILTERRANGESIZE) and
                        (REF[Flame_Filter_Counter][3]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[2] <= REF[Flame_Filter_Counter][3]+FILTERRANGESIZE)): #If statement that will be flagged if the Exon matches the reference, otherwise, keep looping until you reach the end of the reference file.
                        Flame_Filter_EXONCOUNTER += 1
                        Flame_Filter_Counter += 1                        
                        break
                    elif ((not (REF[Flame_Filter_Counter][1]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[0] <= REF[Flame_Filter_Counter][1]+FILTERRANGESIZE)) or
                          (not (REF[Flame_Filter_Counter][2]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[1] <= REF[Flame_Filter_Counter][2]+FILTERRANGESIZE)) or
                          (not (REF[Flame_Filter_Counter][3]-FILTERRANGESIZE <= Flame_Filter_SPLICECOMB[2] <= REF[Flame_Filter_Counter][3]+FILTERRANGESIZE))):
                        Flame_Filter_FAULTYCOUNTER += 1
                        Flame_Filter_Counter += 1
                        continue
            if Flame_Filter_EXONCOUNTER == Flame_Filter_SINGLE_READ_COUNT: #If statement that will determine if the read contains an exon that is not annotated by the reference.
                FILTERCORRECTOUTPUT.append(Flame_Filter_SINGLE_READ)
            elif Flame_Filter_EXONCOUNTER != Flame_Filter_SINGLE_READ_COUNT:
                FILTERFAULTYOUTPUT.append(Flame_Filter_SINGLE_READ)
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    return FILTERCORRECTOUTPUT, FILTERFAULTYOUTPUT #Output Object Type: [List, List]


#Central Translate function that contain the translate function for both the Correct Reads as well as the Faulty Reads.
def TRANSLATEFUNC(REF, TRANSLATERANGESIZE, TRANSLATECORRECTINPUT = 0, TRANSLATEFAULTYINPUT = 0):    
    if TRANSLATECORRECTINPUT != []: #If statement that will activate depending if the input (Input[2]) is not empty.
        print("-----------\tInitiate Translate Function, Corrected\t\t\t\t-----------")
        #-----------Outside counters and variables-----------#
        Flame_Translate_PROGRESSMAX = sum(1 for Flame_Translate_LINE in TRANSLATECORRECTINPUT)
        Flame_Translate_PROGRESSCOUNT1 = 0
        Flame_Translate_PROGRESSCOUNT2 = 0
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
        Flame_Translate_CORRECT = []
        #-----------The Function Itself-----------#
        for Flame_Translate_SINGLE_READ in TRANSLATECORRECTINPUT: #Loop: Going through each correct read. This is reliant on that either you have sorted the reads between "Correct" and "Faulty" from the previous function (FILTERFUNC) or that you have a reliable dataset to use.
            #For ease of use (Below): Singling out all the relevant information for the bed12 read:
            Flame_Translate_STRTPNT = int(Flame_Translate_SINGLE_READ.split("\t")[1])
            Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
            Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
            Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
            Flame_Translate_Counter1 = 0
            Flame_Translate_Counter2 = 0
            Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
            while Flame_Translate_Counter1 < Flame_Translate_SINGLE_READ_COUNT: #Loop: Iterate through every single exon (Flame_Translate_Counter1)
                Flame_Translate_SPLICESTART = (int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_Counter1]) +
                                               Flame_Translate_STRTPNT)
                Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_Counter1])
                Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                #For ease of use (Below): Collapsing the characteristics into a single list.
                Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                              Flame_Translate_SPLICELEN,
                                              Flame_Translate_SPLICESTOP]
                if ((REF[Flame_Translate_Counter2][1]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter2][1]+TRANSLATERANGESIZE) and
                    (REF[Flame_Translate_Counter2][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter2][2]+TRANSLATERANGESIZE) and
                    (REF[Flame_Translate_Counter2][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter2][3]+TRANSLATERANGESIZE)):
                    Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter2][0] + ",") #FIX: Could insert a part that removes the last dash to help downstream but does not matter in the grander scheme. Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element.
                    Flame_Translate_Counter1 += 1
                    continue
                elif ((not (REF[Flame_Translate_Counter2][1]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter2][1]+TRANSLATERANGESIZE)) or
                      (not (REF[Flame_Translate_Counter2][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter2][2]+TRANSLATERANGESIZE)) or
                      (not (REF[Flame_Translate_Counter2][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter2][3]+TRANSLATERANGESIZE))):
                    Flame_Translate_Counter2 += 1
                    continue
            Flame_Translate_CORRECT.append(Flame_Translate_TEMPORARY_STRING) #Function to output the translated read into a list (Flame_Translate_CORRECT) and recycle to next read in the central for-loop.
            #-----------The Progress Bar-----------#
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
            #-----------The Progress Bar-----------#
    elif TRANSLATECORRECTINPUT == []: #If statement that will activate depending if the input (Input[2]) is empty. This is to just create any output for the function.
        print("-----------\tInitiate Translate Function, Corrected\t\t\t\t-----------")
        print("Status: [] Error... Empty Input")
        Flame_Translate_CORRECT = []
        
    if TRANSLATEFAULTYINPUT != []: #If statement that will activate depending if the input (Input[3]) is not empty.
        print("-----------\tInitiate Translate Function, Faulty\t\t\t\t-----------")
        #-----------Outside counters and variables-----------#
        Flame_Translate_PROGRESSMAX = sum(1 for Flame_Translate_LINE in TRANSLATEFAULTYINPUT)
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
        Flame_Translate_FAULTY = []
        #-----------The Function Itself-----------#
        for Flame_Translate_SINGLE_READ in TRANSLATEFAULTYINPUT: #Loop: Going through each faulty read. This is reliant on that either you have sorted the reads between "Correct" and "Faulty" from the previous function (FILTERFUNC) or that you have a reliable dataset to use.
            #For ease of use (Below): Singling out all the relevant information for the bed12 read:
            Flame_Translate_STRTPNT = int(Flame_Translate_SINGLE_READ.split("\t")[1])
            Flame_Translate_SINGLE_READ_STRT_ALL = Flame_Translate_SINGLE_READ.split("\t")[11]
            Flame_Translate_SINGLE_READ_LEN_ALL = Flame_Translate_SINGLE_READ.split("\t")[10]
            Flame_Translate_SINGLE_READ_COUNT = int(Flame_Translate_SINGLE_READ.split("\t")[9])
            Flame_Translate_Counter3 = 0
            Flame_Translate_Counter4 = 0
            Flame_Translate_TEMPORARY_STRING = "" #Reset the temporary string into an empty one.
            for Flame_Translate_COUNT1 in range(Flame_Translate_SINGLE_READ_COUNT): #Loop: Iterate through every single exon (Flame_Translate_COUNT1)
                Flame_Translate_Counter3 = 0
                Flame_Translate_Counter4 = 0
                Flame_Translate_SPLICESTART = (int(Flame_Translate_SINGLE_READ_STRT_ALL.split(",")[Flame_Translate_COUNT1]) +
                                               Flame_Translate_STRTPNT)
                Flame_Translate_SPLICELEN = int(Flame_Translate_SINGLE_READ_LEN_ALL.split(",")[Flame_Translate_COUNT1])
                Flame_Translate_SPLICESTOP = (Flame_Translate_SPLICESTART + Flame_Translate_SPLICELEN)
                #For ease of use (Below): Collapsing the characteristics into a single list.
                Flame_Translate_SPLICECOMB = [Flame_Translate_SPLICESTART,
                                              Flame_Translate_SPLICELEN,
                                              Flame_Translate_SPLICESTOP]
                while Flame_Translate_Counter4 != len(REF): #Loop: Forcing to run through all the Exons specified in the reference unless match in which it breaks the loop in order to optimize the running time.
                    if ((REF[Flame_Translate_Counter3][1]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter3][1]+TRANSLATERANGESIZE) and
                        (REF[Flame_Translate_Counter3][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter3][2]+TRANSLATERANGESIZE) and
                        (REF[Flame_Translate_Counter3][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter3][3]+TRANSLATERANGESIZE)):
                        Flame_Translate_Counter2 += 1
                        Flame_Translate_TEMPORARY_STRING += str(REF[Flame_Translate_Counter3][0] + ",") #Make it more efficient by changing it into a list with ".join" function that also allows for the removal of the last element. list1 = ['1','2','3','4'], s = "-", s = s.join(list1) 
                        Flame_Translate_Counter3 += 1
                        Flame_Translate_BOOLEAN = False 
                        break
                    elif ((not (REF[Flame_Translate_Counter3][1]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[0] <= REF[Flame_Translate_Counter3][1]+TRANSLATERANGESIZE)) or
                          (not (REF[Flame_Translate_Counter3][2]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[1] <= REF[Flame_Translate_Counter3][2]+TRANSLATERANGESIZE)) or
                          (not (REF[Flame_Translate_Counter3][3]-TRANSLATERANGESIZE <= Flame_Translate_SPLICECOMB[2] <= REF[Flame_Translate_Counter3][3]+TRANSLATERANGESIZE))):
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
            #-----------The Progress Bar-----------#
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
            #-----------The Progress Bar-----------#
    elif TRANSLATEFAULTYINPUT == []: #If statement that will activate depending if the input (Input[2]) is empty. This is to just create any output for the function.
        print("-----------\tInitiate Translate Function, Faulty\t\t\t\t-----------")
        print("Status: [] Error... Empty Input")
        Flame_Translate_FAULTY = []
    return Flame_Translate_CORRECT, Flame_Translate_FAULTY #Output Object Type: [List, List]


#Central Quantificiation function that will quantify the combinations of exons.
def QUANTIFYFUNC(QUANTIFYINPUT):
    #-----------Outside counters and variables-----------#
    Flame_Quantify_PROGRESSMAX = sum(1 for Flame_Quantify_LINE in QUANTIFYINPUT)
    Flame_Quantify_PROGRESSCOUNT1 = 0
    Flame_Quantify_PROGRESSCOUNT2 = 0
    Flame_Quantify_DICT = {}
    #-----------The Function Itself-----------#
    if QUANTIFYINPUT != []:
        for Flame_Quantify_COUNT1 in QUANTIFYINPUT: #Loop: Quantifying the Splice Permutaitons by going through a dictionary and  looking for new keys. If there is no matching keys, create a new splice permutation combination. If the are matching keys, just add to it's count.
            if Flame_Quantify_COUNT1 in Flame_Quantify_DICT.keys():
                Flame_Quantify_DICT[Flame_Quantify_COUNT1] += 1
            if Flame_Quantify_COUNT1 not in Flame_Quantify_DICT.keys():
                Flame_Quantify_DICT[Flame_Quantify_COUNT1] = 1
            #-----------The Progress Bar-----------#
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
            #-----------The Progress Bar-----------#
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
    for Flame_Emptyadjmtx_COUNT1 in range(Flame_Emptyadjmtx_REF1): #Simply create an empty adjacency matrix based on the number of elements as the input. This part is responsible for the correct number of columns.
        Flame_Emptyadjmtx_MATRIX.append([0]*Flame_Emptyadjmtx_REF1) #This part is responsible for the correct number of rows.
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    return Flame_Emptyadjmtx_MATRIX #Output Object Type: Nested List.


#ADD COMMENTS!
#Central Adjancency Matrix filling function that will count the connections.
def CORRECTADJMTXFUNC(CORRECTADJMTXINPUT, REF, EMPTYADJMTX):
    #-----------Outside counters and variables-----------#
    Flame_Rightadjmtx_PROGRESSMAX = sum(1 for Flame_Rightadjmtx_LINE in CORRECTADJMTXINPUT)
    Flame_Rightadjmtx_PROGRESSCOUNT1 = 0
    Flame_Rightadjmtx_PROGRESSCOUNT2 = 0
    Flame_Rightadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Rightadjmtx_EXONNR = 0
    Flame_Rightadjmtx_Counter1 = 0
    Flame_Rightadjmtx_Counter2 = 0 #Fix so that the matrix isn't symetric?
    Flame_Rightadjmtx_SPLICECOMB1 = ""
    Flame_Rightadjmtx_SPLICECOMB2 = ""
    #-----------The Function Itself-----------#
    if CORRECTADJMTXINPUT != []:
        for Flame_Rightadjmtx_COUNT1 in CORRECTADJMTXINPUT: #Loop: Iterates through the input (Translated Correct Reads).
            Flame_Rightadjmtx_EXONNR = (len(Flame_Rightadjmtx_COUNT1.split("-")[:-1])) #If one changes the "TMPSTRING" function, then one can avoid the use of this "Remove last element function"
            for Flame_Rightadjmtx_COUNT2 in range(Flame_Rightadjmtx_EXONNR):
                Flame_Rightadjmtx_SPLICECOMB1 = str(Flame_Rightadjmtx_COUNT1.split("-")[Flame_Rightadjmtx_COUNT2])
                Flame_Rightadjmtx_SPLICECOMB2 = str(Flame_Rightadjmtx_COUNT1.split("-")[Flame_Rightadjmtx_COUNT2+1])
                Flame_Rightadjmtx_Counter1 = 0
                Flame_Rightadjmtx_Counter2 = 0
                for Flame_Rightadjmtx_COUNT2 in REF: #Loop: Crosschecks with the reference. If it matches, breaks out of the loop, otherwise, add to a counter to know which position in the Adjacency matrix it needs to add to.
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
            #-----------The Progress Bar-----------#
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
            #-----------The Progress Bar-----------#
    elif CORRECTADJMTXINPUT == []:
        print("Status: [] Error... Empty Input")
    return Flame_Rightadjmtx_EMPTYADJMTX #Output Object Type: Nested List. 

#ADD COMMENTS!
def FAULTYADJMTXFUNC(FAULTYADJMTXINPUT, REF, EMPTYADJMTX, ADJMTXRANGESIZE):
    #-----------Outside counters and variables-----------#
    Flame_Wrongadjmtx_PROGRESSMAX = sum(1 for Flame_Wrongadjmtx_LINE in FAULTYADJMTXINPUT)
    Flame_Wrongadjmtx_PROGRESSCOUNT1 = 0
    Flame_Wrongadjmtx_PROGRESSCOUNT2 = 0
    Flame_Wrongadjmtx_EMPTYADJMTX = EMPTYADJMTX
    Flame_Wrongadjmtx_Counter1 = 0
    Flame_Wrongadjmtx_Counter2 = 0
    tmpcounter = 0
    Flame_Wrongadjmtx_RANGE1 = []
    Flame_Wrongadjmtx_RANGE2 = []
    #-----------The Function Itself-----------#
    for Flame_Wrongadjmtx_COUNT1 in FAULTYADJMTXINPUT: #Loop: Iterates through the input (Potential New Exons)
        Flame_Wrongadjmtx_Counter1 = 0
        Flame_Wrongadjmtx_Counter2 = 0
        for Flame_Wrongadjmtx_COUNT2 in REF: #Loop: Crosschecks with the Reference (Splice Candidates (X > 1% frequency). This loop is for the Row position.
            if (int(Flame_Wrongadjmtx_COUNT1.split("-")[0])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[0])+ADJMTXRANGESIZE):
                break
            elif (not (int(Flame_Wrongadjmtx_COUNT1.split("-")[0])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[0])+ADJMTXRANGESIZE)):
                Flame_Wrongadjmtx_Counter1 += 1
        for Flame_Wrongadjmtx_COUNT2 in REF: #Loop: Crosschecks with the Reference (Splice Candidates (X > 1% frequency). This loop is for the Col position.
            if (int(Flame_Wrongadjmtx_COUNT1.split("-")[1])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[1])+ADJMTXRANGESIZE):
                break
            elif (not (int(Flame_Wrongadjmtx_COUNT1.split("-")[1])-ADJMTXRANGESIZE <= Flame_Wrongadjmtx_COUNT2[0] <= int(Flame_Wrongadjmtx_COUNT1.split("-")[1])+ADJMTXRANGESIZE)):
                Flame_Wrongadjmtx_Counter2 += 1         
        Flame_Wrongadjmtx_EMPTYADJMTX[Flame_Wrongadjmtx_Counter2][Flame_Wrongadjmtx_Counter1] += 1 #Is the two "dimensions" in correct order? Is the columns representing the start?
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    return Flame_Wrongadjmtx_EMPTYADJMTX #Output Object Type: Nested List.

#ADD COMMENTS!
def FAULTYSEPERATORFUNC(FAULTYSEPINPUT, REF):
    #-----------Outside counters and variables-----------#
    Flame_Faultysep_PROGRESSMAX = sum(1 for Flame_Faultysep_LINE in FAULTYSEPINPUT)
    Flame_Faultysep_PROGRESSCOUNT1 = 0
    Flame_Faultysep_PROGRESSCOUNT2 = 0
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
        for Flame_Faultysep_COUNT3 in Flame_Faultysep_STRT_ALL: #Loop: Only if Flame_Faultysep_COUNT3 does not match, append the Exon range to a list that will be returned.
            if Flame_Faultysep_COUNT3 in Flame_Faultysep_NAMESEP:
                pass
            elif Flame_Faultysep_COUNT3 not in Flame_Faultysep_NAMESEP:
                Flame_Faultysep_POTENTIAL_EXON.append(Flame_Faultysep_COUNT3)
        #-----------The Progress Bar-----------#
        Flame_Faultysep_PROGRESSCOUNT1 += 1
        Flame_Faultysep_PROGRESSCOUNT2 += 1
        if Flame_Faultysep_PROGRESSCOUNT2 != Flame_Faultysep_PROGRESSMAX:
            if Flame_Faultysep_PROGRESSCOUNT1 >= int(round(Flame_Faultysep_PROGRESSMAX*0.01, 2)):
                PROGRESSBAR(Flame_Faultysep_PROGRESSCOUNT2 / Flame_Faultysep_PROGRESSMAX)
                Flame_Faultysep_PROGRESSCOUNT1 = 0
                time.sleep(0.001)
        elif Flame_Faultysep_PROGRESSCOUNT2 == Flame_Faultysep_PROGRESSMAX:
            PROGRESSBAR(Flame_Faultysep_PROGRESSCOUNT2/Flame_Faultysep_PROGRESSMAX)
            time.sleep(0.001)
        #-----------The Progress Bar-----------#
    return Flame_Faultysep_POTENTIAL_EXON #Output Object Type: List


#ADD COMMENTS!
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
    for Flame_FrequencySite_COUNT1 in REF:  
        if REF[Flame_FrequencySite_Counter1][3] > Flame_FrequencySite_MAX:
            Flame_FrequencySite_MAX = REF[Flame_FrequencySite_Counter1][3]
            Flame_FrequencySite_Counter1 += 1
        elif REF[Flame_FrequencySite_Counter1][3] <= Flame_FrequencySite_MAX:
            Flame_FrequencySite_Counter1 += 1
    for Flame_FrequencySite_COUNT2 in range((Flame_FrequencySite_MAX) +
                                            FREQUENCYSITERANGESIZE): #Create a smaller list for only the gene range?
        Flame_FrequencySite_GENEEMULATION.append(0) #Make X.append([0]*Y) function instead?
    for Flame_FrequencySite_COUNT2 in FREQUENCYSITEINPUT:
        Flame_FrequencySite_COUNT1RANGE = Flame_FrequencySite_COUNT2.split("-") #This splits the two numbers (Start and End). This makes it so one cannot see the connections.
        if int(Flame_FrequencySite_COUNT1RANGE[1]) < ((Flame_FrequencySite_MAX) +
                                                      FREQUENCYSITERANGESIZE):
            Flame_FrequencySite_START = list(range(int(Flame_FrequencySite_COUNT1RANGE[0])-FREQUENCYSITEWINDOWSIZE,
                                                   int(Flame_FrequencySite_COUNT1RANGE[0])+FREQUENCYSITEWINDOWSIZE)) #Add "+1" as they only have 4 values instead of 5 values (X +- 2)
            Flame_FrequencySite_END = list(range(int(Flame_FrequencySite_COUNT1RANGE[1])-FREQUENCYSITEWINDOWSIZE,
                                                 int(Flame_FrequencySite_COUNT1RANGE[1])+FREQUENCYSITEWINDOWSIZE))
            for Flame_FrequencySite_COUNT3 in Flame_FrequencySite_START:
                Flame_FrequencySite_GENEEMULATION[Flame_FrequencySite_COUNT3] += 1
            for Flame_FrequencySite_COUNT3 in Flame_FrequencySite_END:
                Flame_FrequencySite_GENEEMULATION[Flame_FrequencySite_COUNT3] += 1
        else:
            pass
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    return Flame_FrequencySite_GENEEMULATION #Output Object Type: List

#ADD COMMENTS!
def FREQUENCYTHRESHFUNC(FREQUENCYTRESHINPUT, PERCENTTHRESH, REF):
    #-----------Outside counters and variables-----------#
    Flame_FrequencyThresh_PROGRESSMAX = sum(1 for Flame_FrequencyThresh_LINE in FREQUENCYTRESHINPUT)
    Flame_FrequencyThresh_PROGRESSCOUNT1 = 0
    Flame_FrequencyThresh_PROGRESSCOUNT2 = 0
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
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    return Flame_FrequencyThresh_SPLICECANDIDATES #Output Object Type: List

#ADD COMMENTS!
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
    for Flame_SpliceSignal_COUNT1 in Flame_SpliceSignal_SPLICESIGNALINPUT:
        Flame_SpliceSignal_GU = False
        Flame_SpliceSignal_AG = False
        Flame_SpliceSignal_WINDOW = REF[(int(Flame_SpliceSignal_COUNT1[0])-3): #Make window customizeable?)
                                        (int(Flame_SpliceSignal_COUNT1[0])+4)]
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
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    return Flame_SpliceSignal_SPLICESIGNALINPUT #Output Object Type: List

#ADD COMMENTS!
def SHORTREADFUNC(SHORTREADINPUT1, SHORTREADINPUT2, SHORTREADCANDIDATES, REF):
    #-----------Outside counters and variables-----------#
    Flame_ShortRead_PROGRESSMAX = sum(1 for Flame_ShortRead_LINE in SHORTREADINPUT2)
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
    FLAMEASDF = 0
    #-----------The Function Itself-----------#
    for Flame_ShortRead_COUNT1 in SHORTREADINPUT1:
        Flame_ShortRead_START = int(str(Flame_ShortRead_COUNT1).split("\t")[3])
        Flame_ShortRead_CIGAR = str(Flame_ShortRead_COUNT1).split("\t")[5]
        Flame_ShortRead_CHROM = str(Flame_ShortRead_COUNT1).split("\t")[2]
        if "N" in Flame_ShortRead_CIGAR:
            Fĺame_ShortRead_LENGTH = sum(list(map(int, re.findall(r'\d+',
                                                                  Flame_ShortRead_CIGAR))))
            Flame_ShortRead_COMB = [Flame_ShortRead_START,
                                    Flame_ShortRead_CIGAR,
                                    (Flame_ShortRead_START+Fĺame_ShortRead_LENGTH)]
            if (REF[0][1]) < Flame_ShortRead_COMB[0] < (REF[-1][3]):
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
                else:
                    print("Error, CIGAR-String:", Flame_ShortRead_COMB)
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
                        Flame_ShortRead_CIGARCOUNT = int(Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1-1]) #The Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1-1] is not a mistake.
                        Flame_ShortRead_CIGAROPERATOR = Flame_ShortRead_ITEMS[Flame_ShortRead_COUNT1]
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
        #-----------The Progress Bar-----------#
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
        #-----------The Progress Bar-----------#
    #The function to print out the results from previous funciton. WRITE MORE
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
