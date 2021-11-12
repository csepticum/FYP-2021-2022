import csv

def read_csv(csvfilename):
    rows = []
    with open(csvfilename) as csvfile:
        file_reader = csv.reader(csvfile)
        for row in file_reader:
            rows.append(row)
    return rows

#this arrangemerge just outputs a sorted and merged SLC-H1 library with the y values reset for my plot to be used in R
def arrangemerge(csvv):
    with open("/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-16_MAximum Coverage MAE Clones/nodup_201014merge.csv", mode = 'w') as myfile:
        writer = csv.writer(myfile)
        mergedSLC = read_csv(csvv)
        writer.writerow(mergedSLC[0])
        mergedSLC = mergedSLC[1:]
        mergedSLC.sort(key = lambda x: int(x[5]))
        resety = {}
        y = 1 
        for clones in mergedSLC:
            if clones[15] not in resety.keys():
                resety[clones[15]] = y
                y += 1
        for clones in mergedSLC:
            clones[15] = resety[clones[15]]
            writer.writerow(clones)

toArrange = "/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-16_MAximum Coverage MAE Clones/mergedwith201014.csv"

########
# Now start to get the maximum coverage clones

####
# addendum
# 2 August 2021
# I made this sylMAE3 function to make an R plot which shows my selected clones
# from there, I will visually pick a few more clones out to make the largest span

file1 = "/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-16_MAximum Coverage MAE Clones/201014.csv"
file2 = "/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-16_MAximum Coverage MAE Clones/LsortedMergedSLCH1.csv"

# this function adds a length variable to all the rows in the merged SLCH1 library 
def length(file):
    openfile = read_csv(file)
    with open("/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-16_MAximum Coverage MAE Clones/LsortedMergedSLCH1.csv", mode = 'w') as myfile:
        writer = csv.writer(myfile)
        firstrow = openfile[0]
        firstrow.append("Length of Fragment")
        writer.writerow(firstrow)
        for clones in openfile[1:]:
            clones.append(int(clones[6]) - int(clones[5]) + 1)
            writer.writerow(clones)

#this is for checking if a particular clone in query exists as a subset of any other clones in the list of clones 
# if it is a subset, we return False
# if it is not a subset, we return True 
def subset(clone):
    start = int(clone[5])
    stop = int(clone[6])
    allClones = read_csv(file2)[1:]
    allClones.remove(clone)
    # print("length of the clone set", len(allClones))
    for clones in allClones:
        #return False if start >= int(clones[5]) and stop <= int(clones[6])
        if start >= int(clones[5]) and stop <= int(clones[6]):
            # if clone[0] == 'LOY7-23-4':
            #     print(clones)    
            return False
    return True


def sylMAE3(csvfile):
    allClones = read_csv(csvfile)
    firstrow = allClones[0]
    allClones = allClones[1:]
    allClones2 = allClones[1:]
    desiredClones = []
    desiredClones2 = []
    tempdict = {}
    y = 1
    with open("/Users/sylvester/Dropbox/NUS Mod Stuff/FYP/Programming @ GIS/SLComp1-16_MAximum Coverage MAE Clones/MaxCoverageClones.csv", 'w') as myfile:
        writer = csv.writer(myfile)
        writer.writerow(firstrow)
        while allClones: #keeps iterating as long as we still have stuff in the list
            largest = max(allClones, key = lambda x : int(x[17]))
            # print(largest)
            indextopop = allClones.index(largest)
            popout = allClones.pop(indextopop)
            #print(len(allClones))
            if subset(popout):
                #writer.writerow(popout)
                desiredClones.append(popout)
        # lines 97 - 108 accounts for when multiple fragments are excluded 
        # This can occur because the coordinates are stored per row in the csv spreadsheet
        # it is possible that a shorter fragment is excluded from the csv output but in fact that clone encodes multiple fragments
        for clones in desiredClones:
            for clone in allClones2:
                if clones[0] == clone[0]:
                    desiredClones2.append(clone)
        desiredClones2.sort(key = lambda x : x[5])
        for clones in desiredClones2:
            if clones[0] not in tempdict.keys():
                tempdict[clones[0]] = y
                y += 1
        for clones in desiredClones2:
            clones[15] = tempdict[clones[0]]
            writer.writerow(clones)
            

print(sylMAE3(file2))
    
    
    
