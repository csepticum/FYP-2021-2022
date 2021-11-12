import csv
#import numpy as np
#import matplotlib.pyplot as mpp
#import openpyxl
#from openpyxl import Workbook
#from openpyxl import load_workbook

def read_csv(csvfilename):
    """
    Reads a csv file and returns a list of list
    containing rows in the csv file and its entries.
    """
    rows = []
    with open(csvfilename, encoding='utf-8') as csvfile:
        file_reader = csv.reader(csvfile)
        for row in file_reader:
            rows.append(row)
    return rows

def genomesorter(genlist):
    try:
        outputlist = []
        for ele in genlist:
            if ele == genlist[-1]: #for when we reach final element in list
                outputlist.append([genlist[0],genlist[-1]])
                print("")
                print("final output list:", outputlist)
                return outputlist
            nextnum = genlist[genlist.index(ele) +1]
            if ele == nextnum - 1: #depicts continuous sequence
                print("still ongoing")
                continue
            else:
                outputlist.append([genlist[0], genlist[genlist.index(nextnum)-1]])
                #remove excluded continuous sequence
                genlist = genlist[genlist.index(nextnum):]
                print("current outputlist:", outputlist)
                print("remaining nucleotides to process:", len(genlist))
                continue
    except IndexError as e:
        print("End of list")
        return outputlist


def parser(filename):
    #genome = list(range(1, 4641652+1))
    data = read_csv(filename)
    datas = data[1:]
    #newfilename = filename + "_unhybridized_ranges" + ".xlsx"
    row = 2
    try:
        for clones in datas:
            currentrange = list(range(int(clones[5]),int(clones[6])+1))
            print("current checking nucleotide ranges:", currentrange[0], currentrange[-1])
            print("before filter", len(genome))
            genome1 = list(filter(lambda x: x < int(clones[5]), genome))
            genome2 = list(filter(lambda x: x > int(clones[6]), genome))
            tempgenomelength = len(genome)
            genome.clear()
            genome = genome1+genome2
            nucremoved = len(genome) - tempgenomelength
            print("nucleotides removed:", len(genome) - tempgenomelength)
            print("after filter", len(genome))
            print("")
    except Exception as e:
        print("error encountered, test terminated")
        ranger = genome
        return genomesorter(genome)

    finally:
        print("run ended")
        ranger = genome
        return genomesorter(genome)

def parser2(filename):
    genome = [[1,4641653]]
    data = read_csv(filename)
    datas = data[1:]
    row = 2
    try:
        for clones in datas: #iterate through hybrid list
            print("accessing excel row:,", row)
            print(clones[0])
            for fragments in genome:                
                fragment = list(range(fragments[0], fragments[1]+1))
                print("ongoing", genome[-1][-1])
                genome1 = list(filter(lambda x: x < int(clones[5]), fragment))
                genome2 = list(filter(lambda x: x > int(clones[6]), fragment))
                if len(genome1 + genome2) < len(fragment): #this means there is hybridization found
                    genome.remove(fragments)
                    if len(genome1) > 0:
                        genome.append([genome1[0], genome1[-1]])
                    if len(genome2) > 0:       
                        genome.append([genome2[0], genome2[-1]])
            row+=1
    except KeyboardInterrupt as e:
        print(genome)
    finally:
        print("run ended")
        unhybridized = 0
        for unhybs in genome:
            unhybridized += (int(unhybs[1] - int(unhybs[0]) + 1))
        print("unhybridized percentage =", ((unhybridized/4641652)*100), "%")
        print(genome)
        return genome


#parser2("file.csv")

