from os import listdir
from os.path import isfile, join
import os
import shutil
import pandas
f = "./NewSumTab/full.tsv"


def c_file(f_name):
    f = open(f_name,'r')
    nbrves = 0

    sum_tab = []
    ves_counter =0

    s = []
    t = []
    for line in f:

        if "A" in line:
            line = line.strip(" ")
            line = line.strip('\n')
            line = line.split(',')
            line[0] = line[0].split("=")
            line[0] = line[0][0]

            # print(line[0:5:1])
            for i in [0, 3, 4, 5]:
                line[i] = line[i].strip("// ")
                line[i] = line[i].strip("Arteries[")
                line[i] = line[i].strip("]")
                # print(line[i])

            for i in [3, 4, 5]:
                if line[i] != "0":
                    print(line[0], line[i])
                    s.append(int(line[0]))
                    t.append(int(line[i]))
    return s, t


def sum_tab(f_name):
    data = open(f_name,'r')
    s = []
    t = []

    for line in data:
        line = line.strip("\n")
        line = line.split("\t")
        print(line[0], line[4], line[5], line[6])

        for i in [4, 5, 6]:
            if line[i]!= "0":
                s.append(int(line[0]))
                t.append(int(line[i]))

    return s, t


s, t = sum_tab(f)
print("s = "+str(s)+ ";")
print("t = "+ str(t)+";")