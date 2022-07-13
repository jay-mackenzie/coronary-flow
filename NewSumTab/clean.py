from os import listdir
from os.path import isfile, join
import os
import shutil


mypath = "./"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

done = ["coupled.txt", "700trif.txt", "adjustedfull.txt", "850bifs.txt"]


for n in onlyfiles:

    if "txt" in n and n not in done:
        print(n)

        c = 0
        s = ''
        f = open(n,'r')

        for line in f:
            

            # print(line)
            line = line.strip(" ")
            line = line.strip("\n")
            line = line.split(',')
            line = line[0:7:]
            # print(line)

            for i, _ in enumerate(line):
                if "." in line[i]:
                    # print(line[i])
                    line[i] = round(float(line[i]), 5)
                else:
                    line[i] = int(line[i])

                s += str(line[i]) + " \t"
            s += "\n"
            c += 1

        outname = n.split('.')[0] + '.tsv'
        with open(outname, 'w+') as f:
            f.write(s)
        
        shutil.copyfile(n, "./old/"+n)

        # break




