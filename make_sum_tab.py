f_name = '850full.C'

out_name = "./summaryTables/"+f_name.split(".")[0]+".txt"
print(out_name)
f = open(f_name,'r')
nbrves = 0

sum_tab = []
ves_counter = 0
for line in f:
    line = line.strip(" ")
    if line[0] == "A":
        line = line.strip("Arteries[")
        line = line.split("] = new Tube( ")
        line[1] = line[1].split(",")
        line[0] = int(line[0])

        lst = [[line[0]], line[1][0:9:1]]

        row = []
        
        for item in lst:
            for sub_item in item:
                # print(type(sub_item) == str)
                if type(sub_item) == str:
                    sub_item = sub_item.strip(" ")
                    if sub_item[0].isalpha():
                        if sub_item == "gridPoints":
                            sub_item = 8
                        elif sub_item == "rm":
                            sub_item = 5
                        elif sub_item[0] == "A":
                            sub_item = int(sub_item.strip("Arteries[").strip("]"))
                    else:
                        sub_item = float(sub_item)
                row.append(sub_item)

        # print(row)

        sum_tab.append(row)
        # print(lst)
        # [line[0]].append(line[1])
        # print(line)
        if nbrves < line[0]:
            nbrves = line[0]

        # print(line)

        ves_counter += 1
        if ves_counter > nbrves:
            break
# print(sum_tab)
with open(out_name, 'w') as f:
    for row in sum_tab:
        rw = str(row).strip('[').strip(']')
        rw += "\n"
        print(rw)
        f.write(rw)

