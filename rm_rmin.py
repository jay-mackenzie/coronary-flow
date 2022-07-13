import glob

files = glob.glob("./NewSumTab/*.tsv")

for f_name in files:
	print(f_name)
	f = open(f_name,'r')
	tmp = ""
	for line in f:
		line = line.strip("\n")
		line = line.split("\t")
		if len(line) == 9:
			for i in range(len(line)):
				line[i] = line[i].strip(" ")

			r = list(range(9))
			r.remove(7)
			
			for i in r:
				tmp += line[i] + "\t"
			
			tmp += "\n"
	f.close()

	# f = open(f_name, "w")

	# f.write(tmp)
	
	# f.close()

	print(tmp)
	


