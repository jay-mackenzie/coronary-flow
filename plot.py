import numpy as np
from matplotlib import pyplot as plt

# def

def name(path):
    return path+"/"+path+".2d"

def load(two_D_file):
    return np.loadtxt(two_D_file,delimiter=',')

print("choose name")
# pth = input()
pth = "test2"
nm = name(pth)
# print(nm)
dt = load(nm)

p = []
q = []
# print(dt)


for row in dt:
    if row[0] == 2:
        if row[2] == 0:
            p.append(row[3])
            q.append(row[4])
plt.plot(q)
plt.show()
