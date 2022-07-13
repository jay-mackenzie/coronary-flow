import numpy as np
from matplotlib import pyplot as plt
# import matplotlib

# def

def name(path):
    return path+".2d"

def load(two_D_file):
    return np.loadtxt(two_D_file,delimiter=',')

print("choose name")
# pth = input()
pth = "./Outputs/test2_100.2d"
# nm = name(pth)
# print(nm)
dt = load(pth)

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
