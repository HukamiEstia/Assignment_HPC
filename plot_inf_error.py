import matplotlib.pyplot as plt
import numpy
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("method", help="select the solving method")
args = parser.parse_args()

with open(f"{args.method}.dat") as f:
    data = f.read()
    dataTable = [[float(value) for value in line[:-1].split(' ')] for line in data[:-1].split('\n')]

with open("Analytical.dat") as f:
    data = f.read()
    dataTableAnalytical = [[float(value) for value in line[:-1].split(' ')] for line in data[:-1].split('\n')]


t = [i*0.05 for i in range(len(dataTable))]
infError = [sum([abs(dataTable[i][j] - dataTableAnalytical[i][j]) for j in range(len(dataTable[0]))]) for i in range(len(dataTable))]
plt.rcParams.update({'font.size': 18})

print("infError:", infError[49])

plt.figure()
plt.suptitle(f"{args.method} infinity error")
plt.plot(t, infError, c='r', label='time')
plt.xlabel('t')
plt.ylabel('E')
# # ax[i,j].set_xlim(0, 31)

# plt.subplots_adjust(top=0.85, bottom=0.08,
#                     left=0.10, right=0.95,
#                     hspace=0.60, wspace=0.35)

plt.show()
