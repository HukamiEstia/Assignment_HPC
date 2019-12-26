import matplotlib.pyplot as plt
import numpy
import argparse

# instanciate parser
parser = argparse.ArgumentParser()
parser.add_argument("method", help="select the solving method")

# parse arguments
args = parser.parse_args()

# read the data from file and convert to float
with open(f"{args.method}.dat") as f:
    data = f.read()
    dataTable = [[float(value) for value in line[:-1].split(' ')] for line in data[:-1].split('\n')]

# print(len(dataTable))
# print(dataTable[10])

# construct the x-axis
x = [i*0.05 for i in range(len(dataTable[0]))]

plt.rcParams.update({'font.size': 18})

# create the figure for plottings
fig, ax = plt.subplots(2,3)
fig.figsize=(16,4)
fig.suptitle(f"{args.method} solvings")
for i in range(len(ax)):
    for j in range(len(ax[0])):
        # plot specified timestep
        TimePoint = int(((len(dataTable) - 1)/5)*(i * len(ax[0]) + j))
        print(TimePoint)
        ax[i,j].plot(x, dataTable[TimePoint], c='r', label='time')
        ax[i,j].set_xlabel('x')
        ax[i,j].set_ylabel('T')
        ax[i,j].set_xlim(0, 31)
        ax[i,j].set_ylim(0, 300)
        ax[i,j].set_title(f"t = {round(0.1 * (i * len(ax[0]) + j), 1)}h")

plt.subplots_adjust(top=0.85, bottom=0.08,
                    left=0.10, right=0.95,
                    hspace=0.60, wspace=0.35)

plt.show()
