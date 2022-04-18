import matplotlib.pyplot as plt
import matplotlib
import numpy as np

n_file = 150544*4
for i in range(0, n_file, 200000):
    
    filename = "FR_P" + str(i) + "_process.dat"
    f = open(filename, 'r')
    raw_data = f.read()

    f.close()

    time = raw_data.split(" ")[0:2000:2]
    FR = raw_data.split(" ")[1:2000:2]
    time = list(map(float, time))
    FR = list(map(float, FR))
    
    plt.ylim(1, 60)
    plt.xlim(1000, 2000)
    plt.plot(FR, time, "-")
plt.savefig("wholenode.eps", dpi=300)
plt.savefig("wholenode.png", dpi=300)
plt.show()
