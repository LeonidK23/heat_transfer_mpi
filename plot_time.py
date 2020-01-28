import numpy as np
import matplotlib.pyplot as plt

with open("computation_pool.txt") as f:
    data = f.readlines()

data = [line.strip() for line in data if line.strip() != "ready"]
sizes = []
times = []
for i in range(0, len(data), 2):
    sizes.append(int(data[i].split(" ")[-1]))
    times.append(float(data[i+1].split(" ")[-2]))

plt.plot(sizes, times)
plt.xticks(np.arange(min(sizes), max(sizes)+1, 2.0))
plt.title("Time ~ ghost zone size")
plt.xlabel("Ghost zone size")
plt.ylabel("Computation time")
plt.show()
