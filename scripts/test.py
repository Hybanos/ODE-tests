import matplotlib.pyplot as plt
import numpy as np

# data = np.fromfile("out.txt", dtype=float)
data = np.loadtxt("Euler.txt", delimiter=";").T

steps = data[0]
m1 = data[1]
p1 = data[2:5].T
v1 = data[5:8].T
m2 = data[8]
p2 = data[9:12].T
v2 = data[12:15].T

print(p1)

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

ax.plot(*p1.T)
ax.plot(*p2.T)

plt.show()