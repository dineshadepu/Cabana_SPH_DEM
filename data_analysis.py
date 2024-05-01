import h5py
import numpy as np
import matplotlib.pyplot as plt

f = h5py.File("build/particles_10.h5", "r")
# print(np.array(f["radius"]))
x = np.array(f["rb_position"][:, 0])
y = np.array(f["rb_position"][:, 1])
print(x)
print(y)
# plt.scatter(x, y, label="SPH appr")
# plt.legend()
# plt.savefig("colliding_fluid_blocks.png")
# plt.show()
# plt.plot()
