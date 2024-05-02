import h5py
import numpy as np
import matplotlib.pyplot as plt
import os


# start: Get the files
files = [filename for filename in os.listdir('.') if filename.startswith("particles") and filename.endswith("xmf") ]
files.sort()
files_num = []
for f in files:
    f_last = f[10:]
    files_num.append(int(f_last[:-4]))
files_num.sort()

sorted_files = []
for num in files_num:
    sorted_files.append("particles_" + str(num) + ".h5")
print(sorted_files)
# ends: Get the files
files = sorted_files

contact_count = []
for f in files:
    f = h5py.File(f, "r")
    # print(np.array(f["radius"]))
    x = np.array(f["rb_position"][:, 0])
    y = np.array(f["rb_position"][:, 1])
    contact_count.append(f["rb_contact_count"][0])

    cnt = f["rb_contact_count"][0]
    if cnt > 0:
        # print(f["rb_contact_idx"][1])
        print("----------------------------")
        print(f["rb_contact_tng_frc"][0])
        print(f["rb_contact_tng_frc"][1])
        print("----------------------------")

# print(x)
# print(y)
print(contact_count)

# plt.scatter(x, y, label="SPH appr")
# plt.legend()
# plt.savefig("colliding_fluid_blocks.png")
# plt.show()
# plt.plot()
