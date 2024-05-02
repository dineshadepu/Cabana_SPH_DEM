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
    idx_check = 1
    contact_count.append(f["rb_contact_count"][idx_check])

    cnt = f["rb_contact_count"][idx_check]
    if cnt > 0:
        # print(f["rb_contact_idx"][1])
        print("----------------------------")
        # print(f["rb_contact_tng_frc"][0])
        # print(f["rb_contact_tng_frc"][1])
        print(f["rb_contact_idx"][idx_check])
        # print(f["rb_contact_tng_frc"][idx_check])
        print(f["rb_contact_fn_magn"][idx_check])
        print("----------------------------")

# print(x)
# print(y)
print(contact_count)

# plt.scatter(x, y, label="SPH appr")
# plt.legend()
# plt.savefig("colliding_fluid_blocks.png")
# plt.show()
# plt.plot()

# // auto rb_position = Cabana::slice<0>     ( rb,    "rb_position");
# // auto rb_limits = Cabana::slice<1>       ( rb,      "rb_limits");
# // auto rb_ids = Cabana::slice<2>          ( rb,         "rb_ids");
# // auto rb_velocity = Cabana::slice<3>     ( rb,    "rb_velocity");
# // auto rb_force = Cabana::slice<4>        ( rb,       "rb_force");
# // auto rb_torque = Cabana::slice<5>       ( rb,      "rb_torque");
# // auto rb_lin_acc = Cabana::slice<6>      ( rb,     "rb_lin_acc");
# // auto rb_ang_acc = Cabana::slice<7>      ( rb,     "rb_ang_acc");
# // auto rb_ang_mom = Cabana::slice<8>      ( rb,     "rb_ang_mom");
# // auto rb_ang_vel = Cabana::slice<9>      ( rb,     "rb_ang_vel");
# // auto rb_rotation_matrix = Cabana::slice<10>     ( rb,     "rb_rotation_matrix");
# // auto rb_mass = Cabana::slice<11>        ( rb,        "rb_mass");
# // auto rb_density = Cabana::slice<12>     ( rb,     "rb_density");
# // auto rb_body_moi = Cabana::slice<13>    ( rb,    "rb_body_moi");
# // auto rb_inv_body_moi = Cabana::slice<14>( rb,    "rb_inv_body_moi");
# // auto rb_global_moi = Cabana::slice<15>  ( rb,    "rb_global_moi");
# // auto rb_inv_global_moi = Cabana::slice<16>( rb,    "rb_inv_global_moi");
# // auto rb_rotation_angle = Cabana::slice<17>  ( rb,    "rb_rotation_angle");
# // auto rb_I_zz = Cabana::slice<18>( rb,    "rb_I_zz");
# // auto rb_rad_s = Cabana::slice<19>( rb,    "rb_rad_s");
# // auto rb_moi_1 = Cabana::slice<20>( rb,    "rb_moi_1");
# // auto rb_E = Cabana::slice<21>( rb,    "rb_E");
# // auto rb_nu = Cabana::slice<22>( rb,    "rb_nu");
# // auto rb_contacts_count = Cabana::slice<23>( rb,    "rb_contact_count");
# // auto rb_contact_idx = Cabana::slice<24>( rb,    "rb_contact_idx");
# // auto rb_contact_tng_frc = Cabana::slice<25>( rb,    "rb_contact_tng_frc");
# // auto rb_contact_tng_disp = Cabana::slice<26>( rb,    "rb_contact_tng_disp");
# // auto rb_contact_fn_magn = Cabana::slice<27>( rb,    "rb_contact_fn_magn");
# // auto rb_contact_normal_overlap = Cabana::slice<28>( rb,    "rb_contact_normal_overlap");
# // auto rb_G = Cabana::slice<29>( rb,    "rb_G");
