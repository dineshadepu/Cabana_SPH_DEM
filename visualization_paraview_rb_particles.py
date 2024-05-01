# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *

import os, sys
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# =====================================
# start: get the files and sort
# =====================================
files = [filename for filename in os.listdir('.') if filename.startswith("particles") and filename.endswith("xmf") ]
files.sort()
files_num = []
for f in files:
    f_last = f[10:]
    files_num.append(int(f_last[:-4]))
files_num.sort()

sorted_files = []
for num in files_num:
    sorted_files.append("particles_" + str(num) + ".xmf")
print(sorted_files)
files = sorted_files
# =====================================
# end: get the files and sort
# =====================================



# create a new 'XDMF Reader'
particles_0xmf = XDMFReader(registrationName='particles_0.xmf*', FileNames=files)
particles_0xmf.PointArrayStatus = [
		 "rb_position",
		 "rb_ids",
		 "rb_limits",
		 "rb_velocity",
		 "rb_force",
		 "rb_torque",
		 "rb_lin_acc",
		 "rb_ang_acc",
		 "rb_ang_mom",
		 "rb_ang_vel",
		 "rb_rotation_matrix",
		 "rb_mass",
		 "rb_density",
		 "rb_body_moi",
		 "rb_inv_body_moi",
		 "rb_global_moi",
		 "rb_inv_global_moi",
		 "rb_rotation_angle",
		 "rb_I_zz",
		 "rb_rad_s",
		 "rb_moi_1",
		 "rb_E",
		 "rb_nu",
		 "rb_contacts_count",
		 "rb_contact_idx",
		 "rb_contact_tng_frc",
		 "rb_contact_tng_disp",
		 "rb_contact_fn_magn",
		 "rb_contact_normal_overlap"]

particles_0xmf.GridStatus = ['points']

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
particles_0xmfDisplay = Show(particles_0xmf, renderView1, 'UnstructuredGridRepresentation')

particles_0xmfDisplay.SetRepresentationType('Point Gaussian')

particles_0xmfDisplay.GaussianRadius = 0.025

ColorBy(particles_0xmfDisplay, ('POINTS', 'rb_mass'))

# particles_0xmfDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
particles_0xmfDisplay.SetScalarBarVisibility(renderView1, True)
