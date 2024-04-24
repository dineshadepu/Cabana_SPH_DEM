# trace generated using paraview version 5.12.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os


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
files = sorted_files
# =====================================
# end: get the files and sort
# =====================================

# =====================================
# start: get the files and sort
# =====================================
files_rb = [filename for filename in os.listdir('.') if filename.startswith("rigid_bodies") and filename.endswith("xmf") ]
files_rb.sort()
# print(files_rb)
files_rb_num = []
for f in files_rb:
    f_last = f[13:]
    files_rb_num.append(int(f_last[:-4]))
files_rb_num.sort()

sorted_files_rb = []
for num in files_rb_num:
    sorted_files_rb.append("rigid_bodies_" + str(num) + ".xmf")
# print(sorted_files_rb)
files_rb = sorted_files_rb
# =====================================
# end: get the files and sort
# =====================================
# print(files)
# print(files_rb)

# create a new 'XDMF Reader'
rigid_bodies_0xmf = XDMFReader(registrationName='rigid_bodies_0.xmf*', FileNames=files_rb)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# set active source
SetActiveSource(rigid_bodies_0xmf)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
rigid_bodies_0xmfDisplay = Show(rigid_bodies_0xmf, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
rigid_bodies_0xmfDisplay.Representation = 'Surface'

# show color bar/color legend
rigid_bodies_0xmfDisplay.SetScalarBarVisibility(renderView1, True)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.15, 4.855, 0.0]
renderView1.CameraFocalPoint = [1.15, 0.5, 0.0]
renderView1.CameraViewUp = [1.0, 0.0, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# get color transfer function/color map for 'rb_ids'
rb_idsLUT = GetColorTransferFunction('rb_ids')

# get opacity transfer function/opacity map for 'rb_ids'
rb_idsPWF = GetOpacityTransferFunction('rb_ids')

# get 2D transfer function for 'rb_ids'
rb_idsTF2D = GetTransferFunction2D('rb_ids')

# change representation type
rigid_bodies_0xmfDisplay.SetRepresentationType('Point Gaussian')

# set scalar coloring
ColorBy(rigid_bodies_0xmfDisplay, ('POINTS', 'rb_force', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(rb_idsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
rigid_bodies_0xmfDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
rigid_bodies_0xmfDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'rb_force'
rb_forceLUT = GetColorTransferFunction('rb_force')

# get opacity transfer function/opacity map for 'rb_force'
rb_forcePWF = GetOpacityTransferFunction('rb_force')

# get 2D transfer function for 'rb_force'
rb_forceTF2D = GetTransferFunction2D('rb_force')

# create a new 'XDMF Reader'
particles_0xmf = XDMFReader(registrationName='particles_0.xmf*', FileNames=files)

# set active source
SetActiveSource(particles_0xmf)

# show data in view
particles_0xmfDisplay = Show(particles_0xmf, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
particles_0xmfDisplay.Representation = 'Surface'

# show color bar/color legend
particles_0xmfDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'ids'
idsLUT = GetColorTransferFunction('ids')

# get opacity transfer function/opacity map for 'ids'
idsPWF = GetOpacityTransferFunction('ids')

# get 2D transfer function for 'ids'
idsTF2D = GetTransferFunction2D('ids')

# change representation type
particles_0xmfDisplay.SetRepresentationType('Point Gaussian')

# set scalar coloring
ColorBy(particles_0xmfDisplay, ('POINTS', 'frc_dem', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(idsLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
particles_0xmfDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
particles_0xmfDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'frc_dem'
frc_demLUT = GetColorTransferFunction('frc_dem')

# get opacity transfer function/opacity map for 'frc_dem'
frc_demPWF = GetOpacityTransferFunction('frc_dem')

# get 2D transfer function for 'frc_dem'
frc_demTF2D = GetTransferFunction2D('frc_dem')

renderView1.ResetActiveCameraToNegativeZ()

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1542, 751)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.15, 0.5, 4.845059295778354]
renderView1.CameraFocalPoint = [1.15, 0.5, 0.0]
renderView1.CameraParallelScale = 1.2539936203984452


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------
