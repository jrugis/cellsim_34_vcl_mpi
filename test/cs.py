#from mayavi import mlab
import numpy as np
import struct
import subprocess

##################################################################
# functions
##################################################################

##################################################################
# get the mesh node coordinates
def get_mesh_coords(fname):
  f1 = open(fname, 'r')
  for line in f1: 
    if line.startswith("$Nodes"): break
  for line in f1:
    pcount = int(line)
    break
  x = np.empty((pcount))
  y = np.empty((pcount))
  z = np.empty((pcount))
  t = 0
  for line in f1:
    v = line.split()
    x[t] = float(v[1])
    y[t] = float(v[2])
    z[t] = float(v[3])
    t += 1
    if t == pcount: break
  f1.close
  return x, y, z
##################################################################
# read in a bin data file
def get_data(fname):
  f1 = open(fname, "rb")
  rows = struct.unpack('l', f1.read(8))[0]
  cols = struct.unpack('l', f1.read(8))[0]
  data = np.zeros((rows, cols), dtype=np.float32)
  for j in range(0, cols):     # the data is in column order
    for i in range(0, rows):
      t = f1.read(4)
      data[i,j] = struct.unpack('f', t)[0]
      pass
  f1.close()
  return data
##################################################################
# read in part of a simulation data file
def get_part_data(fname):
  f1 = open(fname, "rb")
  rows = struct.unpack('l', f1.read(8))[0]
  cols = struct.unpack('l', f1.read(8))[0]
  v = [212, 1796, 6465] # apical, basal, middle (in numerical order!)
  data = np.zeros((np.shape(v)[0], cols), dtype=np.float32)
  f1.seek(v[0] * 4, 1)
  for i in range(0, cols):     # the data is in column order
    data[0, i] = struct.unpack('f', f1.read(4))[0] # apical
    f1.seek(((v[1] - v[0]) - 1) * 4, 1)
    data[2, i] = struct.unpack('f', f1.read(4))[0] # basal
    f1.seek(((v[2] - v[1]) - 1) * 4, 1)
    data[1, i] = struct.unpack('f', f1.read(4))[0] # middle
    f1.seek(((rows + v[0] - v[2]) - 1) * 4, 1)
  f1.close()
  return data
##################################################################
# get simulation time data
def get_sim_time():
  #fname = subprocess.check_output("ls *.dat", shell=True).split()[0]
  fname = "cs1.dat"
  f1 = open(fname, "r")
  for line in f1:
    if line.startswith("%! delT totalT"):
      v = next(f1).rstrip().split()
  f1.close()
  return float(v[0]), float(v[1])
##################################################################
# plot mesh with scalar value
def plot_mesh(x, y, z, v, wtitle, ptitle):
  mlab.figure(figure = wtitle)
  p = mlab.points3d(x, y, z, v, colormap='RdYlGn')
  p.module_manager.scalar_lut_manager.reverse_lut = True
  mlab.outline(p)
  mlab.colorbar(p, orientation='vertical')
  mlab.text(0.05, 0.9, ptitle)
##################################################################
##################################################################
##################################################################
##################################################################


