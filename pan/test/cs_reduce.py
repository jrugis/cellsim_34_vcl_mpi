import numpy as np
import os
import struct
import sys

##################################################################
# functions
##################################################################
def read_write_data(fname):
  f = open(fname + ".bin", "rb")
  rows = struct.unpack('l', f.read(8))[0]
  cols = struct.unpack('l', f.read(8))[0]
  v = [212, 1796, 6465] # apical, basal, middle (in numerical order!)
  data = np.zeros((np.shape(v)[0], cols), dtype=np.float32)
  f.seek(v[0] * 4, 1)
  for i in range(0, cols):     # the data is in column order
    data[0, i] = struct.unpack('f', f.read(4))[0] # apical
    f.seek(((v[1] - v[0]) - 1) * 4, 1)
    data[2, i] = struct.unpack('f', f.read(4))[0] # basal
    f.seek(((v[2] - v[1]) - 1) * 4, 1)
    data[1, i] = struct.unpack('f', f.read(4))[0] # middle
    f.seek(((rows + v[0] - v[2]) - 1) * 4, 1)
  f.close()
  f = open(fname + "R.bin", "wb")
  f.write(struct.pack('l', np.shape(v)[0]))
  f.write(struct.pack('l', cols))
  rows = np.shape(v)[0]
  for i in range(0, cols):     # the data is in column order
    for j in range(0, rows):
      f.write(struct.pack('f', data[j, i]))
  f.close() 
  return

##################################################################
# main program
##################################################################

fname = "c"
read_write_data(fname)
fname = "ip3"
read_write_data(fname)
