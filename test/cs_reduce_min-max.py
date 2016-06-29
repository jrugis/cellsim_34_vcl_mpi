import numpy as np
import os
import struct
import sys

import cs

##################################################################
# functions
##################################################################
def read_write_data(fname):
  data = cs.get_data(fname + ".bin")
  sys.stdout.flush()
  max_per_row = np.amax(data, axis=1)
  rows = [np.argmax(max_per_row), np.argmin(max_per_row)]
  cols = data.shape[1]
  f = open(fname + "R.bin", "wb")
  f.write(struct.pack('l', len(rows)))
  f.write(struct.pack('l', cols))
  for i in range(0, cols):     # the data is in column order
    for j in rows:
      f.write(struct.pack('f', data[j, i]))
  f.close() 
  return

##################################################################
# main program
##################################################################

fname = "c"
print "reducing " + fname + ".bin"
read_write_data(fname)

fname = "ip3"
print "reducing " + fname + ".bin"
read_write_data(fname)

