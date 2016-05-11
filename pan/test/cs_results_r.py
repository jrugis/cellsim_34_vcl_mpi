
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import os
import struct
import sys

import cs

##################################################################
# main program
##################################################################

# get the simulation time scale
t_delta, t_end = cs.get_sim_time()
y = np.linspace(0.0, t_end, t_end/t_delta, endpoint=True)
print(len(y))

# plot the calcium and ip3 simulation data
plt.rcParams['axes.color_cycle'] = ['r', 'g', 'b']
fig, plots = plt.subplots(2, sharex=True)
plots[0].set_title('calcium')
plots[0].plot(y, np.transpose(cs.get_data('cR.bin')), lw=0.5)
plots[1].set_title('ip3')
plots[1].set_ylabel("concentration (uM)")
plots[1].plot(y, np.transpose(cs.get_data('ip3R.bin')), lw=0.5)
plt.xlabel('time (s)')

open('results.pdf', 'w').close()
plt.savefig('results.pdf')
plt.show()

