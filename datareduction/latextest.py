#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

x=np.arange(0,1,0.01)
y=np.cos(x)
plt.rc('text',usetex=True)
plt.plot(x,y)
plt.xlabel(r'$\rm \mu_{\delta}~(mas~yr^{-1})$')
plt.show()
