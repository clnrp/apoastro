#!/usr/bin/python
# -*- coding: UTF-8 -*-

# Author Cleoner S. Pietralonga
# e-mail: cleonerp@gmail.com
# Apache License

import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import *

from apoastro.keplerian_elements import *
from apoastro.utils import *

# Earth
K = KeplerianElements(1.0, 0.01673163, 0.00054346, 100.46691572, 102.93005885, 354.88739611)
K.setRates(-0.000704431, -0.00003661, -0.01337178, 35999.37306329, 0.31795260, -0.24123856)

position = K.getPosition(toJ2000(time.time()))
print(position)