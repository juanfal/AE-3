#!/usr/bin/env python3
# encoding: utf-8
# normaldist.py
# juanfc 2019-05-03

import os
import numpy as np

mu, sigma = 0, 10 # mean and standard deviation

with open(os.path.expanduser("~/Desktop/salida.txt"), 'w') as outf:
    for i in range(10000):
        print(np.random.normal(mu, sigma), file=outf)
