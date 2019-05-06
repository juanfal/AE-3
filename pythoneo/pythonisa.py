#!/usr/bin/env python3
# encoding: utf-8
# prueba.py
# juanfc 2019-02-07

# https://stackoverflow.com/questions/30399534/shift-elements-in-a-numpy-array#30534478


# class MyClass:
#     """A simple example class"""
#     i = 12345

#     def f(self):
#         return 'hello world'

# x = MyClass()

# x.f()

# import numpy as np
# def test():
#     a[0,0,0] = 34
#     print(a)

# a = np.zeros( (2, 3, 4), dtype=np.int)
# test()
# print(a)


# import numpy as np

# nitems = 100
# ncells = 3
# cells = np.zeros((ncells), dtype=np.int)
# for _ in range(nitems):
#     dest = np.random.randint(ncells)
#     cells[dest] += 1
# print(cells)



# import numpy as np
# nitems = 100
# ncells = 3
# values = np.random.randint(0,ncells,nitems)
# cells  = np.array_split(values,3)
# lengths= [ len(cell) for cell in cells ]
# print(lengths,np.sum(lengths))

# import numpy as np
# nitems = 100
# ncells = 3

# values = np.random.randint(0,ncells,nitems)
# # print(values, len(values), values.sum())
# ind_split = [ np.random.randint(0,nitems) ]
# print(ind_split)

# ind_split.append(np.random.randint(ind_split[-1],nitems))
# print(ind_split)

# cells  = np.array_split(values,ind_split)
# print(cells)
# # lengths= [ len(cell) for cell in cells ]
# # print(lengths)
# # print(lengths,np.sum(lengths))

# import numpy as np
# nitems = 800
# ncells = 60

# values = np.sort(np.random.randint(0, nitems+1, ncells-1, dtype=np.int))
# print(values)
# # TEST
# r = np.zeros(ncells, dtype=np.int)
# r[0] = values[0]
# r[-1] = nitems-values[-1]
# for i in range(1, ncells-1):
#     r[i] = values[i]-values[i-1]
# print(r)
# print(r.sum())


# import numpy as np
# nitems = 800
# ncells = 29

# values = np.sort(np.random.randint(0, nitems+1, ncells, dtype=np.int))
# values[-1] = nitems
# print(values)
# # TEST
# prev = 0
# for i in range(0, ncells):
#     prev, values[i] = values[i], values[i]-prev
# print(values)
# print(values.sum())

# from numpy import *
# nitems = 800
# ncells = 10

# values = sort(random.randint(0, nitems+1, ncells, dtype=int))
# values[-1] = nitems
# print(values)
# values = delete(insert(values, -1, nitems)-insert(values,0,0), -1)
# print(values)
# print(values.sum())





# # --------------------------------
# from numpy import *
# from time import sleep

# g_nitems = 1000
# g_ncells = 10
# g_nsamples = 100000


# # range_array = [np.random.randint(nitems) for _ in range(ncells - 1)] + [0, nitems]
# # range_array.sort()
# # cells = [range_array[i + 1] - range_array[i] for i in range(ncells)]
# # print(cells)


# def genDist(nitems, ncells):
#     r = sort(array([random.randint(nitems+1) for _ in range(ncells - 1)] + [0, nitems]))
#     cells = array([r[i + 1] - r[i] for i in range(ncells)])
#     return cells

# # Some stats

# test = zeros(g_ncells, dtype=int)
# Max = zeros(g_ncells, dtype=int)
# for _ in range(g_nsamples):
#     tmp = genDist(g_nitems, g_ncells)
#     # print(tmp.sum(), tmp, end='\r')
#     # sleep(0.5)
#     test += tmp
#     for i in range(g_ncells):
#         if tmp[i] > Max[i]:
#             Max[i] = tmp[i]

# print("\n", Max)
# print(test//g_nsamples)


# --------------------------------

# # https://gist.github.com/3608e423b2cbf218ace3d10f5d8b062c
# # https://gist.github.com/e8371d4b84dfb62a6fbf9b92be90436d
# # https://stackoverflow.com/questions/54773159/simple-random-distribution-of-n-items-on-n-cells/54882138#54882138
# from numpy import *
# from time import sleep

# g_nitems =    10000
# g_ncells =   10
# g_nsamples = 10000

# def genDist(nitems, ncells):
#     r = sort(random.randint(0, nitems+1, ncells-1))
#     return concatenate((r,[nitems])) - concatenate(([0],r))

# # Some stats

# test = zeros(g_ncells, dtype=int)
# Max = zeros(g_ncells, dtype=int)
# for _ in range(g_nsamples):
#     tmp = genDist(g_nitems, g_ncells)
#     # print(tmp.sum(), tmp, end='\r')
#     # print(_, end='\r')
#     # sleep(0.5)
#     test += tmp
#     for i in range(g_ncells):
#         if tmp[i] > Max[i]:
#             Max[i] = tmp[i]

# print("\n", Max)
# print(test//g_nsamples)

# # Enumerados, malamente:
# #
# def enum(*args):
#     enums = dict(zip(args, range(len(args))))
#     return type('Enum', (), enums)

# gForms = enum("INDIVIDUAL",       # they consume df
#               "RECIPIENT",       # they consume df (it doesn't consume directly)
#               "ACTOR",           # they consume df+if of recipient
#               "RECIPROCAL")      # they consume (intra 2x df+if) (inter df1+if1+df2+if2)


# print(gForms.INDIVIDUAL)
# print(gForms.__name__)
# print(gForms.__bases__)
# print(gForms.__dict__)





# # con Class
# from enum import IntEnum
# class Form(IntEnum):
#     INDIVIDUAL = 0      # they consume df
#     RECIPIENT = 1       # they consume df (it doesn't consume directly)
#     ACTOR = 2           # they consume df+if of recipient
#     RECIPROCAL = 3      # they consume (intra 2x df+if) (inter df1+if1+df2+if2)
# gNforms = 4



# from scipy.stats import norm
# print(norm.ppf(0.764401401, 10, 3))



# import re
# def is_number_regex(s):
#     """ Returns True is string is a number. """
#     if re.match("^\d+?\.\d+?$", s) is None:
#         return s.isdigit()
#     return True



# def isInt(anyNumberOrString):
#     try:
#         int(anyNumberOrString)
#         return True
#     except ValueError :
#         return False

# def nextCol(col):
#     return chr(ord(col)+1)

# def next_string(s):
#     strip_zs = s.rstrip('z')
#     if strip_zs:
#         return strip_zs[:-1] + chr(ord(strip_zs[-1]) + 1) + 'a' * (len(s) - len(strip_zs))
#     else:
#         return 'a' * (len(s) + 1)



# print(nextCol('a'))
# print(next_string('A1avezz'))


# from scipy.stats import norm

# def frange(x, y, jump):
#   while x < y:
#     yield x
#     x += jump

# def tukeyppf(prob, mean, stddev):
#     # return 0
#     return mean+((prob**0.14) - (1 - prob)**0.14)

# def ppf(prob, mean, stddev):
#     return 0
#     # return norm.ppf(prob, mean, stddev)


# for p in frange(-1.0, 1.0, 0.2):
#     print("%.1f: %4.1f  %4.1f" % (p, tukeyppf(p, 10, 5), ppf(p, 10, 5)))



import numpy as np
import matplotlib.pyplot as plt

mu, sigma = 10, 3 # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)

count, bins, ignored = plt.hist(s, 30, density=True)
plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
               np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
         linewidth=2, color='r')
plt.show()

