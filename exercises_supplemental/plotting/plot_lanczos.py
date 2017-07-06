#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:58:16 2017

@author: thierry
"""

# system functions
import os
import sys
# efficient matrix and array operations
import numpy as np
# for reading the data
import csv
# handling command line arguments
import argparse
# this is the python package for plotting scientific things
import matplotlib.pyplot as pyplot

''' data in the file should be tsv. tabular seperated values '''
def read_data(f):
  with f as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter="\t")
    for line in tsv_reader:
      pass
    #now line is the last line
  return int(line[0]),float(line[1]),float(line[2]),float(line[3])

''' data in the file should be tsv. tabular seperated values '''
def read_data_3(fn):
  x_array = []
  y_array = []
  z_array = []
  with open(fn,'r') as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter="\t")
    for line in tsv_reader:
      x = float(line[0])
      y = float(line[1])
      z = float(line[2])
      x_array.append(x)
      y_array.append(y)
      z_array.append(z)
  return x_array,y_array,z_array
  
''' x should be 2xn array'''
def print_correlation(x):
  print(x)
  pyplot.plot(x[0],x[1])
  pyplot.show()

def print_kinetic_energy(x):
  print(x)
  pyplot.scatter(x[0],-1*np.asarray(x[1]))
  pyplot.ylabel(r'$-E_K / V_0 $')
  pyplot.xlabel(r'$\rho$')
  pyplot.show()
  
def print_phase_diagram(x):
  print(x)
  pyplot.scatter(x[0],x[1])
  pyplot.ylabel(r'$\mu$ / $V_0$')
  pyplot.xlabel(r'$t$ / $V_0$')
  pyplot.show()
  
def print_mu(x):
  print(x)
  # data comes
  # mu, \rho, E
  fig1 = pyplot.figure()
  pyplot.plot(x[1],x[0])
  pyplot.ylabel(r'$\mu$')
  pyplot.xlabel(r'$\rho$')
  fig2 = pyplot.figure()
  pyplot.plot(x[1],x[2])
  pyplot.xlabel(r'$\rho$')
  pyplot.ylabel(r'$E_{N_{b}}$')
  pyplot.show()
  
if __name__ == '__main__':
  ''' this takes care of handling command line arguments.
  here: it is use to tell the program thie list of files
  it should convert to images.
  Try the help flag in the command line to see the options.
  usage: python2 visualize.py -h
  '''
  parser = argparse.ArgumentParser(description='Plot data of bose hubbard model.')  
  parser.add_argument('dataFiles', nargs='+', type=argparse.FileType('r'), help='Valid data files must be provided.')
   
  goodArguments, otherArguments = parser.parse_known_args()
  
  for afile in goodArguments.dataFiles:
    print("Found file: %s" % afile)
    W,E0,E1,Ediff = read_data(afile)
#    if( afile.name == "correlation.dat"):
#      x = read_data(afile.name)
#      print_correlation(x)
#    if( afile.name == "kineticEnergy.dat"):
#      x = read_data(afile.name)
#      print_kinetic_energy(x);
#    if( afile.name == "phasediagram.dat"):
#      x = read_data(afile.name)
#      print_phase_diagram(x);
#    if( afile.name == "mu.dat"):
#      x = read_data_3(afile.name)
#      print_mu(x);
