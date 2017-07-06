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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pyplot

interactiveView=False
''' data in the file should be tsv. tabular seperated values '''
def read_data(fn):
  if not fn == 'mu.dat':
    x_array = []
    y_array = []
    with open(fn,'r') as tsv_file:
      tsv_reader = csv.reader(tsv_file, delimiter="\t")
      for line in tsv_reader:
        x = float(line[0])
        y = float(line[1])
        x_array.append(x)
        y_array.append(y)
    return x_array,y_array
  else:
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

''' data in the file should be tsv. tabular seperated values '''
#def read_data_3(fn):
#  x_array = []
#  y_array = []
#  z_array = []
#  with open(fn,'r') as tsv_file:
#    tsv_reader = csv.reader(tsv_file, delimiter="\t")
#    for line in tsv_reader:
#      x = float(line[0])
#      y = float(line[1])
#      z = float(line[2])
#      x_array.append(x)
#      y_array.append(y)
#      z_array.append(z)
#  return x_array,y_array,z_array

def print_data(pdf):
  x = read_data(afile)
  if afile == 'correlation.dat':
    print_correlation(x,pdf)
  elif afile == 'kineticEnergy.dat':
    print_kinetic_energy(x,pdf)
  elif afile == 'mu.dat':
    print_mu(x,pdf)
  elif afile == 'phasediagram.dat':
    print_phase_diagram(x,pdf)
    
''' x should be 2xn array'''
def print_correlation(x,pdf):
  fig_correl = pyplot.figure()
  pyplot.plot(x[0],x[1])
  if interactiveView:
    pyplot.show()
  else:
    pdf.savefig(fig_correl)

def print_kinetic_energy(x,pdf):
  fig_kinetic = pyplot.figure()
  pyplot.scatter(x[0],-1*np.asarray(x[1]))
  pyplot.ylabel(r'$-E_K / V_0 $')
  pyplot.xlabel(r'$\rho$')
  pyplot.grid()
  if interactiveView:
    pyplot.show()
  else:
    pdf.savefig(fig_kinetic)
  
def print_mu(x,pdf):
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
  if interactiveView:
    pyplot.show()
  else:
    pdf.savefig(fig1)
    pdf.savefig(fig2)
    
def print_phase_diagram(x,pdf):
  fig_phase = pyplot.figure()
  #pyplot.scatter(x[0][0::2],x[1][0::2])
  pyplot.scatter(x[0],x[1])
  pyplot.ylabel(r'$\mu$ / $V_0$')
  pyplot.xlabel(r'$t$ / $V_0$')
  if interactiveView:
    pyplot.show()
  else:
    pdf.savefig(fig_phase)
    
if __name__ == '__main__':
  ''' this takes care of handling command line arguments.
  here: it is use to tell the program thie list of files
  it should convert to images.
  Try the help flag in the command line to see the options.
  usage: python2 visualize.py -h
  '''
  
  #parser = argparse.ArgumentParser(description='Plot data of bose hubbard model.')  
  #parser.add_argument('dataFiles', nargs='+', type=argparse.FileType('r'), help='Valid data files must be provided.')
  #goodArguments, otherArguments = parser.parse_known_args()
  
  outputFilesOfBosehubbard = ["correlation.dat","kineticEnergy.dat","mu.dat","phasediagram.dat",]
  with PdfPages('bosehubbard_result.pdf') as pdf:
    for afile in outputFilesOfBosehubbard:
      if os.path.isfile(afile):
        print("found file: %s" % afile)
        print_data(pdf)
      
      else:
        print("file %s not found" % afile)
    
#    if( afile.name == "correlation.dat"):
#      x = read_data(afile.name)
#      print_correlation(x)
#    if( afile.name == "kineticEnergy.dat"):
#      x = read_data(afile.name)
#      print_kinetic_energy(x);
#    if( afile.name == "phasendigram.dat"):
#      x = read_data(afile.name)
#      print_phase_diagram(x);
#    if( afile.name == "mu.dat"):
#      x = read_data_3(afile.name)
#      print_mu(x);
