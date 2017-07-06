#!/usr/bin/python
# -*- coding: utf-8 -*-

# system functions
import os
import sys
# efficient matrix and array operations
import numpy as np
# for reading the data
import csv
import cStringIO
# this is the python package for plotting scientific things
import matplotlib.pyplot as pyplot


def PlotConcentrations(filename, predArr, preyArr, xbounds, ybounds, time):
  # make a page with two image plots for predator and prey concentrations respectively
  fig, axes = pyplot.subplots(1, 2, figsize = (8, 4), dpi = 100) # size is in 'inches'
  # axes refers to figures, in matplotlib slang
  axes[0].imshow(predArr, extent = (xbounds[0], xbounds[1], ybounds[0], ybounds[1]))
  axes[1].imshow(preyArr, extent = (xbounds[0], xbounds[1], ybounds[0], ybounds[1]))
  for ax in axes:
    ax.set(xlabel = 'x', ylabel = 'y')
  axes[0].set(title = 'predators')
  axes[1].set(title = 'prey')
  fig.text(0.01, 0.99, 'time = %0.1f' % time, verticalalignment = 'top')
  plotfilename = os.path.splitext(filename)[0]+'.png'
  fig.savefig(plotfilename, dpi = 100, facecolor = 'white')
  pyplot.close(fig) # release resources


def PlotAverageConc(filenames, globalData):
  # plot average concentrations in one figure
  time, pred, prey = globalData # unpacking first dimension      
  fig, ax = pyplot.subplots(1,1, figsize = (8, 4))
  ax.plot(time, pred, color = 'r', label = 'predators')
  ax.plot(time, prey, color = 'b', label = 'prey')
  ax.legend()
  ax.set(xlabel = 'time', ylabel = 'conc.', title = 'predators & prey')
  plotfilename = os.path.commonprefix(files)+'-avg.pdf'  
  fig.savefig(plotfilename)
  pyplot.close(fig) # release resources


files = sys.argv[1:] # [1:] is a slice operation which gets the list elements from index 1 to the end. Index 0 is skipped because it contains the name of the script itself.


globalData = []  # a list
# lets process all the files
for fn in files:
  print fn,' ...',
  f = open(fn, 'r') # open file in read mode

  # extract data from the header (with #'es in front).
  # We know there are four such lines and we
  # want the numbers after the '='.
  time = float(f.readline().split('=')[1].strip())
  pred = float(f.readline().split('=')[1].strip())
  prey = float(f.readline().split('=')[1].strip())
  grid_size = int(f.readline().split('=')[1].strip())
  print 'time = %s, pred = %f, prey = %f' % (time, pred, prey)

  # Keep this for later. This line adds a 3-tuple to the globalData list.
  globalData.append((time, pred, prey))
  

  # setup a memory file stream for the gridded data which makes up the rest of the file
  io = cStringIO.StringIO(f.read())
  # and give it to this csv reader. It translates tabular text files into a lists of individual table entries
  reader = csv.reader(io, delimiter=' ', quotechar='\n')  
  
  processedLines = []
  for row in reader:
    if not row: continue # filter empty lines
    floats = map(float, row) # typename also acts as constructor function.
    processedLines.append(floats)
  data = np.asarray(processedLines) # turn the list into a numpy array for its multidimensional slicing operations ..

  # get x-y coordinate bounds
  xbounds = np.amin(data[:,0]), np.amax(data[:,0]) 
  ybounds = np.amin(data[:,1]), np.amax(data[:,1])
  # get the concentrations
  predArr = data[:,2].reshape(grid_size, grid_size)
  preyArr = data[:,3].reshape(grid_size, grid_size)
  # plot it
  PlotConcentrations(fn, predArr, preyArr, xbounds, ybounds, time)
  
globalData = zip(*globalData) # transpose operation on builtin lists
PlotAverageConc(files, globalData)
