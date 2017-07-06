# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 14:04:00 2017

This file was created as part of the problem sheet 
within the computerphysics lecture at Saarland University.

Feel free to copy and modify for whatever purpuses, but
mainly try to please your lecturer.

@author: thierry

@dependencies: best as scientist to have the complete scipy
package at stack. Usually shiped with every major distribution.
"""
import numpy as np
import scipy
import scipy.integrate
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

#systemparameter, here static
r=1;
s=1./2;
p=5;
q=5;
m=1; #only for Volterra-Lotka relevant

#timevector
tmax=5;
tmin=0;
Dt=0.01;
tspan = np.arange(tmin, tmax, Dt)


#initial condition
#diffetnet initial conditions
firstCoordinate = 1./9.
secondCoordinate = np.arange(0,1,0.1)
X0 = [np.repeat(firstCoordinate,len(secondCoordinate)), secondCoordinate]
X0 = np.asarray(X0)

def equationODE(x,t,IndODE='rb'):
  u,v = x
  if IndODE == 'rb':
    dxdt = [r*u-p*u*v,-s*v+q*u*v];
  elif IndODE == 'vl':
    print("excercise 2 (Hopf-Bifurkation): implement this")
  else:
    print("no know option!")
  #technical detail: conversation to numpy array
  return np.asarray(dxdt)

#solution due to Euler
#use equations explicitly stated in equationODE
def euler(Dt,tspan,x0,IndODE='rb'):
  print("excercise: implement this");
  #return dx;

#solution due to Runge-Kutta scheme
#use equations explicitly stated in equationODE
def rk(Dt,tspan,x0,IndODE='rb'):
  print("excercise: implement this")
#  return dx;

def print_example_phase_space():
  list_of_solutions=[]
  list_of_labels=[]
  for aInitialCondition in X0.transpose():
    sol = scipy.integrate.odeint(equationODE,aInitialCondition,tspan, args=('rb',))
    #solEuler = euler(Dt,tspan,aInitialCondition,IndODE='rb')
    #solRK = rk(Dt,tspan,aInitialCondition,IndODE='rb')
    
    if 0: #some output  
      print('aInitialCondition: %s' % aInitialCondition)
      print('solution         : %s' % sol)
    list_of_solutions.append(sol)
    list_of_labels.append(aInitialCondition)
    
    if 0:
      list_of_solutions.append(solEuler)
      list_of_labels.append("euler %s" % aInitialCondition)
      
    if 0:
      list_of_solutions.append(solRK)
      list_of_labels.append("rk %s" % aInitialCondition)
      
  
  fig1,ax = plt.subplots()
  fig1.suptitle('Phasespace trajectories')
  for sol,label in zip(list_of_solutions, list_of_labels):
    #phase space
    ax.plot(sol[:,0],sol[:,1], label=label)
  
  ax.legend(prop={'size':6})
  ax.grid()
  ax.set_xlabel('u')
  ax.set_ylabel('v')
  
  viewLive=True
  if viewLive:
    plt.show() #this opens an iteractive window
  else:
    #this writes the figure to pdf. note: other formats are also possible
    #see documentation
    with PdfPages('example_solution.pdf') as pdf:
      pdf.savefig(fig1)
  


if __name__ == '__main__':
  print_example_phase_space();

