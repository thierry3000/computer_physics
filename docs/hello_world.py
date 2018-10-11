import matplotlib.pyplot as plt
import numpy as np

def plot_histo():
  aArrayOfNumbers = np.random.randn(1000)
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  ax1.hist(aArrayOfNumbers)
  ax1.set_xlabel(r'numbers drawn from $p_{normal}$') #note the latex
  ax1.set_ylabel(r'counts')
  plt.show()
  fig1.savefig("normal_dist.pdf")

if __name__ == '__main__':
  print('hello world')
  #plot_histo()


