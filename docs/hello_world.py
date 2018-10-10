import matplotlib.pyplot as plt
import numpy as np

def plot_histo():
	aArrayOfNumbers = np.random.randn(1000)
	fig1 = plt.figure()
    	ax1 = fig1.add_subplot(111)
    	ax1.hist(aArrayOfNumbers)
	plt.show()
    	ax1.set_xlabel(r'distance from nearest vessel / $\mu m')
    	ax1.set_ylabel(r'probability')
if __name__ == '__main__':
	print('hello world')
	plot_histo()


