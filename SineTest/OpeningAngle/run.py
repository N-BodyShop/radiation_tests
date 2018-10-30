import sys, os
import numpy as np

dThetaRad = np.linspace(0,1,11)

for i in dThetaRad:
	out = open('sine_test.param', 'w')
	newLine = 'dThetaRad           = ' + str(i) + '\n'
	with open('default.param', 'r') as f:
		for line in f:
			if line[:9] == 'dThetaRad':
				out.write(newLine)
				print(newLine)
			else:
				out.write(line)
	out.close()

	os.system('/home/grondjj/Code/gasoline/gasoline.radiation.pvol sine_test.param > output.txt')
	os.system('mkdir ./'+str(i)+'/')
	os.system('mv ./sine_064^3.00001* ./output.txt  ./'+str(i)+'/')
	print(str(i) + ' done\n')
