import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
from glob import glob

path = "./"
folders = glob(path+"[0-9].[0-9]")
N = 64**3
rays = []
cells = []
theta = []
rmse = []

t0f = np.array(pyn.load(path+"0.0/sine_064^3.00001").g['radFlux'])

for folder in folders:
	theta.append(float(folder.split('/')[-1]))
	tf = np.array(pyn.load(folder + "/sine_064^3.00001").g['radFlux'])
	rmse.append(np.sqrt(np.sum(((tf - t0f)/t0f)**2)/len(tf)))
	switch = 0
	with open(folder+"/output.txt", 'r') as f:
		for line in f:
			if line[:7] == "nActive":
				if switch == 0:
					cells.append(int(line.split(',')[-2].split(' ')[-1]))
					rays.append(int(line.split(',')[-1].split(' ')[-1]))
					switch += 1

theta = np.array(theta)
rmse = np.array(rmse)
rays = np.array(rays)
cells = np.array(cells)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot([-0.1,1.1],[N,N],c='k',linestyle='-.', label = r"$N_{\rm source}$", zorder = -1)
ax2.scatter(theta, rmse, c = "C3", marker = "o")
ax1.set_xlabel(r"$\theta_{\rm open}$")
ax1.set_xlim(-0.1,1.1)
ax2.set_ylabel(r"${\rm RMS}\left[\left(F_{\theta_{\rm open}} - F_{0}\right)/F_{0}\right]$", color = 'C3', labelpad = 10)
ax2.tick_params('y', colors = 'C3')
ax2.set_yscale('log')
ax2.set_ylim(1e-4,1e-1)
ax1.scatter(theta, rays/N, c = "C0", marker = "s")
ax1.set_ylabel(r"$N_{\rm rays}/N$", color = 'C0')
ax1.tick_params('y', colors = 'C0')
ax1.set_ylim(1e1,1e6)
ax1.set_yscale('log')
ax1.legend(loc = 'lower right')
ax1.tick_params(axis='y', which='both', colors='C0')
ax2.tick_params(axis='y', which='both', colors='C3')
plt.tight_layout()
#plt.savefig("opening_angle.png")
#plt.savefig("opening_angle.pdf")
plt.show()
plt.close()
