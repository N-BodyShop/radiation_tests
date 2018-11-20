import pynbody as pyn
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os


path = 'simple.00001'
def plotflux(path): 
    # plot flux as a function of distance from the stellar center of mass
    snap = pyn.load('simple.00001')
    pf = snap.g['radFlux']

    # center of mass of the stars
    star_cm = np.sum(snap.s['mass']*snap.s['pos'].transpose(), axis=1)/snap.s['mass'].sum()
    star_cm.units = snap['pos'].units

    # distance of gas from the star cm
    gdx = snap.g['pos'] - star_cm
    gr = np.sqrt(sum((gdx*gdx).transpose()))
    
    plt.plot(gr, pf*gr*gr*4*3.14159, '.')
    plt.xlabel('Distance')
    plt.ylabel('Flux times r^2')

prog = 'gasoline.radiation.pvol'

def thetacmp(prog, output) :
    '''run a series of thetas and compare them'''
    thetas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    colors = ['k', 'b', 'g', 'r', 'c', 'm']
    for i in range(5,-1, -1) :
        os.system(prog + ' -thetaRadiation ' + str(thetas[i]) + ' simple.param')

        snap = pyn.load(output)
        pf = snap.g['radFlux']

        # center of mass of the stars
        star_cm = np.sum(snap.s['mass']*snap.s['pos'].transpose(), axis=1)/snap.s['mass'].sum()
        star_cm.units = snap['pos'].units

        # distance of gas from the star cm
        gdx = snap.g['pos'] - star_cm
        gr = np.sqrt(sum((gdx*gdx).transpose()))
    
        plt.plot(gr, pf*gr*gr*4*3.14159, colors[i] + '.')
    plt.xlabel('Distance')
    plt.ylabel('Flux times r^2')
        
