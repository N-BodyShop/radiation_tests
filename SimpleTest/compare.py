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

prog_gas = 'gasoline.radiation.pvol'
output_gas = 'simple.00001'
prog_changa = 'ChaNGa'
output_changa = 'simple.000001'

def thetacmp(prog, output) :
    '''run a series of thetas and compare them;
    arguments are the progam to run (e.g. gasoline) and the output basename
    (e.g. simple.00001)'''
    thetas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    colors = ['k', 'b', 'g', 'r', 'c', 'm']
    # Do comparision in reverse order so that 0.0 theta is plotted on top
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
    
        # Note normalization: assumes that all sources have luminousity "1",
        # and therefore a flux of 1 per steradian/r^2
        plt.plot(gr, pf*gr*gr*4*3.14159, colors[i] + '.')
    plt.xlabel('Distance')
    plt.ylabel('Flux times r^2')
        
global df
global pfc, pfg, snapc

def cha_vs_gas(theta) :
    '''Point by point and summary comparison of ChaNGa and gasoline at
    a given theta.'''

    os.system(prog_gas + ' -thetaRadiation ' + str(theta) + ' simple.param > DIAG.gas')
    snapg = pyn.load(output_gas)
    pfg = snapg.g['radFlux']

    os.system(prog_changa + ' -thetaRadiation ' + str(theta) + ' simple.param > DIAG.cha')
    snapc = pyn.load(output_changa)
    pfc = snapc.g['radFlux']
    df = pfc - pfg
    print('rms offset: ',  np.sqrt(sum(df*df)/len(df)))
    print('mean offset: ',  sum(df)/len(df))
    print('max offset: ',  np.sqrt(max(df*df)), ' at ', np.argmax(df*df))
    
    
def read_direct() :
    pfd = np.loadtxt('direct.radFlux', skiprows=1)
    return pfd

def cha_vs_direct(theta) :
    '''Point by point and summary comparison of ChaNGa with direct summation
    a given theta.'''

    os.system(prog_changa + ' -thetaRadiation ' + str(theta) + ' simple.param > DIAG.cha')
    snapc = pyn.load(output_changa)
    pfc = snapc.g['radFlux']
    pfd = read_direct()
    df = (pfc - pfd)/pfd
    print('rms offset: ',  np.sqrt(sum(df*df)/len(df)))
    print('mean offset: ',  sum(df)/len(df))
    print('max offset: ',  np.sqrt(max(df*df)), ' at ', np.argmax(df*df))

def gas_vs_direct(theta) :
    '''Point by point and summary comparison of gasoline with direct summation
    a given theta.'''

    os.system(prog_gas + ' -thetaRadiation ' + str(theta) + ' simple.param > DIAG.cha')
    snapg = pyn.load(output_gas)
    pfg = snapg.g['radFlux']
    pfd = read_direct()
    df = (pfg - pfd)/pfd
    print('rms offset: ',  np.sqrt(sum(df*df)/len(df)))
    print('mean offset: ',  sum(df)/len(df))
    print('max offset: ',  np.sqrt(max(df*df)), ' at ', np.argmax(df*df))
