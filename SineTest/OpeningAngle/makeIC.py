import numpy as np
import pynbody as pyn
from pytipsy import wtipsy
import matplotlib.pyplot as plt
import sys, os

glass = pyn.load("glass16.std")
ave_smooth = np.mean(glass.g['smooth'])/2.0


A = 275.

K = [
    [-3.9183979312,    1.7277434311,    -4.4760946235],
    [-3.6818210124,   -4.6196879507,     4.8650067951],
    [-4.8318010154,    3.7694701799,     0.5674512078],
    [-2.2982788483,    1.5017569101,     4.7169462550],
    [-0.2899744545,   -3.0979576042,     1.2700275201],
    [ 1.2629427377,   -1.6617264291,    -2.6004131525],
    [ 1.5882242887,    4.0722593957,     0.6164443955],
    [-2.2533940969,   -2.8064777632,     2.7491546521],
    [-1.4325694148,    3.3247096735,     4.8429906219],
    [ 1.2877420072,   -4.5755170864,    -4.0017233831],
    [ 4.7697038055,    0.5400955240,    -4.2038391810],
    [-3.0132004027,   -1.8712514970,    -2.5144156413],
    [-4.5886202920,    4.3842240224,     1.2468493771],
    [-0.3728171127,    0.1952433788,     4.0740556758],
    [-1.8422321041,    0.9015977953,    -4.4536128228],
    [ 1.9869365065,   -1.0376502859,     1.9588883209],
    [-1.7484848428,   -1.3860287110,     3.7558326553],
    [ 4.8524061986,   -3.2725060415,     0.8265039621],
    [ 3.6632929739,   -4.5975984328,    -0.8901345438],
    [-1.7209027499,    2.7260105431,     3.1924274482],
    [ 4.9733320159,    4.7771822996,    -2.5157917802],
    [ 0.0572383510,   -2.9724262614,    -1.8285497490],
    [ 0.9382397149,   -0.4870226783,    -2.7550973533],
    [ 1.9433613886,    0.3881782794,    -3.7839527412]
    ]

P =[0.8297757341, 3.8911574503, 3.6687303166,
    1.5283482786, 4.1130005370, 4.4817986610,
    2.9719652876, 0.4422411225, 2.8719892948,
    1.7278095211, 5.8721174853, 1.5740076593,
    1.9857146466, 6.2487393885, 6.2733361932,
    2.1777831978, 0.5326044174, 5.5254701709,
    4.5288698416, 3.8756095744, 0.4067369319,
    4.1252582497, 1.3352986404, 4.7749378604]

dims = np.array([4])
sizes = dims**3*len(glass)
names = []
for i in sizes:
    cube = int(np.cbrt(i))
    if cube < 100:
        names.append('0' +str(cube)+'^3')
    else:
        names.append(str(cube)+'^3')

for d, s, n in zip(dims, sizes, names):
    print(d, s, n)
    N = d
    
    s = 0
    f = len(glass.g['x'])
    inc = len(glass.g['x'])
    pos = np.arange(-N/2,N/2)
    
    gx = np.empty(N**3*inc)
    gy = np.empty(N**3*inc)
    gz = np.empty(N**3*inc)
    sx = np.empty(N**3*inc)
    sy = np.empty(N**3*inc)
    sz = np.empty(N**3*inc)
    
    for i in pos:
        for j in pos:
            for k in pos:
                gx[s:f] = glass.g['x'][0:inc]+i
                gy[s:f] = glass.g['y'][0:inc]+j
                gz[s:f] = glass.g['z'][0:inc]+k
                sx[s:f] = glass.g['z'][0:inc]+i
                sy[s:f] = glass.g['y'][0:inc]+j
                sz[s:f] = glass.g['x'][0:inc]+k
                s+=inc
                f+=inc

    gx += .5
    gy += .5
    gz += .5
    sx += .5
    sy += .5
    sz += .5
    
    gx /= N 
    gy /= N 
    gz /= N 
    sx /= N 
    sy /= N 
    sz /= N 

    gdx = np.zeros(N**3*inc)
    gdy = np.zeros(N**3*inc)
    gdz = np.zeros(N**3*inc)
    sdx = np.zeros(N**3*inc)
    sdy = np.zeros(N**3*inc)
    sdz = np.zeros(N**3*inc)
                  
    for k, p in zip(K,P):
        gdx += 1/A*np.sin(2*np.pi*(k[0]*gx + k[1]*gy + k[2]*gz + p))
        gdy += 1/A*np.sin(2*np.pi*(k[0]*gx + k[1]*gy + k[2]*gz + p))
        gdz += 1/A*np.sin(2*np.pi*(k[0]*gx + k[1]*gy + k[2]*gz + p))
        sdx += 1/A*np.sin(2*np.pi*(k[0]*sx + k[1]*sy + k[2]*sz + p))
        sdy += 1/A*np.sin(2*np.pi*(k[0]*sx + k[1]*sy + k[2]*sz + p))
        sdz += 1/A*np.sin(2*np.pi*(k[0]*sx + k[1]*sy + k[2]*sz + p))
    
    gx += gdx
    gy += gdy
    gz += gdz
    sx += sdx
    sy += sdy
    sz += sdz
    
    gsl = np.where(np.fabs(gz) < ave_smooth)
    ssl = np.where(np.fabs(sz) < ave_smooth)
    gx_sl = gx[gsl]
    gy_sl = gy[gsl]
    sx_sl = sx[ssl]
    sy_sl = sy[ssl]

    #plt.figure(figsize = (12,12))
    #plt.title(r'$'+n+'$')
    #plt.xlabel(r'$x$')
    #plt.ylabel(r'$z$')
    #plt.scatter(gx_sl, gy_sl,s=1, c = 'k')
    #plt.scatter(sx_sl, sy_sl,s=1, c = 'orange')
    #plt.show()
    #plt.close()

    ng = len(gx)
    gm = np.ones(ng)/ng
    ns = len(sx)
    nd = 0
    nt = ng + ns + nd

    header = {'time':0,
              'n':nt,
              'ndim':3,
              'ngas':ng,
              'ndark':nd,
              'nstar':ns,}

    gas_particles = {'mass':gm,
                     'x':gx, 'y':gy,'z':gz,
                     'vx':np.zeros(ng),'vy':np.zeros(ng),'vz':np.zeros(ng),
                     'dens':np.zeros(ng),
                     'tempg':np.ones(ng),
                     'h':np.zeros(ng),
                     'phi':np.zeros(ng),
                     'eps':np.ones(1),
                     'zmetal':np.zeros(ng)}
    star_particles = {'mass':np.ones(ns),
                     'x':sx, 'y':sy,'z':sz,
                     'vx':np.zeros(ns),'vy':np.zeros(ns),'vz':np.zeros(ns),
                     'tform':np.ones(ns),
                     'phi':np.zeros(ns),
                     'eps':np.ones(ns),
                     'metals':np.zeros(ns)}
    icName = 'sine_'+n+'.std'
    paramName = 'sine_'+n+".param"
    dirName = "./" + n + "/"
    wtipsy(icName,header,gas_particles,[],star_particles)
