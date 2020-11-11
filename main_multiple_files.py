import numpy as np
import tecplot as tp
from raytrace import raytrace
import sys
import matplotlib.pyplot as plt

# sys.argv[1] is the full path to a directory which contains mutiple 3d*plt file
import glob
filenames = sorted(glob.glob(sys.argv[1] + '/3d*plt'))
N = len(filenames)
NorthFlux = np.zeros(N)
SouthFlux = np.zeros(N)
NorthArea = np.zeros(N)
SouthArea = np.zeros(N)
DiffFlux = np.zeros(N)

for ifile, filename in enumerate(filenames):
    Area, Flux = raytrace(filename)
    # print('File {:2d}: Northern - Southern open flux is {:.2e} Wb'.format(ifile, Flux['North'] - Flux['South']))
    NorthFlux[ifile] = Flux['North']
    SouthFlux[ifile] = Flux['South']
    NorthArea[ifile] = Area['North']
    SouthArea[ifile] = Area['South']
    DiffFlux[ifile]  = NorthFlux[ifile] - SouthFlux[ifile]

np.savez(sys.argv[1] + '/pymat.npz', NorthFlux=NorthFlux, SouthFlux=SouthFlux, NorthArea=NorthArea, SouthArea=SouthArea, DiffFlux=DiffFlux)

f = plt.figure(figsize=(8,4))
plt.xlabel(r'$time (iteration)$', fontsize=16)
plt.ylabel(r'$Flux Diff$', fontsize=16)
plt.plot(DiffFlux, color='blue', lw=2, label='Fluxdifference(north - south)')
plt.legend(fontsize = 12, loc=0, borderaxespad=0.1)
f.tight_layout(); plt.show(block=False)
f.savefig(sys.argv[1] + '/FluxDiff.png')
plt.close(f)

f = plt.figure(figsize=(8,4))
plt.xlabel(r'$time (iteration)$', fontsize=16)
plt.ylabel(r'$Flux$', fontsize=16)
plt.plot(NorthFlux, color='red',  lw=2, label='Northern Flux')
plt.plot(SouthFlux, color='blue', lw=2, label='Southern Flux')
plt.legend(fontsize = 12, loc=0, borderaxespad=0.1)
f.tight_layout(); plt.show(block=False)
f.savefig(sys.argv[1] + '/Flux.png')
plt.close(f)

f = plt.figure(figsize=(8,4))
plt.xlabel(r'$time (iteration)$', fontsize=16)
plt.ylabel(r'$Area$', fontsize=16)
plt.plot(NorthArea, color='red',  lw=2, label='Northern Area')
plt.plot(SouthArea, color='blue', lw=2, label='Southern Area')
plt.legend(fontsize = 12, loc=0, borderaxespad=0.1)
f.tight_layout(); plt.show(block=False)
f.savefig(sys.argv[1] + '/Area.png')
plt.close(f)

# -----------------------------------------------------------------------------------
# sys.argv[1] is the full path to a 3d*plt file
# filename = sys.argv[1]
# Area, Flux = raytrace(filename)
# Area, Flux = raytrace('/Volumes/DATA_disk/Frontera_run/SW_test/RUN_46/3d__var_7_t00000118_n00692607.plt')
# print('Northern polar cap area is {:.2e} m*m'.format(Area['North']))
# print('Southern polar cap area is {:.2e} m*m'.format(Area['South']))
# print('Northern open flux is {:.2e} Wb'.format(Flux['North']))
# print('Southern open flux is {:.2e} Wb'.format(Flux['South']))
