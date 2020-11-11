import numpy as np
import tecplot as tp
from raytrace import raytrace
import sys

# sys.argv[1] is the full path to a directory which contains mutiple 3d*plt file
# import glob
# filenames = sorted(glob.glob(sys.argv[1] + '/3d*plt'))

# for ifile, filename in enumerate(filenames):
#     Area, Flux = raytrace(filename)
#     print('File {:2d}: Northern - Southern open flux is {:.2e} Wb'.format(ifile, Flux['North'] - Flux['South']))

# -----------------------------------------------------------------------------------
# sys.argv[1] is the full path to a 3d*plt file
filename = sys.argv[1]
Area, Flux = raytrace(filename)
# Area, Flux = raytrace('/Volumes/DATA_disk/Frontera_run/SW_test/RUN_46/3d__var_7_t00000118_n00692607.plt')
print('Northern polar cap area is {:.2e} m*m'.format(Area['North']))
print('Southern polar cap area is {:.2e} m*m'.format(Area['South']))
print('Northern open flux is {:.2e} Wb'.format(Flux['North']))
print('Southern open flux is {:.2e} Wb'.format(Flux['South']))
