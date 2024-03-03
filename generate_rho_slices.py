import h5py
import sys
import numpy as np
from scipy.constants import epsilon_0
from ECIX_tools import filemanager as exfm
from rich.progress import track

fname = sys.argv[1]
Nd = 10 #N_discard

grid = exfm.h5_to_dict(fname, group='grid')
AA = h5py.File(fname, 'a')
del AA['grid']
AA.close()

dh = grid['xg'][1] - grid['xg'][0]
grid['xg'] = grid['xg'][Nd:-Nd]
grid['yg'] = grid['yg'][Nd:-Nd]
N = len(grid['zg'])

for ii in track(range(N)):
    phi = exfm.h5_to_dict(fname, group=f'slices/slice{ii}')['phi']
    AA = h5py.File(fname, 'a')
    del AA[f'slices/slice{ii}']
    AA.close()
    rho = np.zeros_like(phi)
    rho[1:-1, 1:-1] = -epsilon_0*(phi[:-2,1:-1] + phi[2:,1:-1] + phi[1:-1,2:] + phi[1:-1,:-2] - 4 * phi[1:-1, 1:-1])/dh**2
    exfm.dict_to_h5({'phi' : phi[Nd:-Nd,Nd:-Nd], 
                     'rho' : rho[Nd:-Nd,Nd:-Nd]},
                    fname, group=f'slices/slice{ii}', readwrite_opts='a')