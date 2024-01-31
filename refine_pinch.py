import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os
import time
import shutil
from collections import deque
import sys
#sys.path.append('..')
import ecloud_xsuite_filemanager as exfm
import refinement_helpers as rh
import volume_helpers as vh

from TricubicInterpolation import pyTricubic as pyTI

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename")
parser.add_argument("--MTI", nargs="?", type=float, default=1.)
parser.add_argument("--MLI", nargs="?", type=float, default=1.)
parser.add_argument("--DTO", nargs="?", type=float, default=1.)
parser.add_argument("--DLO", nargs="?", type=float, default=1.)
parser.add_argument("--symmetric2D", action="store_true")
args = parser.parse_args()


start_time = time.time()
N_nodes_discard = 10
magnify_transverse_in = args.MLI
magnify_longitudinal_in = args.MLI
demagnify_transverse_out = args.DTO
demagnify_longitudinal_out = args.DLO
compression_opts = 0
do_symmetric2D = args.symmetric2D
debug = False

symm_str = ''
if do_symmetric2D: symm_str = '_symm2D'


fname = args.filename
pinch_in = args.filename[:-3]

print(pinch_in)
# pinch_out = 'refined_'+pinch_in+symm_str+'_MTI%.1f_MLI%.1f_DTO%.1f_DLO%.1f.h5'%(magnify_transverse_in, magnify_longitudinal_in, demagnify_transverse_out, demagnify_longitudinal_out)
pinch_out = f"refined_{pinch_in}{symm_str}_MTI{args.MTI:.1f}_MLI{args.MLI:.1f}_DTO{args.DTO:.1f}_DLO{args.DLO:.1f}.h5"
out_fname = pinch_out
out_efname = 'e_' + pinch_out

print('Pinch in: ' + pinch_in)
print('Magnify = %f'%magnify_transverse_in)
pic_out, pic_in, zg = rh.setup_pic(fname, magnify=magnify_transverse_in, N_nodes_discard=N_nodes_discard, symmetric_slice_2D=do_symmetric2D)



dz_original = zg[1]-zg[0]
dz_inside = dz_original/magnify_longitudinal_in
zg_inside = np.linspace(zg[0],zg[-1],int((zg[-1]-zg[0])/dz_inside)+1)

dz_new = dz_original*demagnify_longitudinal_out
zg_new = np.linspace(zg[3],zg[-4],int((zg[-4]-zg[3])/dz_new)+1) # skip three first and three last slices
#dz_new = demagnify_longitudinal_out*dz_original


Nd = N_nodes_discard + 3

dx_new = pic_out.dx*demagnify_transverse_out
nx_new = int(np.rint((pic_out.xg[-Nd-1] - pic_out.xg[Nd])/dx_new)) + 1
xg_new = np.linspace( pic_out.xg[Nd], pic_out.xg[-Nd-1], nx_new)

dy_new = pic_out.dy*demagnify_transverse_out
ny_new = int(np.rint((pic_out.yg[-Nd-1] - pic_out.yg[Nd])/dy_new)) + 1
yg_new = np.linspace( pic_out.yg[Nd], pic_out.yg[-Nd-1], ny_new)

#print(zg[:5])
#print(zg_new[:5])
#print(zg[-5:])
#print(zg_new[-5:])
#print(zg.shape, zg_inside.shape)

#for ii,zz in enumerate(zg_inside):
#    rh.get_slice(pic_out, pic_in, fname, zz, symmetric_slice_2D = do_symmetric2D)
#    print('%d/%d'%(ii,len(zg_inside)))

phi_slices = deque()
ex_slices = deque()
ey_slices = deque()
ez_slices = deque()
z_slices = deque()

for i in range(7):
    z_slices.append(zg_inside[i])
    phi_slices.append(rh.get_slice(pic_out,pic_in,fname, zg_inside[i], symmetric_slice_2D=do_symmetric2D))
    ex_slices.append(rh.ex_from_phi_slice(phi_slices[i], pic_in.Dh))
    ey_slices.append(rh.ey_from_phi_slice(phi_slices[i], pic_in.Dh))
ez_slices.append(np.zeros_like(phi_slices[0]))
for i in range(1,6):
    ez_slices.append(rh.ez_from_two_phi_slices(phi_slices[i-1], phi_slices[i+1], dz_inside))

Sx = nx_new//2 if do_symmetric2D else 0
Sy = ny_new//2 if do_symmetric2D else 0

grid_dict = {'xg': xg_new[Sx:],
             'yg': yg_new[Sy:],
             'zg': zg_new,
             'x0': xg_new[Sx],
             'y0': yg_new[Sy],
             'z0': zg_new[0],
             'dx': xg_new[Sx+1] - xg_new[Sx],
             'dy': yg_new[Sy+1] - yg_new[Sy],
             'dz': zg_new[1] - zg_new[0],
             'Nx': len(xg_new),
             'Ny': len(yg_new),
             'Nz': len(zg_new)
            }

settings_dict = {'magnify_transverse_in': magnify_transverse_in,  
                 'magnify_longitudinal_in': magnify_longitudinal_in, 
                 'demagnify_transverse_out': demagnify_transverse_out,
                 'demagnify_longitudinal_out': demagnify_longitudinal_out,
                 'symmetric2D': do_symmetric2D
                }


exfm.dict_to_h5(grid_dict, out_fname, compression_opts=compression_opts, group='grid', readwrite_opts='w')
exfm.dict_to_h5(settings_dict, out_fname, compression_opts=compression_opts, group='settings', readwrite_opts='a')
exfm.dict_to_h5(grid_dict, out_efname, compression_opts=compression_opts, group='grid', readwrite_opts='w')
exfm.dict_to_h5(settings_dict, out_efname, compression_opts=compression_opts, group='settings', readwrite_opts='a')


slice_index = 7
kk = 0
while kk < len(zg_new):
    z_new = zg_new[kk]
    z_up = (z_slices[3]+z_slices[4])/2
    z_lo = (z_slices[2]+z_slices[3])/2
    if z_new >= z_lo and z_new < z_up:
        ## extract transverse slices and check transverse new grid
        phi_out_slice = rh.exact_slice(phi_slices, pic_in.xg, pic_in.yg, z_slices, xg_new, yg_new, z_new)
        ex_out_slice = rh.exact_slice(ex_slices, pic_in.xg, pic_in.yg, z_slices, xg_new, yg_new, z_new)
        ey_out_slice = rh.exact_slice(ey_slices, pic_in.xg, pic_in.yg, z_slices, xg_new, yg_new, z_new)
        ez_out_slice = rh.exact_slice(ez_slices, pic_in.xg, pic_in.yg, z_slices, xg_new, yg_new, z_new)

#        if do_symmetric2D: rh.fix_phi(phi_out_slice, Sx, Sy)

        exfm.dict_to_h5({'phi' : phi_out_slice[Sx:,Sy:,:]}, out_fname, compression_opts=compression_opts, group='slices/slice%d'%kk, readwrite_opts='a')
        exfm.dict_to_h5({'ex' : ex_out_slice[Sx:,Sy:,:]}, out_efname, compression_opts=compression_opts, group='slices/ex_slice%d'%kk, readwrite_opts='a')
        exfm.dict_to_h5({'ey' : ey_out_slice[Sx:,Sy:,:]}, out_efname, compression_opts=compression_opts, group='slices/ey_slice%d'%kk, readwrite_opts='a')
        exfm.dict_to_h5({'ez' : ez_out_slice[Sx:,Sy:,:]}, out_efname, compression_opts=compression_opts, group='slices/ez_slice%d'%kk, readwrite_opts='a')
        kk += 1
    elif z_new >= z_up:
        z_slices.append(zg_inside[slice_index])
        phi_slices.append(rh.get_slice(pic_out,pic_in,fname, z_slices[-1], symmetric_slice_2D=do_symmetric2D))
        ex_slices.append(rh.ex_from_phi_slice(phi_slices[-1], pic_in.Dh))
        ey_slices.append(rh.ey_from_phi_slice(phi_slices[-1], pic_in.Dh))
        ez_slices.append(rh.ez_from_two_phi_slices(phi_slices[-3], phi_slices[-1], dz_inside))

        z_slices.popleft()
        phi_slices.popleft()
        ex_slices.popleft()
        ey_slices.popleft()
        ez_slices.popleft()
        slice_index += 1
    else:
        print('Lost synchronization')
    print('%d/%d, %d/%d'%(kk,len(zg_new),slice_index,len(zg_inside)))

del phi_slices
del ex_slices
del ey_slices
del ez_slices
del pic_out
del pic_in

print('Max x = %f'%max(xg_new))
print('Max y = %f'%max(yg_new))

perc = 0.
print('Reading... ')

ti_method = 'Exact-Mirror2' if do_symmetric2D else 'Exact'

plot_hist_str = 'plt.hist(bins[:-1], bins, weights=counts, histtype=\'step\', log=True)'
hist_counts, bins = np.histogram([], bins=2000, range=(-13,0))
log10_vx_hist = np.zeros_like(hist_counts)
log10_vy_hist = np.zeros_like(hist_counts)
log10_vz_hist = np.zeros_like(hist_counts)
max_log10_vx = -1000
max_log10_vy = -1000
max_log10_vz = -1000
max_index_x = -1
max_index_y = -1
max_index_z = -1


kk=0
grid = exfm.h5_to_dict(out_fname, group='grid')
phi0 = exfm.h5_to_dict(out_fname, group='slices/slice%d'%kk)['phi']
ex0 = exfm.h5_to_dict(out_efname, group='slices/ex_slice%d'%kk)['ex']
ey0 = exfm.h5_to_dict(out_efname, group='slices/ey_slice%d'%kk)['ey']
ez0 = exfm.h5_to_dict(out_efname, group='slices/ez_slice%d'%kk)['ez']
for kk in range(1,len(zg_new)):
    if (1.*kk)/len(zg_new) >= perc:
        print('%d%%... '%(int(kk/slice_index*100)))
        perc += 0.1
    phi1 = exfm.h5_to_dict(out_fname, group='slices/slice%d'%kk)['phi']
    ex1 = exfm.h5_to_dict(out_efname, group='slices/ex_slice%d'%kk)['ex']
    ey1 = exfm.h5_to_dict(out_efname, group='slices/ey_slice%d'%kk)['ey']
    ez1 = exfm.h5_to_dict(out_efname, group='slices/ez_slice%d'%kk)['ez']
    TIphi = pyTI.Tricubic_Interpolation(A=np.array([phi0,phi1]).transpose(1,2,0,3), dx=grid['dx'], dy=grid['dy'], dz=grid['dz'],method=ti_method)
    TIex = pyTI.Tricubic_Interpolation(A=np.array([ex0,ex1]).transpose(1,2,0,3), dx=grid['dx'], dy=grid['dy'], dz=grid['dz'], method=ti_method)
    TIey = pyTI.Tricubic_Interpolation(A=np.array([ey0,ey1]).transpose(1,2,0,3), dx=grid['dx'], dy=grid['dy'], dz=grid['dz'], method=ti_method)
    TIez = pyTI.Tricubic_Interpolation(A=np.array([ez0,ez1]).transpose(1,2,0,3), dx=grid['dx'], dy=grid['dy'], dz=grid['dz'], method=ti_method)
    iz = 0
    vx = np.empty([TIphi.ix_bound_up+1 - TIphi.ix_bound_low, TIphi.iy_bound_up+1 - TIphi.iy_bound_low])
    vy = np.empty([TIphi.ix_bound_up+1 - TIphi.ix_bound_low, TIphi.iy_bound_up+1 - TIphi.iy_bound_low])
    vz = np.empty([TIphi.ix_bound_up+1 - TIphi.ix_bound_low, TIphi.iy_bound_up+1 - TIphi.iy_bound_low])
    for ix in range(TIphi.ix_bound_low, TIphi.ix_bound_up+1):
        for iy in range(TIphi.iy_bound_low, TIphi.iy_bound_up+1):
            vx[ix,iy] = vh.var_x(ix, iy, iz, TIphi, TIex, dx_new)    
            vy[ix,iy] = vh.var_y(ix, iy, iz, TIphi, TIey, dy_new)    
            vz[ix,iy] = vh.var_z(ix, iy, iz, TIphi, TIez, dz_new)    

    lvx = np.log10(vx)
    lvy = np.log10(vy)
    lvz = np.log10(vz)

    log10_vx_hist += np.histogram(lvx.flatten(), bins=2000, range=(-13,0))[0]
    log10_vy_hist += np.histogram(lvy.flatten(), bins=2000, range=(-13,0))[0]
    log10_vz_hist += np.histogram(lvz.flatten(), bins=2000, range=(-13,0))[0]

    if max_log10_vx < np.max(lvx):
        max_index_x = np.unravel_index(np.argmax(lvx, axis=None), lvx.shape) + (kk,)
        max_log10_vx = np.max(lvx)

    if max_log10_vy < np.max(lvy):
        max_index_y = np.unravel_index(np.argmax(lvy, axis=None), lvy.shape) + (kk,)
        max_log10_vy = np.max(lvy)

    if max_log10_vz < np.max(lvz):
        max_index_z = np.unravel_index(np.argmax(lvz, axis=None), lvz.shape) + (kk,)
        max_log10_vz = np.max(lvz)

    

    phi0, ex0, ey0, ez0 = phi1, ex1, ey1, ez1

print()
print('Comparison complete.')

print('Log10 Vx max: ', max_log10_vx, ' at ', max_index_x)
print('Log10 Vy max: ', max_log10_vy, ' at ', max_index_y)
print('Log10 Vz max: ', max_log10_vz, ' at ', max_index_z)

stats_dict = {'max_log10_vx': max_log10_vx,
              'max_log10_vy': max_log10_vy,
              'max_log10_vz': max_log10_vz,
              'max_index_x': max_index_x,
              'max_index_y': max_index_y,
              'max_index_z': max_index_z,
              'bins': bins,
              'log10_vx_hist': log10_vx_hist,
              'log10_vy_hist': log10_vy_hist,
              'log10_vz_hist': log10_vz_hist,
              'plot_hist_str': plot_hist_str,
              'ti_method': ti_method
             }

exfm.dict_to_h5(stats_dict, out_fname, compression_opts=compression_opts, group='stats', readwrite_opts='a')
exfm.dict_to_h5(stats_dict, out_efname, compression_opts=compression_opts, group='stats', readwrite_opts='a')

end_time = time.time()
print('Running time: %f mins'%((end_time-start_time)/60.))

# shutil.copyfile(out_fname, pinches_folder+out_fname)
#shutil.copyfile(out_efname, pinches_folder+out_efname)

