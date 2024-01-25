import numpy as np

import sys
#sys.path.append('..')
#sys.path.append('../PyFRIENDS')
import ecloud_xsuite_filemanager as exfm
from TricubicInterpolation import cTricubic as ti
import PyPIC.PyPIC_Scatter_Gather as PyPICSC
import PyPIC.geom_impact_poly as PyPICpoly
import PyPIC.MultiGrid as PyPICMG

def symmetrize(A):
    shape = A.shape
    assert len(shape) == 2
    assert shape[0] % 2 == 1
    assert shape[1] % 2 == 1
    ix0 = shape[0] // 2
    iy0 = shape[1] // 2

    newA = np.empty([ix0+1,iy0+1])
    newA[0,0] = A[ix0,iy0]
    newA[1:,0] = 0.5*(A[0:ix0,iy0][::-1] + A[ix0+1:,iy0])
    newA[0,1:] = 0.5*(A[ix0,0:iy0][::-1] + A[ix0,iy0+1:])
    newA[1:,1:] = 0.25*( A[ix0+1:,iy0+1:] +
                         A[0:ix0,iy0+1:][::-1,:] + 
                         A[ix0+1:,0:iy0][:,::-1] + 
                         A[0:ix0,0:iy0][::-1,::-1]
                       )

    retA=np.empty_like(A)

    retA[ix0,iy0] = newA[0,0]

    retA[ix0+1:,iy0] = newA[1:,0]
    retA[0:ix0,iy0] = newA[1:,0][::-1]
    retA[ix0,iy0+1:] = newA[0,1:]
    retA[ix0,0:iy0] = newA[0,1:][::-1]

    retA[ix0+1:,iy0+1:] = newA[1:,1:]
    retA[0:ix0,iy0+1:] = newA[1:,1:][::-1,:]
    retA[ix0+1:,0:iy0] = newA[1:,1:][:,::-1]
    retA[0:ix0,0:iy0] = newA[1:,1:][::-1,::-1]

    return retA

def setup_pic(fname, magnify=2., N_nodes_discard=10, symmetric_slice_2D=True):

    grid = exfm.h5_to_obj(fname, group='grid')
    ob0 = exfm.h5_to_obj(fname, group='slices/slice0')
    Dh_magnify = (grid.xg[1]-grid.xg[0])/magnify
    x_magnify_min = grid.xg[N_nodes_discard]
    x_magnify_max = grid.xg[-N_nodes_discard-1]
    y_magnify_min = grid.yg[N_nodes_discard]
    y_magnify_max = grid.yg[-N_nodes_discard-1]
    # x_magnify = -grid.xg[N_nodes_discard]
    # y_magnify = -grid.yg[N_nodes_discard]
    
    if symmetric_slice_2D:
        pic_rho = symmetrize(ob0.rho)
        pic_phi = symmetrize(ob0.phi)
    else:
        pic_rho = ob0.rho.copy()
        pic_phi = ob0.phi.copy()
    xg_out = grid.xg.copy()
    yg_out = grid.yg.copy()
    zg_out = grid.zg.copy()
    del grid
    del ob0

    chamb = PyPICpoly.polyg_cham_geom_object({'Vx':np.array([xg_out[-1], xg_out[0], xg_out[0], xg_out[-1]]),
                                       'Vy':np.array([yg_out[-1], yg_out[-1], yg_out[0], yg_out[0]]),
                                       'x_sem_ellip_insc':1e-3,
                                       'y_sem_ellip_insc':1e-3})

    pic = PyPICSC.PyPIC_Scatter_Gather(xg = xg_out, yg = yg_out)
    pic.phi = pic_phi
    pic.rho = pic_rho
    pic.chamb = chamb

    # Internal pic
    picdg = PyPICMG.AddInternalGrid(pic, 
            x_min_internal=x_magnify_min,
            x_max_internal=x_magnify_max, 
            y_min_internal=y_magnify_min, 
            y_max_internal=y_magnify_max, 
            Dh_internal=Dh_magnify, 
            N_nodes_discard = N_nodes_discard)
    # picdg = PyPICMG.AddInternalGrid(pic, 
    #         x_min_internal=-x_magnify,
    #         x_max_internal=x_magnify, 
    #         y_min_internal=-y_magnify, 
    #         y_max_internal=y_magnify, 
    #         Dh_internal=Dh_magnify, 
    #         N_nodes_discard = N_nodes_discard)
    picinside = picdg.pic_internal 

    picinside.rho = np.reshape(pic.gather_rho(picinside.xn, picinside.yn),
            (picinside.Nxg, picinside.Nyg))
    picinside.solve(flag_verbose = True, pic_external=pic)

    return pic, picinside, zg_out

def get_slice(picoutside, picinside, fname, zslice, symmetric_slice_2D=True):

    grid = exfm.h5_to_obj(fname, group='grid')

    dz = grid.zg[1] - grid.zg[0]
    iz_f = (zslice-grid.zg[0])/dz

    i_left = int(np.floor(iz_f))
    i_right = i_left + 1

    if zslice == grid.zg[-1]:
        ilast = len(grid.zg)-1
        oblast = exfm.h5_to_obj(fname, group=f'slices/slice{ilast}')
        if symmetric_slice_2D:
            phi = symmetrize(oblast.phi)
            rho = symmetrize(oblast.rho)
        else:
            phi = oblast.phi.copy()
            rho = oblast.rho.copy()
        del oblast
    else:
        obleft = exfm.h5_to_obj(fname, group=f'slices/slice{i_left}')
        obright = exfm.h5_to_obj(fname, group=f'slices/slice{i_right}')
        if symmetric_slice_2D:
            rho0 = symmetrize(obleft.rho)
            phi0 = symmetrize(obleft.phi)
            rho1 = symmetrize(obright.rho)
            phi1 = symmetrize(obright.phi)
        else:
            rho0 = obleft.rho.copy()
            phi0 = obleft.phi.copy()
            rho1 = obright.rho.copy()
            phi1 = obright.phi.copy()
        z0 = grid.zg[i_left]
        phi = phi0 + (phi1 - phi0)*((zslice-z0)/dz)
        rho = rho0 + (rho1 - rho0)*((zslice-z0)/dz)
        del phi1, phi0, rho1, rho0
        del obleft, obright

    picoutside.phi = phi
    picoutside.rho = rho


    picinside.rho = np.reshape( picoutside.gather_rho(picinside.xn, picinside.yn),
                               (picinside.Nxg, picinside.Nyg))
    picinside.solve(flag_verbose = True, pic_external = picoutside)

    phi_refined = picinside.gather_phi(picinside.xn, picinside.yn)
    phi_refined = phi_refined.reshape(picinside.Nxg, picinside.Nyg)

    return phi_refined


def ex_from_phi_slice(phi_slice, Dh):
    ex_slice = np.zeros_like(phi_slice)
    ex_slice[1:-1,:] = (0.5*(phi_slice[2:,:] - phi_slice[0:-2,:]))/Dh
    return ex_slice

def ey_from_phi_slice(phi_slice, Dh):
    ey_slice = np.zeros_like(phi_slice)
    ey_slice[:,1:-1] = (0.5*(phi_slice[:,2:] - phi_slice[:,0:-2]))/Dh
    return ey_slice

def ez_from_two_phi_slices(phi_slice0, phi_slice2, dz_in):
    ez_slice = (0.5*(phi_slice2 - phi_slice0))/dz_in
    return ez_slice

#def phi_n_e_slices(pic_out, pic_in, fname, islice, dz, max_slice, symmetric_slice_2D=True):
#    slices = np.zeros([3,pic_in.Nxg, pic_in.Nyg])
#    slices[0,:,:] = get_slice(pic_out, pic_in, fname, islice - 1, symmetric_slice_2D) 
#    slices[1,:,:] = get_slice(pic_out, pic_in, fname, islice, symmetric_slice_2D) 
#    if islice + 1 < max_slice:
#        slices[2,:,:] = get_slice(pic_out, pic_in, fname, islice + 1, symmetric_slice_2D) 
#
#    ex_slice = np.zeros([pic_in.Nxg,pic_in.Nyg])
#    ey_slice = np.zeros([pic_in.Nxg,pic_in.Nyg])
#    ez_slice = np.zeros([pic_in.Nxg,pic_in.Nyg])
#
#    ex_slice[1:-1,:] = (0.5*(slices[1,2:,:] - slices[1,0:-2,:]))/pic_in.Dh
#    ey_slice[:,1:-1] = (0.5*(slices[1,:,2:] - slices[1,:,0:-2]))/pic_in.Dh 
#    ez_slice[:,:] = (0.5*(slices[2,:,:] - slices[0,:,:]))/dz
#    return slices[1,:,:], ex_slice, ey_slice, ez_slice
#
def exact_slice(slices, xg_in, yg_in, zg_in, xg_out, yg_out, z_out):
    slices_array = np.empty([len(xg_in), len(yg_in), len(zg_in)])
    for i in range(len(zg_in)-1):
        slices_array[:,:,i] = slices[i]

    tinterp = ti.Tricubic_Interpolation(A=slices_array, 
                                        x0=xg_in[0],
                                        y0=yg_in[0],
                                        z0=zg_in[0],
                                        dx=xg_in[1]-xg_in[0],
                                        dy=yg_in[1]-yg_in[0],
                                        dz=zg_in[1]-zg_in[0]
                                       )

    the_slice = np.empty([len(xg_out), len(yg_out),8])
    for ii, x_out in enumerate(xg_out):
        for jj, y_out in enumerate(yg_out):
            the_slice[ii,jj,0] = tinterp.val(     x_out, y_out, z_out)
            the_slice[ii,jj,1] = tinterp.ddx(     x_out, y_out, z_out)
            the_slice[ii,jj,2] = tinterp.ddy(     x_out, y_out, z_out)
            the_slice[ii,jj,3] = tinterp.ddz(     x_out, y_out, z_out)
            the_slice[ii,jj,4] = tinterp.ddxdy(   x_out, y_out, z_out)
            the_slice[ii,jj,5] = tinterp.ddxdz(   x_out, y_out, z_out)
            the_slice[ii,jj,6] = tinterp.ddydz(   x_out, y_out, z_out)
            the_slice[ii,jj,7] = tinterp.ddxdydz( x_out, y_out, z_out)
    return the_slice

def fix_phi(phi_slice, Sx, Sy):
    phi_slice[Sx,:,1] = 0
    phi_slice[:,Sy,2] = 0

    phi_slice[Sx,:,4] = 0
    phi_slice[:,Sy,4] = 0

    phi_slice[Sx,:,5] = 0

    phi_slice[:,Sy,6] = 0

    phi_slice[Sx,:,7] = 0
    phi_slice[:,Sy,7] = 0

    return


