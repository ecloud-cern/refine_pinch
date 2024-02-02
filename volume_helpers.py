import numpy as np

def pint(i):
    return (1./(i+1) if i >= 0 else 0)

def vol(A,B):
    return np.matmul(np.matmul(B,PP),A)

PP = np.empty([64,64])
for i in range(4):
    for j in range(4):
        for k in range(4):
            for l in range(4):
                for m in range(4):
                    for n in range(4):
                        PP[i+4*j+16*k,l+4*m+16*n] = pint(i+l)*pint(j+m)*pint(k+n)


def var_x(ix, iy, iz, TI, TIdx, dx):

    coef1 = TI.get_coefs(TI.construct_b(ix,iy,iz))
    a_coef = np.empty_like(coef1)
    for k in range(4):
        for j in range(4):
            for i in range(3):
                a_coef[i+4*j+16*k] = (i+1)*coef1[(i+1)+4*j+16*k]/dx
            a_coef[3+4*j+16*k] = 0

    coef2 = TIdx.get_coefs(TIdx.construct_b(ix,iy,iz))
    v1 = vol(a_coef,a_coef)
    v2 = vol(a_coef,coef2)
    v3 = vol(coef2,coef2)
    vv = (v1 - 2*v2 + v3)
    return vv/v3, vv

def var_y(ix, iy, iz, TI, TIdy, dy):

    coef1 = TI.get_coefs(TI.construct_b(ix,iy,iz))
    a_coef = np.empty_like(coef1)
    for k in range(4):
        for i in range(4):
            for j in range(3):
                a_coef[i+4*j+16*k] = (j+1)*coef1[i+4*(j+1)+16*k]/dy
            a_coef[i+4*3+16*k] = 0

    coef2 = TIdy.get_coefs(TIdy.construct_b(ix,iy,iz))
    v1 = vol(a_coef,a_coef)
    v2 = vol(a_coef,coef2)
    v3 = vol(coef2,coef2)
    vv = (v1 - 2*v2 + v3)
    return vv/v3, vv

def var_z(ix, iy, iz, TI, TIdz, dz):

    coef1 = TI.get_coefs(TI.construct_b(ix,iy,iz))
    a_coef = np.empty_like(coef1)
    for i in range(4):
        for j in range(4):
            for k in range(3):
                a_coef[i+4*j+16*k] = (k+1)*coef1[i+4*j+16*(k+1)]/dz
            a_coef[i+4*j+16*3] = 0

    coef2 = TIdz.get_coefs(TIdz.construct_b(ix,iy,iz))
    v1 = vol(a_coef,a_coef)
    v2 = vol(a_coef,coef2)
    v3 = vol(coef2,coef2)
    vv = (v1 - 2*v2 + v3)
    return vv/v3, vv



