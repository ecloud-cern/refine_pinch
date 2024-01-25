import h5py
import sys

#fname = "LHC6.8TeV_v1_Q1R1_0_sey1.35_1.20e11ppb.h5"
fname = sys.argv[1]

AA = h5py.File(fname, "a")
N = len(AA['slices'].keys())

for ii in range(N):
    AA[f'slices/tslice{ii}'] = AA[f'slices/slice{ii}']
    del AA[f'slices/slice{ii}']

for ii in range(N):
    AA[f'slices/slice{N - ii - 1}'] = AA[f'slices/tslice{ii}']
    del AA[f'slices/tslice{ii}']

# AA['grid/zg'][()] = AA['grid/zg'][()][::-1]
# print(AA['slices'].keys())
# print(AA['grid/zg'][()])