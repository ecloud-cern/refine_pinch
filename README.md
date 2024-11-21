# refine_pinch

# Requires:

```
pip install git+https://github.com/kparasch/TricubicInterpolation.git@v1.1.0
```

## Typical usage:

If the file doesn't include the charge density (rho), e.g. when lumping many Inner Triplet ecloud slices to get one effective slice, it must first be generated:

```
python generate_rho.py Pinch.h5
```

```
python reorder_slices.py Pinch.h5
python refine_pinch.py Pinch.h5
```

### arguments in refine_pinch.py

- MTI (magnify tranverse input) refers to ratio between input grid and auxilliary grid
- MLI (magnify longitudinal input) refers to ratio between number of slices in the input grid and and the auxilliary grid
- DTO (demagnify transverse output) refers to ratio between input grid and output grid
- DLO (demagnify longitudinal output) refers to ratio between number of slices in the input grid and and the output grid
