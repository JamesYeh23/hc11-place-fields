# hc-11 `sessInfo.mat` schema

**Status: not yet written.** This document is filled in during step 2, by walking the
HDF5 tree of an actual session file — *not* by copying field names out of the CRCNS data
description PDF. Anything the PDF claims that the file contradicts gets recorded here and
raised before the loader is written around it.

## Method

The `.mat` files are MATLAB v7.3, which is HDF5 underneath, so `scipy.io.loadmat` cannot
read them; `h5py` is used instead. Two things to check for every field:

- **Orientation.** HDF5 stores MATLAB arrays transposed. An `N × 2` position array in
  MATLAB arrives as `2 × N` in `h5py`.
- **References.** MATLAB cell arrays and structs-of-arrays come through as
  `h5py.Reference` objects that must be dereferenced.

## Fields

_To be populated in step 2: full path in the HDF5 tree, dtype, shape, units, value range,
and what the field actually appears to contain._

| Path | dtype | Shape (h5py) | Shape (logical) | Units | Notes |
|------|-------|--------------|-----------------|-------|-------|

## Discrepancies with the data description PDF

_To be populated in step 2._
