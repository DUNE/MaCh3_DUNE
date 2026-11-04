import uproot
import awkward as ak
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

file_name = 'MCEventInfoForBinningOptimisation.root'

with uproot.open(file_name) as f:
    # Load some data from the file
    df = f['MCEvent'].arrays([
                                    'TrueENu',
                                    'TrueCZ',
                                    'TrueEHad',
                                    'TrueELep',
                                    'RecoENu',
                                    'RecoEHad',
                                    'RecoELep',
                                    'RecoCZ',
                                    'Sample',
                                    'EventWeight'
                                    ],
                                     library='pd')

df['RecoY'] = df['RecoEHad'] / df['RecoENu']

display(df)

def create_energy_bins(cp_step, m32_step, bsm_step):

    delta_cp_bins = np.arange(0., 1., cp_step)
    mass_ordering_bins = np.arange(1, 10, m32_step)
    bsm_bins = np.arange(10, 100+bsm_step, bsm_step)

    E_nu_bins = np.concatenate([delta_cp_bins,mass_ordering_bins,bsm_bins])
    return E_nu_bins

E_nu_bins = create_energy_bins(0.2, 0.05, 0.5)
print(f"{len(E_nu_bins)}: {E_nu_bins}")

def optimise_edges(axis, bad_bins):
    new_bins = histo_axes_digi_edges[axis].copy()
    if (len(bad_bins) == len(new_bins)-2):
        raise RuntimeError("all bins are bad!")
    print(f"modifying: {axis}")

    edges_to_remove=set()

    for bad_bin in bad_bins:
        if bad_bin == 0:
            raise RuntimeError(f"underflow bin index received. should not happen!")
        elif bad_bin == len(histo_axes_digi_edges[axis])-1:
            raise RuntimeError(f"manual overflow bin index received. should not happen!")
        elif bad_bin == len(histo_axes_digi_edges[axis]):
            raise RuntimeError(f"np.digi overflow bin index received. should not happen!")
        else:
            bad_edge_index=bad_bin
            if bad_bin == len(new_bins)-1:
                bad_edge=bad_edge-1

            bad_edge=histo_axes_digi_edges[axis][bad_edge_index]
            #print(f"bin . edge to be removed: {bad_bin} {bad_edge}")
            edges_to_remove.update([bad_edge])
    print(f"edges to remove: {edges_to_remove}")
    for bad_edge in edges_to_remove:
            new_bins=np.delete(new_bins, np.where(new_bins == bad_edge))

    histo_axes_digi_edges[axis]=new_bins
    print(f"new {axis} edges: {histo_axes_digi_edges[axis]}")
###

###
def get_bin_indeces(df, axis):
    mc_counts=[]
    data_counts=[]

    data=df[axis].to_numpy()
    weights=df['EventWeight'].to_numpy()
    edges=histo_axes_digi_edges[axis]

    indeces = np.digitize(data, edges, right=False)

    for bin in range(1,len(edges)):
        mc_c = np.count_nonzero(indeces == bin)
        mc_counts.append(mc_c)

        bin_weigths=weights[indeces == bin]
        data_c=np.sum(bin_weigths)*exposure_scaling
        data_counts.append(data_c)
    return indeces, mc_counts, data_counts
###

###
def check_bins_are_ok_axis(axis, min_mc, min_data, df, idx_sel=None):

    df_sel=None
    if idx_sel is not None:
        df_sel=df.iloc[idx_sel]
    else:
        df_sel=df

    indeces, mc_counts, data_counts = get_bin_indeces(df_sel, axis)

    if idx_sel is not None:
        df.iloc[idx_sel, df.columns.get_loc(axis+'binIndex')]=indeces
    else:
        df[axis+'binIndex']=indeces

    n_bins=len(mc_counts[:-1])

    bad_bins=[]
    bad_bins_data_counts=[]
    bad_bins_mc_counts=[]
    ok=False

    bin_idx = 0
    while bin_idx < n_bins:
        mc_count=mc_counts[bin_idx]
        data_count=data_counts[bin_idx]
        while True:
            if mc_count < min_mc or data_count < min_data:
                #print(f"bad bin: {bin_idx} {mc_count}, {data_count}")

                # last bin, special
                if bin_idx == n_bins-1:
                    #print(f"last bin: {bin_idx}, merging with previous")
                    # actually need to merge with last good bin
                    bad_bins.append(bin_idx)
                    bin_idx += 1
                    break
                else:
                    bin_idx += 1
                    #print(f"merging with {bin_idx}, {mc_counts[bin_idx]}, {data_counts[bin_idx]}")
                    bad_bins.append(bin_idx)
                    mc_count=mc_count+mc_counts[bin_idx]
                    data_count=data_count+data_counts[bin_idx]
                    #print(f"merged {mc_count}, {data_count}")

            else:
                bin_idx += 1
                break
    if (len(bad_bins)==0):
        ok=True

    if not ok:
        print(f"bad bins: {bad_bins}")
    #print(f"mc counts: {mc_counts[:-1]}")
    print(f"data counts: {data_counts[:-1]}")

    return ok, bad_bins, bad_bins_mc_counts, bad_bins_data_counts
###

###
def optimise_bins(axis, hdim_bin, min_mc, min_data, df, idx_sel=None):
    while True:
        ok, bad_bins, bad_bins_mc_counts, bad_bins_data_counts = check_bins_are_ok_axis(axis,
                                                      min_mc,
                                                      min_data,
                                                      df,
                                                      idx_sel)
        if ok:
            print(f"{axis} bin {hdim_bin} ok")
            break
        else:
            print(f"{axis} bin {hdim_bin} to be optimised")
            optimise_edges(axis, bad_bins)
    return True
###

###
def optimise_bins_1d(df, axis):
    ok = optimise_bins(axis, 0, min_mc_events, min_data_events, df, idx_sel=None)
    return ok
###

###
def optimise_bins_2d(df, axis1, axis2):

    dim_1_edges=histo_axes_digi_edges[axis1]

    # ignore overflow bin
    n_bins=len(dim_1_edges)-1

    for bin in range(1,n_bins):

            print("")
            sel_idx = df.index[df[axis1+'binIndex'] == bin]
            print(f"checking {axis1} bin {bin} if {axis2} bins are ok with total events: {len(sel_idx)}")

            ok = optimise_bins(axis2,
                               bin,
                                min_mc_events,
                                min_data_events*dim_2_factor,
                                df,
                                sel_idx)
    print("")
    return True
###

###
def optimise_bins_3d(df, axis1, axis2, axis3):

    dim1_edges=histo_axes_digi_edges[axis1]
    n_bins_dim1=len(dim1_edges)-1

    dim2_edges=histo_axes_digi_edges[axis2]
    n_bins_dim2=len(dim2_edges)-1

    for i in range(1,n_bins_dim1):
            for j in range(1,n_bins_dim2):

                print("")
                sel_idx = df.index[ (df[axis1+'binIndex'] == i) & (df[axis2+'binIndex'] == j) ]
                print(f"checking {axis1},{axis2} bin {i},{j} if {axis3} bins are ok  with total: {len(sel_idx)}")

                ok = optimise_bins(axis3,
                                   j,
                                              min_mc_events,
                                              min_data_events,
                                              df,
                                              sel_idx)
    print("")
    return True
###

# order of axes here determines optimisation order
histo_axes_bins = {
                        'RecoCZ': np.arange(-1., 1.01, 0.2),
                        'RecoENu': E_nu_bins,
                        'RecoY': np.arange(0., 1.01, 0.1),
                        }

axes=list(histo_axes_bins.keys())

# add overflow bin to reproduce numpy histogram behaviour
histo_axes_digi_edges=dict.fromkeys(histo_axes_bins.keys())
for axis, edges in histo_axes_bins.items():
    histo_axes_digi_edges[axis] = np.concatenate([edges, [float('inf')]])
    print(f"Initial {axis} edges: {histo_axes_digi_edges[axis]}")
    df[axis+'binIndex'] = np.nan

histo_axes_are_optimized = dict.fromkeys(histo_axes_bins.keys(), False)

# algo parameters
min_mc_events=100
min_data_events=3
# scale required data counts when opt. 2nd dimension
dim_2_factor=10

# 1. is 400 kt-y
exposure_scaling=1.

print("")

for i, axis in enumerate(histo_axes_bins.keys()):
    dim=i+1
    while True:
        if dim == 1:
            histo_axes_are_optimized[axis] = optimise_bins_1d(df, axes[0])
        elif dim == 2:
            histo_axes_are_optimized[axis] = optimise_bins_2d(df, axes[0], axes[1])
        elif dim == 3:
            histo_axes_are_optimized[axis] = optimise_bins_3d(df, axes[0], axes[1], axes[2])
        else:
            print(f"dim {dim} too high")
            break

        if histo_axes_are_optimized[axis] == True:
            print(f"{axis} done")
            print("")
            break
        else:
            raise RuntimeError(f"{axis} optimisation failed")

print(f"exposure {400*exposure_scaling} kt-y")
for axis, edges in histo_axes_digi_edges.items():
    print(f"Final {axis} edges: {edges[:-1]}")
