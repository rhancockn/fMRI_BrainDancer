# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created/Last Edited: October 21, 2019

@author: Rajat Kumar
@maintainer: Rajat Kumar
Notes:
Script for executing the analysis routine for denoising.

To do:
Include the option to generate outer mask and corresponding time_series.

"""

# %% All imports
import nibabel as nib, pandas as pd, numpy as np, matplotlib.pyplot as plt
from display import display_all_slices
from cartridge import findcartridge,inner_mask, outer_mask, cen_rotation
from quadrants import quadrant_mask_T2
from t2 import T2Star
from epi import phantom_motion, create_mean_slices, simulate_inner, scanner_output
import os

# %%
# =============================================================================
# Main Routine
# =============================================================================


def create_denoising_dataset(epi_path,log_path,acqtimes_path,rot_dir=-1, interactive=True, img_dir=None, slice_indices=None, inner_mask_level=.004):
    """Generates masks and timeseries for analysis.

    Parameters
    ----------
    epi_path : str
        Path to the phantom data.
    log_path : str
        Path to the BrainDancer movement log file.
    acqtimes_path : str
        Path to the slice timing file.
    interactive : bool, optional
        If True (default), prompt user for decision inputs
    img_dir : str, optional
        If specified, displays will be saved to file. Otherwise, display onscreen.
    slice_indices : tuple, optional
        start and end indices of slices to include for processing.
    inner_mask_level : float
        Threshold for segmenting the phantom.
    """

    data_read = nib.load(epi_path)
    
    if interactive:
        display_all_slices(data_read,0)
    
        start = int(input('Enter the first good slice: '))
        end = int(input('Enter the last good slice: '))
    
    # non-interactive. select specified slices
    elif slice_indices is not None:
        start = slice_indices[0]
        end = slice_indices[1]
        slice_img_path = os.path.join(img_dir, 'slices.png')
        display_all_slices(data_read,0, file=slice_img_path, subset=np.arange(start, end+1))
    # non-interactive, but empty slice list
    else:
        raise TypeError('slice_indices cannot be None in non-interactive mode')

    with open(log_path, 'r') as fp:
        line = fp.readline()
        if line.startswith('Sequence'):
            # skip lines
            log = pd.read_csv(log_path, header=2)
        else:
            log = pd.read_csv(log_path)

    motion_time = np.max(log['Tmot'].values)
    acq_times = pd.read_csv(acqtimes_path)
    motionfree = acq_times[acq_times['Time']>motion_time]['Slice'].values
    total_slices = []    
    for i in list(motionfree):
        if start<= i <= end:
            total_slices.append(i)
    print('Selected Slices for Analysis are: ', total_slices)
    
    imask = []
    cen = []
    imask_metrics = []
    center_rotation_all = []
    omask = []
    detect_remove = []
    updated_total_slices = []
    good_slices = []
    for i in range(len(total_slices)):
        if interactive:
            level_img_path = None
            center_img_path = None
        else:
            level_img_path = os.path.join(img_dir, f'contours_{i:03d}.png')    
            center_img_path = os.path.join(img_dir, f'centers_{i:03d}.png')

        img_complete,cy_complete,cx_complete, radii_complete = inner_mask(epi_path,total_slices[i],volume_num=0,lvl=inner_mask_level,rad1=7,rad2=50,step=1, img_path=level_img_path)
        center_rotation  = cen_rotation(epi_path,total_slices[i],img_complete,cy_complete,cx_complete,radii_complete, canny_sgm=1, img_path=center_img_path)

        if interactive:
            detect = int(input('Enter 1 if this slice is good'))
            good_slices.append(detect)
        
        center_rotation_all.append(center_rotation)
        imask.append(img_complete)
        updated_total_slices.append(total_slices[i])

        # TO DO - Include the option to generate outer mask and corresponding time_series, with something like below:
        #out_mask = outer_mask(data_read,findcartridge(data_read,total_slices[i],0),total_slices[i],0)
        #omask.append(out_mask)

    # update good slices
    if not interactive:
        row_med = np.median([x[0] for x in center_rotation_all])
        col_med = np.median([x[1] for x in center_rotation_all])

        for row_cor,col_cor in center_rotation_all:
            if np.all([row_cor <= row_med+1, row_cor >= row_med-1, col_cor <= col_med+1, col_cor >= col_med-1]):
                good_slices.append(1)
            else:
                good_slices.append(0)
    
    print(good_slices)
    print(center_rotation_all)
    center_rotation_all = [x for good,x in zip(good_slices, center_rotation_all) if good==1 ]
    imask = [x for good,x in zip(good_slices, imask) if good==1 ]
    updated_total_slices = [x for good,x in zip(good_slices, updated_total_slices) if good==1 ]
    
    print(good_slices)
    print(center_rotation_all)
    
    if img_dir is not None:
        motion_img = os.path.join(img_dir, 'motion.png')
    else:
        motion_img = None
    positions = phantom_motion(log_path, img_path=motion_img)
    synth = create_mean_slices(data_read,updated_total_slices,imask,200)
    simulated_data = simulate_inner(synth,positions,updated_total_slices,imask,center_rotation_all,rot_dir)
    scanner_inner = scanner_output(data_read,positions,updated_total_slices,imask,200) # add omask in future for outer cylinder
    
    return simulated_data, scanner_inner, imask, center_rotation_all, updated_total_slices


# %%
