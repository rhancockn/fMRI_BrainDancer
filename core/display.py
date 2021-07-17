# %%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created/Last Edited: September 3, 2019

@author: Rajat Kumar
@maintainer: Rajat Kumar
Notes:
Functions to help with display.  

"""

# %% All imports
import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np

# %%
# =============================================================================
# Display EPI slices
# =============================================================================
def display_all_slices(data,volume_number,file=None,subset=None):
    """Displays slices of a NIFTI volume.

    Parameters
    ----------
    data : obj
        A nibabel nifti object to be displayed.
    volume_number : int
        Volume index to be displayed.
    file : str, optional
        If specified, the slice display will be saved to the given filepath. Otherwise, display onscreen.
    subset : list, optional
        If specified, show only specified slices indices. Otherwise, show all slices
    """


    data_read = data
    count = 1
    
    if subset is None:
        subset = np.arange(data_read.header['dim'][3])
    
    n_slices = len(subset)

    fig = plt.figure(figsize=(15,n_slices+10))
    
    for slice_i in subset:
        fig.subplots_adjust(hspace=0, wspace=0.0005)      
        ax = fig.add_subplot(n_slices/4 + 1, 4, count)
        im = plt.imshow(data_read.get_data()[:,:,slice_i,volume_number])
        ax.set_title(f'Slice {slice_i}',fontsize=18)
        ax.set_axis_off()
        count += 1

    if file is None:
        plt.show()
    else:
        plt.savefig(file)
        plt.close()

# %%
# =============================================================================
# Display T2* maps
# =============================================================================
def display_T2starmaps(image,vmax=55,vmin=40,colormap='seismic'):
    
        plt.imshow(np.nan_to_num(image),cmap=colormap,vmax=vmax, vmin=vmin)
        plt.colorbar()
        plt.xticks([])
        plt.yticks([])
        
# %%