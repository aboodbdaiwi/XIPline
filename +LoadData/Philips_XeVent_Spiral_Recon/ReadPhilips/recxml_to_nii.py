# -*- coding: utf-8 -*-
"""
Title: 
    recxml_to_nii.py

Agenda: 
    A .py script to convert .rec/.xml files to .nii format.

Author: 
    Joseph Plummer - joseph.plummer@cchmc.org 

Creation date: 
    2023-03-13
    
Modification date: 
    
    
"""

# %% Import packages

# Base
from readphilips.file_io import io
import readphilips as rp
import numpy as np
import nibabel as nib
import os
import matplotlib.pyplot as plt
plt.style.use('dark_background')

# %% Import data

# Data folder location
location = io(
    initial_dir="C:\\Users\\PLUMD2\\Desktop\\Images for Compressed Sensing").selectDirectory()

# Load data
filename = io(initial_dir=location, f_type="rec files",
              ext_type="*.rec").selectFile()
rec = rp.PhilipsData(filename)
rec.compute()
rec_convert = np.array(rec.data)
rec_convert = np.squeeze(rec_convert)
print("Data shape = " + str(np.shape(rec_convert)))

# %% Save images as Nifti files

rec_convert_temp = np.array(rec_convert)
rec_convert = np.rot90(rec_convert_temp, k=1,  axes=[1, 3])
rec_convert = np.flip(rec_convert, axis=1)
rec_convert = np.flip(rec_convert, axis=3)

# Build an array using matrix multiplication
scaling_affine = np.array([[1, 0, 0, 0],
                           [0, 1, 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, 1]])

# Rotate gamma radians about axis i
cos_gamma = np.cos(0)
sin_gamma = np.sin(0)
rotation_affine_1 = np.array([[1, 0, 0, 0],
                              [0, cos_gamma, -sin_gamma,  0],
                              [0, sin_gamma, cos_gamma, 0],
                              [0, 0, 0, 1]])
cos_gamma = np.cos(np.pi)
sin_gamma = np.sin(np.pi)
rotation_affine_2 = np.array([[cos_gamma, 0, sin_gamma, 0],
                              [0, 1, 0, 0],
                              [-sin_gamma, 0, cos_gamma, 0],
                              [0, 0, 0, 1]])
cos_gamma = np.cos(0)
sin_gamma = np.sin(0)
rotation_affine_3 = np.array([[cos_gamma, -sin_gamma, 0, 0],
                              [sin_gamma, cos_gamma, 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]])
rotation_affine = rotation_affine_1.dot(
    rotation_affine_2.dot(rotation_affine_3))

# Apply translation
translation_affine = np.array([[1, 0, 0, 0],
                               [0, 1, 0, 0],
                               [0, 0, 1, 0],
                               [0, 0, 0, 1]])

# Multiply matrices together
aff = translation_affine.dot(rotation_affine.dot(scaling_affine))

# %% Make path to save image results

# Check whether a specified save data path exists
results_exist = os.path.exists(location + "/results")

# Create a new directory because the results path does not exist
if not results_exist:
    os.makedirs(location + "/results")
    print("A new directory inside: " + location +
          " called 'results' has been created.")

# %% Save data as a .nii file depending on data size
if np.ndim(rec_convert) == 4:
    print("Loaded data has 4 dimensions.")
    print("First dimension assumed to be multiple bins used in reconstruction.")
    N_bins = np.shape(rec_convert)[0]
    for i in range(N_bins):
        img = rec_convert[i, ...]
        ni_img = nib.Nifti1Image(abs(img), affine=aff)
        nib.save(ni_img, location + '/results/img_recon_bin_' +
                 str(i) + '.nii.gz')
elif np.ndim(rec_convert) == 3:
    print("Loaded data has 3 dimensions.")
    print("Single reconstruction dataset to be saved.")
    img = rec_convert
    ni_img = nib.Nifti1Image(abs(img), affine=aff)
    nib.save(ni_img, location + '/results/img_recon_full' + '.nii.gz')

# %%
