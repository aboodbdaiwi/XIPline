
# -*- coding: utf-8 -*-
"""
Title: 
    predict_mask.py
    
Agenda: 
    Functions to read xenon images and predict mask using pretrained models

Author: 
    Abdullah S. Bdaiwi - work: abdullah.bdaiwi@cchmc.org 
                       - personal: abdiawi89@gmail.com
    
Creation date: 
    June 21 2023
    
Modification date: 

    
"""

# % load packeges 
import numpy as np
from matplotlib import pyplot as plt
import os
import random
from skimage.io import imread, imshow
import nibabel as nib 
from scipy.io import savemat, loadmat
from keras.models import load_model
from tkinter import Tk
from tkinter.filedialog import askopenfilename

folder_path =  os.getcwd()
script_path = __file__
print('getcwd:      ', os.getcwd())
print('__file__:    ', __file__)

# Create an instance of Tkinter
root = Tk()
# Hide the main window
root.withdraw()
# Open the file selection dialog
file_path = askopenfilename()
# Print the selected file path
print("Selected file:", file_path)

# load mat files
X_test = loadmat(file_path)
X_test = X_test["Images"]
#X_test = X_test/np.max(X_test)

Vent_wout_H = False
Vent_w_H = False
GaxExchnage_wout_H = False
GaxExchnage_w_H = True

# % load model
if Vent_wout_H:
    model = load_model(folder_path+'\AutoSegment_2DVent_Xe_2000e.hdf5',compile=False) 
elif Vent_w_H:
    model = load_model(folder_path+'\AutoSegment_2DVent_Xe_H_1000e.hdf5',compile=False) 
elif GaxExchnage_wout_H:
    model = load_model(folder_path+'\\AutoSegment_3DGasExchange_Xe_200e.hdf5',compile=False) 
elif GaxExchnage_w_H:
    model = load_model(folder_path+'\\AutoSegment_3DGasExchange_Xe_H_1000e.hdf5',compile=False) 
        
#% predict mask for each slice
if Vent_wout_H or Vent_w_H:
    gen_masks = np.zeros((X_test.shape[1],X_test.shape[2],X_test.shape[0]))
    for i in range(0, X_test.shape[0]):   
        test_img = X_test[i]
        test_img_input=np.expand_dims(test_img, 0)
        prediction = model.predict(test_img_input)
        gen_masks[:,:,i] = prediction[0,:,:,0]
    gen_masks = gen_masks > 0.9
elif GaxExchnage_wout_H or GaxExchnage_w_H:
    IMG_SIZE = 112
    if len(X_test.shape) == 4:
        X_test_tmp = np.zeros((1,X_test.shape[1],X_test.shape[2],X_test.shape[3],1))
        X_test_tmp[:,:,:,:,0] = X_test
        X_test = X_test_tmp
    gen_masks = model.predict(X_test)
    gen_masks = gen_masks > 0.2
gen_masks = gen_masks.astype(float)

# # view images
# test_img_number = random.randint(0, len(X_test)-1)
# plt.figure(figsize=(16, 8))
# plt.subplot(231)
# plt.title('Testing Image')
# plt.imshow(X_test[0,:,:,test_img_number,0], cmap='gray')
# plt.subplot(232)
# plt.title('Prediction on test image')
# plt.imshow(gen_masks[0,:,:,test_img_number,0], cmap='gray')
# plt.show()


save_path = os.path.dirname(file_path)

# save .mat file
savemat(save_path+'/auto_segmented_mask.mat', {"auto_segmented_mask":gen_masks} )
# save .nii file
if Vent_wout_H or Vent_w_H:
    preNii_gen_mask = gen_masks
elif GaxExchnage_wout_H or GaxExchnage_w_H:
    preNii_gen_mask = gen_masks[0,:,:,:,0]
    
preNii_gen_mask = np.flip(preNii_gen_mask,0)
preNii_gen_mask = np.rot90(preNii_gen_mask)
final_generated_mask = nib.Nifti1Image(preNii_gen_mask.astype(np.uint8), affine=np.eye(4))
nib.save(final_generated_mask, save_path+'/auto_segmented_mask.nii.gz') # Here you put the path + the extionsion 'nii' or 'nii.gz'

#% Create montages
''' 
# Create the Xe montage
hpxe_Xe_montage = np.concatenate(X_test[:,:,:,0], axis=1)
# Display the montage
plt.imshow(hpxe_Xe_montage, cmap='gray')
plt.axis('off')
plt.title('Xenon image')
plt.show()

# Create the H montage
hpxe_H_montage = np.concatenate(X_test[:,:,:,1], axis=1)
# Display the montage
plt.imshow(hpxe_H_montage, cmap='gray')
plt.axis('off')
plt.title('Proton image')
plt.show()

# Create the generated mask montage
temp_masks = np.expand_dims(gen_masks, 3)
hpxe_genM_montage = np.concatenate(temp_masks[:,:,:,0], axis=1)
# Display the montage
plt.imshow(hpxe_genM_montage, cmap='gray')
plt.axis('off')
plt.title('Generated Mask')
plt.show()
'''

'''
create .exe file to run in matlab 
•	compile your .exe file.
    o	To create the exe file:
        	pip install pyinstaller
        	In terminal, run:
            pyinstaller --onefile --name="predict_mask_2DVent_wout_H.exe" .\predict_mask.py
This will create a single .exe file under \dist. 

'''