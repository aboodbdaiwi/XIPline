# -*- coding: utf-8 -*-
"""
Created on Tue May 24 12:53:28 2022

@author: WESMD7
"""

#%
import numpy as np
import sys
import nibabel as nib 
from scipy.io import savemat, loadmat
from keras.models import load_model
#from tkinter import Tk
#from tkinter.filedialog import askopenfilename
#from tkinter.filedialog import askdirectory

#import pandas as pd
#import os
# folder_path =  os.getcwd()
# script_path = __file__
# print('getcwd:      ', os.getcwd())
# print('__file__:    ', __file__)

# # % test Xe -Unet - Resnet 
# ''' '''
# # Create an instance of Tkinter
# root = Tk()
# # Hide the main window
# root.withdraw()
# # Open the file selection dialog
# file_path = askopenfilename()
# # Print the selected file path
# print("Selected file:", file_path)

# # load mat files
# X_test = loadmat(file_path)
# X_test = X_test["VentImage"]

#%
 



def getLungMask3D(imagearray,model,IMG_SIZE):
    X = np.empty((1,IMG_SIZE,IMG_SIZE,IMG_SIZE,1))
    X[0,:,:,:,0] = imagearray
    y = model.predict(X)
    #y = np.flip(y,2)  #Kept here in case of future troubleshooting
    return y

def save_NIBnifti(file,output):
    preNii_gen_mask = file
    #preNii_gen_mask = np.flip(preNii_gen_mask,0)
    #preNii_gen_mask = np.rot90(preNii_gen_mask)
    img = nib.Nifti1Image(preNii_gen_mask, np.eye(4))
    nib.save(img,output)
    
def Segment3D(InputIMG,modelinput,threshold):
    IMG_SIZE =112
    img_size = (IMG_SIZE,IMG_SIZE,IMG_SIZE)
    
    InputIMG = loadmat(InputIMG)
    imgarray = InputIMG["VentImage"]

    #imgarray=np.empty((img_size))
    #imgarray[:,:,:] = get_NIBnifti(InputIMG);
    #imgarray = pd.read_csv(InputIMG,header=None)
    #imgarray = pd.DataFrame.to_numpy(imgarray)
    #imgarray = imgarray.reshape(112,112,112)
    
    mx = imgarray.max()
    mn = imgarray.min()
    mx-=mn
    imgarray[:,:,:]= ((imgarray[:,:,:]-mn)/mx)
    
    
    MODEL = modelinput
    model = load_model(MODEL,compile=False)
    
    lungmask = getLungMask3D(imgarray,model,IMG_SIZE)
    img = lungmask[0,:,:,:,0]
    
    #img = img > float(threshold)
    for tt in range(IMG_SIZE):
        for uu in range(IMG_SIZE):
            for vv in range(IMG_SIZE):
                if img[tt,uu,vv] >= float(threshold):
                    img[tt,uu,vv] = 1
                else:
                    img[tt,uu,vv] = 0
    return img

if __name__ == '__main__':
    InputImg = (sys.argv[1])    #'VentImage.mat'   #X_test #(sys.argv[1])
    modelinput = (sys.argv[2])  #'D:\Github\HP_Xe_Analysis_App\+Segmentation\AutoSegment_3DGasExchange_Xe_200e.hdf5' #(sys.argv[2])
    thresh = (sys.argv[3])      #0.2  #float(sys.argv[3])
    outname = (sys.argv[4])     #'AutoMask.nii.gz' # (sys.argv[4])
    MaskOut = Segment3D(InputImg,modelinput,thresh)
    save_NIBnifti(MaskOut,outname)
        
