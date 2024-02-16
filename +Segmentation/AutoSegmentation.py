# -*- coding: utf-8 -*-
"""
Created on Tue August 24 12:53:28 2023

@author: ASB
"""

import sys
import numpy as np
from keras.models import load_model
import nibabel as nib
from scipy.io import savemat, loadmat
#import os
#from matplotlib import pyplot


def save_NIBnifti(file,outpath):
    # save .nii file
    img = nib.Nifti1Image(file, np.eye(4))
    nib.save(img,outpath+'AutoMask.nii.gz')
    # save .mat file
    savemat(outpath+'AutoMask.mat', {"AutoMask":file} )
    
def Segment3D(SegmentType):
    inoutImgFolder = 'C:/XIPline/'
    modelFolder = 'C:/XIPline/models/'
    InputImg = 'InputImage.mat'
    imgarray = loadmat(inoutImgFolder+InputImg)
    imgarray = imgarray["Images"]

    X_test = imgarray / np.max(imgarray)
    # % load model
    if SegmentType == 'vent_2D_1ch_cor':
        model = load_model(modelFolder+'AutoSegment_2DVent_Xe_coronal_2000e.hdf5',compile=False) 
    elif SegmentType == 'vent_2D_2ch_cor':
        model = load_model(modelFolder+'AutoSegment_2DVent_Xe_H_coronal_1000e.hdf5',compile=False) 
    elif SegmentType == 'vent_2D_1ch_axi':
        model = load_model(modelFolder+'AutoSegment_2DVent_Xe_axial_2000e.hdf5',compile=False) 
    elif SegmentType == 'diff_2D_1ch_axi':
        model = load_model(modelFolder+'AutoSegment_2DDiff_Xe_axial_2000e.hdf5',compile=False)         
    elif SegmentType == 'gx_3D_1ch_iso':
        model = load_model(modelFolder+'AutoSegment_3DGasExchange_Xe_200e.hdf5',compile=False) 
    elif SegmentType == 'gx_3D_2ch_iso':
        model = load_model(modelFolder+'AutoSegment_3DGasExchange_Xe_H_1000e.hdf5',compile=False) 
            
    #% predict mask for each slice
    if SegmentType == 'vent_2D_1ch_cor' or SegmentType == 'vent_2D_2ch_cor' or SegmentType == 'vent_2D_1ch_axi' or SegmentType == 'diff_2D_1ch_axi':
        gen_masks = np.zeros((X_test.shape[1],X_test.shape[2],X_test.shape[0]))
        for i in range(0, X_test.shape[0]):   
            test_img = X_test[i]
            test_img_input=np.expand_dims(test_img, 0)
            prediction = model.predict(test_img_input)
            gen_masks[:,:,i] = prediction[0,:,:,0]
        gen_masks = gen_masks > 0.9
    elif  SegmentType == 'gx_3D_1ch_iso' or SegmentType == 'gx_3D_2ch_iso':
        IMG_SIZE = 112
        if len(X_test.shape) == 4:
            X_test_tmp = np.zeros((1,X_test.shape[1],X_test.shape[2],X_test.shape[3],1))
            X_test_tmp[:,:,:,:,0] = X_test
            X_test = X_test_tmp
        gen_masks = model.predict(X_test)
        gen_masks = gen_masks > 0.2
        gen_masks = gen_masks[0]
    gen_masks = gen_masks.astype(float)
    lungmask = gen_masks #[0]

    return lungmask

if __name__ == '__main__':
    #InputImg = (sys.argv[1])
    #modelinput = (sys.argv[2])
    SegmentType = (sys.argv[1]) # vent_2D_1ch_cor, vent_2D_2ch_cor, vent_2D_1ch_axi, gx_3D_1ch_iso, gx_3D_2ch_iso
    #denoise = (sys.argv[4])
    #outname = 'AutoMask.nii.gz' #(sys.argv[4])
    MaskOut = Segment3D(SegmentType)
    outpath = 'C:/XIPline/' #os.path.dirname(InputImg)
    save_NIBnifti(MaskOut,outpath)
    
#InputImg = 'Mat2Py_MLvars.csv'
#modelinput = 'D:/Github/Segmentation3D_Xe_model1_1.hdf5'
#threshold = 0.2
#outname = 'AutoMask.nii.gz'
#var = Segment3D(InputImg,modelinput,threshold)
#pyplot.imshow(var[:,:,50])



















