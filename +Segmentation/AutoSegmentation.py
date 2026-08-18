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
import os

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# from matplotlib import pyplot
# Set XIPline root based on operating system
if sys.platform.startswith('win'):
    XIPlineRoot = 'C:/XIPline'
else:
    # macOS / Linux
    XIPlineRoot = os.path.join(os.path.expanduser('~'), 'XIPline')

def save_NIBnifti(file, outpath):
    # Save .nii.gz file
    nii_path = os.path.join(outpath, 'AutoMask.nii.gz')
    img = nib.Nifti1Image(file, np.eye(4))
    nib.save(img, nii_path)
    # Save .mat file
    mat_path = os.path.join(outpath, 'AutoMask.mat')
    savemat(mat_path, {"AutoMask": file})

def Segment3D(SegmentType):
    inoutImgFolder = XIPlineRoot
    modelFolder = os.path.join(XIPlineRoot, 'models')
    InputImg = 'InputImage.mat'
    InputImgPath = os.path.join(inoutImgFolder, InputImg)
    imgarray = loadmat(InputImgPath)
    imgarray = imgarray["Images"]
    X_test = imgarray / np.max(imgarray)
    # Load model
    if SegmentType == 'vent_2D_1ch_cor':
        # model = load_model(
        #     os.path.join(modelFolder, '2DVent_Xe_coronal_1000e_20250509.hdf5'),
        #     compile=False
        # )
        # model = load_model(
        #     os.path.join(modelFolder, '2DVent_XeCTC_20250603_1000epochs.hdf5'),
        #     compile=False
        # )
        model = load_model(
            os.path.join(modelFolder, '2DVent_XeCTC_20251003_1000epochs.hdf5'),
            compile=False
        )
    elif SegmentType == 'vent_2D_2ch_cor':
        model = load_model(
            os.path.join(modelFolder, '2DVent_Xe_H_coronal_1000e_20230528.hdf5'),
            compile=False
        )
    elif SegmentType == 'vent_2D_1ch_axi':
        model = load_model(
            os.path.join(modelFolder, '2DVent_Xe_axial_1000e_20250509.hdf5'),
            compile=False
        )
    elif SegmentType == 'vent_anat_2D_1ch_cor':
        model = load_model(
            os.path.join(modelFolder, '2DVent_H_coronal_2000e_20230818.hdf5'),
            compile=False
        )
    elif SegmentType == 'diff_2D_1ch':
        # model = load_model(
        #     os.path.join(modelFolder, '2DDiff_Xe_axial_2000e_20240118.hdf5'),
        #     compile=False
        # )
        model = load_model(
            os.path.join(modelFolder, '2DDiff_XeCTC_20251020_1000epochs.hdf5'),
            compile=False
        )
    elif SegmentType == 'gx_3D_1ch_iso':
        # model = load_model(
        #     os.path.join(modelFolder, '3DGasExchange_Xe_100e_20250324.hdf5'),
        #     compile=False
        # )
        model = load_model(
            os.path.join(modelFolder, '3DGasExchange_Xe_200e_20230623.hdf5'),
            compile=False
        )
    elif SegmentType == 'gx_3D_2ch_iso':
        # model = load_model(
        #     os.path.join(modelFolder, '3DGasExchange_Xe_HLR_100e_20250324.hdf5'),
        #     compile=False
        # )
        model = load_model(
            os.path.join(modelFolder, '3DGasExchange_Xe_HLR_1000e_20230623.hdf5'),
            compile=False
        )
    else:
        raise ValueError(f'Unknown SegmentType: {SegmentType}')
    # Predict mask
    if SegmentType in [
        'vent_2D_1ch_cor',
        'vent_2D_2ch_cor',
        'vent_2D_1ch_axi',
        'vent_anat_2D_1ch_cor',
        'diff_2D_1ch'
    ]:
        gen_masks = np.zeros(
            (X_test.shape[1], X_test.shape[2], X_test.shape[0])
        )
        for i in range(X_test.shape[0]):
            test_img = X_test[i]
            test_img_input = np.expand_dims(test_img, 0)
            prediction = model.predict(test_img_input)
            gen_masks[:, :, i] = prediction[0, :, :, 0]
        gen_masks = gen_masks > 0.95
    elif SegmentType in ['gx_3D_1ch_iso', 'gx_3D_2ch_iso']:
        if len(X_test.shape) == 4:
            X_test_tmp = np.zeros(
                (
                    1,
                    X_test.shape[1],
                    X_test.shape[2],
                    X_test.shape[3],
                    1
                )
            )
            X_test_tmp[:, :, :, :, 0] = X_test
            X_test = X_test_tmp
        gen_masks = model.predict(X_test)
        gen_masks = gen_masks > 0.6
        gen_masks = gen_masks[0]
    gen_masks = gen_masks.astype(float)
    lungmask = gen_masks
    return lungmask

if __name__ == '__main__':
    # Segment types:
    # vent_2D_1ch_cor
    # vent_2D_2ch_cor
    # vent_2D_1ch_axi
    # vent_anat_2D_1ch_cor
    # diff_2D_1ch
    # gx_3D_1ch_iso
    # gx_3D_2ch_iso
    SegmentType = sys.argv[1]
    MaskOut = Segment3D(SegmentType)
    save_NIBnifti(MaskOut, XIPlineRoot)
    
'''
SegmentType = 'diff_2D_1ch'
MaskOut = Segment3D(SegmentType)
save_NIBnifti(MaskOut, XIPlineRoot)
'''