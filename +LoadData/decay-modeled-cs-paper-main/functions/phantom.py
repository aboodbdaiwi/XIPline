# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 11:04:55 2022

@author: WESMD7
"""

import numpy as np
import matplotlib.pyplot as plt


def lungphant_2d(ImSize):
    """_summary_

    Args:
        ImSize (int): Number of voxels along length of image.

    Returns:
        Array: Hyperpolarized xenon lung phantom.
    """

    cent1x = 5*ImSize/16
    cent2x = 11*ImSize/16
    centy = ImSize/2

    x = np.array([i for i in range(ImSize)])+1
    y = np.array([i for i in range(ImSize)])+1
    X, Y = np.meshgrid(x, y)

    # Make 2 elipses to be lungs - we'll make one have sharply varying signal
    # intensities, and one to have smoothly varying signal intensities

    el1 = (X-cent1x)**2/(3*ImSize/16)**2 + (Y-centy)**2/(6*ImSize/16)**2
    # el1 = (X-cent1x)**2/(3*ImSize/32)**2 + (Y-centy)**2/(6*ImSize/32)**2
    # el1 = (X-cent1x)**2/(4*ImSize/16)**2 + (Y-centy)**2/(6*ImSize/16)**2

    el1[el1 <= 1] = 1
    el1[el1 > 1] = 0

    el2 = (X-cent2x)**2/(3*ImSize/16)**2 + (Y-centy)**2/(6*ImSize/16)**2
    # el2 = (X-cent2x)**2/(3*ImSize/32)**2 + (Y-centy)**2/(6*ImSize/32)**2
    # el2 = (X-cent2x)**2/(4*ImSize/16)**2 + (Y-centy)**2/(6*ImSize/16)**2

    el2[el2 <= 1] = 1
    el2[el2 > 1] = 0

    # Make left lung sharply varying:
    sharp = np.ones([ImSize, ImSize])
    sharp[np.where(X < cent1x)] = 0.5
    sharp[np.where(Y < ImSize/4)] = sharp[np.where(Y < ImSize/4)]*0.4

    sharp[np.where(Y >= ImSize/4)] = sharp[np.where(Y >= ImSize/4)]*0.6
    sharp[np.where(Y >= ImSize/2)] = sharp[np.where(Y >= ImSize/2)]*0.8/0.6
    sharp[np.where(Y >= ImSize/4*3)] = sharp[np.where(Y >= ImSize/4*3)]/0.8

    el1 = el1*sharp

    # Make right lung smoothly varying:
    sigx = ImSize/8
    sigy = ImSize/4

    # smooth = np.exp(-((X-cent2x)**2/(2*sigx**2) + (Y-centy)**2/(2*sigy**2)))
    smooth = np.exp(-((X-cent2x)**2/(4*sigx**2) + (Y-centy)**2/(2*sigy**2)))

    # figure; imagesc(smooth)
    el2 = el2*smooth

    im = el1 + el2
    return im
