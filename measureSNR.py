import os
import sys
import numpy as np
import math
from scipy.signal import get_window
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import sineModelModified as SMM
import utilFunctions as UF
eps = np.finfo(float).eps


"""
With the input signal and the obtained output, compute two different SNR values for the following cases:

1) SNR1: Over the entire length of the input and the output signals.
2) SNR2: For the segment of the signals left after discarding M samples from both the start and the 
end, where M is the analysis window length. Note that this computation is done after STFT analysis 
and synthesis.

The input arguments to the function are the wav file name including the path (inputFile), window 
type (window), window length (M), FFT size (N), and hop size (H). The function should return a python 
tuple of both the SNR values in decibels: (SNR1, SNR2). Both SNR1 and SNR2 are float values. 

"""

def computeSNR(originalFile, rebuildFile):
    """
    Input:
            inputFile (string): wav file name including the path 
            window (string): analysis window type (choice of rectangular, triangular, hanning, hamming, 
                    blackman, blackmanharris)
            M (integer): analysis window length (odd positive integer)
            N (integer): fft size (power of two, > M)
            H (integer): hop size for the stft computation
    Output:
            The function should return a python tuple of both the SNR values (SNR1, SNR2)
            SNR1 and SNR2 are floats.
    """
    ## your code here
    fs, x = UF.wavread(originalFile)
    #y = SMM.sineModel(x, fs, w, N, t)
    fs, y = UF.wavread(rebuildFile)

    #Compute SNR1
    energySignal = math.fsum((abs(x[0:264000]))**2)
    energyNoise = math.fsum((abs(x[0:264000]-y[0:264000]))**2)
    SNR1 = 10*np.log10(energySignal / energyNoise) 

    '''
    #Compute SNR2
    energySignal2 = math.fsum((abs(x[M:-M]))**2)
    energyNoise2 = math.fsum((abs(x[M:-M]-y[M:-M]))**2)
    SNR2 = 10*np.log10(energySignal2 // energyNoise2)
    '''
    return SNR1
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
