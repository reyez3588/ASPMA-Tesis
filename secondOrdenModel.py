import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal
#import control
from control.matlab import *  # MATLAB-like functions


def soModel(valStable, valPeak, numSamples, delaySamples, peakValIndex, fs, H, initCondition):
    tPeak = peakValIndex * H / fs
    K = valStable
    d = (valPeak - valStable) / valStable #valor delta
    dampingFactor = math.sqrt(((math.log(d))**2) / ((math.pi**2) + ((math.log(d))**2)))
    freqNatural = (math.pi / tPeak)/ (math.sqrt(1-(dampingFactor**2))) #Frecuencia natural no amortiguada
    w0 = freqNatural
    #Descomentar las lineas siguientes para mostrar los parametros calculados
    '''
    print('El valor de la ganancia es: ', K)
    print('El valor del factor de amortiguamiento es: ', dampingFactor)
    print('El valor de la frecuencia natural no amortiguadaes: ', freqNatural)
    '''

    '''
    Realizar la grafica e la funcion de transferencia mapeado en el plano z
    '''
    sysS = TransferFunction([K*(w0**2)], [1, 2*dampingFactor*w0, w0**2])
    sysD = sample_system(sysS, H/fs, method='bilinear')

    zerosArray = np.zeros(delaySamples)
    onesArray = np.ones(numSamples - delaySamples)
    inputSys = np.concatenate((zerosArray, onesArray), axis=None)
    inputSys[-30:-1] = 0 #para agregar un release al final de la nota
    #print(zArray.size)
    #print(inputSys.size)

    #x = np.arange(numSamples) * H / fs
    #T, yout, xout = step(sysD, T=numSamples*H/fs,X0=initCondition/1., input=inputSys
    yout, T, xout = lsim(sysD, U=inputSys)
    #print("len: ", len(yout))

    '''
    timePlot = np.arange(numSamples-delaySamples) * H / fs
    zArray = np.ones(delaySamples)*initCondition
    #zArray[delaySamples:]= timePlot
    yout, T = step(sysD, timePlot, initCondition)
    x = np.arange(numSamples) * H / fs
    y = np.concatenate((zArray, yout), axis=None)

    '''

    #Descomentar las lineas siguientes para graficar las armonicas parametrizadas

    #plt.plot(x, y)
    #plt.show()
    #plt.grid()
    #print(sysD)

    return T, yout

#x = soModel(valStable=200, valPeak=210, numSamples=1764, delaySamples= 0, peakValIndex=1000,  fs=44100, H=150)