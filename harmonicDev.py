"""En este programa se realiza la modificación del modelo sinusoidal
se puede ver como las armónicas tienen un filtrado y también son ordenadas"""

import numpy as np
import scipy.io.wavfile as waves #para guardar archivos dewav
from scipy.signal import get_window
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'software/models/'))
import utilFunctions as UF
import sineModel as SM
import sineModelModified as SMM #Archivo modificado
import stft

import time


#import matplotlib
#matplotlib.use('Qt5Agg') #con esto mostrara las graficas en ventanas
import matplotlib.pyplot as plt

from scipy.signal import medfilt
import measureSNR as SNR




H = 150
M = 605
N = 1024
maxSines = 20 #Numero de armonicos
t = -58 #threshold para modelo sinusoidal
inputFile = 'sounds/faSax_2.wav'
w = get_window('blackman',M,False)

fs, x = UF.wavread(inputFile)
mX, pX = stft.stftAnal(x, w, N, H)

plt.figure(1, figsize=(9.5, 6))

plt.subplot(311)
numFrames = int(mX[:,0].size)
frmTime = H*np.arange(numFrames)/float(fs)
binFreq = np.arange(N/2+1)*float(fs)/N
plt.pcolormesh(frmTime, binFreq, np.transpose(mX))
plt.title('Spectrogram (%s), M=%d, N=%d, H=%d' % (inputFile,M,N,H))
plt.autoscale(tight=True)



plt.subplot(312)

xtfreq, xtmag, xtphase = SM.sineModelAnal(x, fs, w, N, H, t=-58, maxnSines = 20, minSineDur=0.001, freqDevOffset=10, freqDevSlope=0.01)
numFrames = int(xtfreq[:,0].size)
frmTime = H*np.arange(numFrames)/float(fs)
xtfreq[xtfreq<=0] = np.nan
plt.plot(frmTime, xtfreq)
#y = SM.sineModelSynth(xtfreq, xtmag, xtphase, N, H, fs)

''' 
Del analisis sinusoidal lo que importa es la magnitud y el tiempo en que se desarrollan los armonicos, ya que los
armonicos no son necesarios y se pueden generar.
la funcion sineModelAnal() regresar los valores sin orden por lo que se ordenaran iniciando por la fundamental
y posteriormente las armonicas
'''
plt.subplot(313)
xtfreq3, xtmag3, xtphase3 = SMM.sineModelAnal(x, fs, w, N, H, t=-58, maxnSines = 20, minSineDur=0.01, freqDevOffset=1, freqDevSlope=0.01)
#tfreq4, tmag4, tphase4 = SMM.armonicOrder(xtfreq3, xtmag3, xtphase3) #funcion para ordenar armonicas
numFrames3 = int(xtfreq3[:,0].size)
frmTime3 = H*np.arange(numFrames3)/float(fs)
xtfreq3[xtfreq3<=0] = np.nan
plt.plot(frmTime3, xtfreq3)
plt.show()

#Realizar sintesis
#Sntesis utilizando modelo sinusoidal original
y = SM.sineModelSynth(xtfreq, xtmag, xtphase, 600, 150, fs) #el valor de N debe cuatro veces H
waves.write('sinusoidalModel.wav', fs, y)
SNR = SNR.computeSNR(inputFile, 'sinusoidalModel.wav')
print("SRN modelo sinusoidal= ", SNR)

#Sintesis utilizando modelo sinusoidal modificado
print(xtfreq3.shape)
y3 = SMM.sineModelSynth(xtfreq3, xtmag3, xtphase3, 600, 150, fs) #el valor de N debe cuatro veces H
waves.write('test1.wav', fs, y3)
#time.sleep(5) # espera en segundos
#SNR = SNR.computeSNR(inputFile, 'test1.wav')
#print("SNR Modelo modificado = ", SNR)


#Se filtran las magnitudes de cada una de las armonicas para obtener una mejor visualizacion
xtmag3[xtmag3==0] = np.min(xtmag3)    # if zeros add epsilon to handle log



#En esta seccion se realiza la graficacion de las magnitudes de las armonicas
#plt.figure(1, figsize=(9.5, 6))

#plt.subplot(211)
plt.title('Amplitud de las armónicas sin filtrado')
linesNoFilt = []
for n in range(maxSines):
    linesNoFilt.append(plt.plot(frmTime3, xtmag3[:, n], label='%i' % (n + 1)))  # sin filtrado
plt.legend(loc='upper right')
plt.xlabel("Tiempo")
plt.ylabel("dB")
plt.show()

#plt.subplot(212)
plt.title('Amplitud de las armónicas con filtrado')
linesFilt = []
for n in range(maxSines):
    linesFilt.append(   plt.plot(frmTime3, medfilt(xtmag3[:,n],61), label='%i' %(n + 1))   ) #con filtrado
    xtmag3[:,n] = medfilt(xtmag3[:, n], 61) #aplicar filtro y guardar en la misma variable
plt.legend(loc='upper right')
plt.xlabel("Tiempo")
plt.ylabel("dB")
plt.show()

#Sintesis de audio aplicando filtrado
y4 = SMM.sineModelSynth(xtfreq3, xtmag3, xtphase3, 600, 150, fs)  # el valor de N debe cuatro veces H
waves.write('audioFiltrado.wav', fs, y4)

"""
#graficar solo una sección el audio sintetizado
plt.subplot(211)
plt.plot(x[44100:44200])
plt.subplot(212)
plt.plot(y4[44100:44200])
plt.show()
"""