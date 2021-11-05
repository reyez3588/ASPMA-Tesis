import numpy as np
import scipy.io.wavfile as waves #para guardar archivos de wav
from scipy.signal import get_window
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'software/models/'))
import utilFunctions as UF
import sineModel as SM
import sineModelModified as SMM #Archivo modificado
import measureSNR as SNR # Para calculo de SNR archivo propio
import secondOrdenModel as sOM
import stft

import matplotlib.pyplot as plt

from scipy.signal import medfilt #para poder utilizar el filtro de media


H = 150
M = 605
N = 1024
maxSines = 20 #Numero de armonicos
t = -58 #threshold para modelo sinusoidal
inputFile = 'sounds/faSax_2.wav'
w = get_window('blackman',M,False)

#Aplicar modelo sinusoidal modificado
fs, x = UF.wavread(inputFile)
xtfreq3, xtmag3, xtphase3 = SMM.sineModelAnal(x, fs, w, N, H, t=-58, maxnSines = 20, minSineDur=0.01, freqDevOffset=1, freqDevSlope=0.01)
numFrames3 = int(xtfreq3[:,0].size)
frmTime3 = H*np.arange(numFrames3)/float(fs)
print('numero de muestras de la señal: ', frmTime3.size)
#xtfreq3[xtfreq3<=0] = np.nan #Comentar esto para obtener un audio de salida de lo contrario solo se tiene ruido, con esto la grafica se ve mucho mejor
plt.plot(frmTime3, xtfreq3)
plt.title('Modelo Sinosoidal Modificado (%s), M=%d, N=%d, H=%d' % (inputFile,M,N,H))
plt.show()


'''
#Realizar sintesis
#Sintesis utilizando modelo sinusoidal modificado
y3 = SMM.sineModelSynth(xtfreq3, xtmag3, xtphase3, 600, 150, fs) #el valor de N debe cuatro veces H
waves.write('test1.wav', fs, y3)
#time.sleep(5) # espera en segundos
#SNR = SNR.computeSNR(inputFile, 'test1.wav')
#print("SNR Modelo modificado = ", SNR)
'''



#Se filtran las magnitudes de cada una de las armonicas para obtener una mejor visualizacion
plt.subplot(211)
plt.ylabel("dB")
xtmag3[xtmag3==0] = np.min(xtmag3)    # if zeros add epsilon to handle log
plt.title('Amplitud de las armónicas con filtrado')
linesFilt = []
for n in range(maxSines):
    linesFilt.append(   plt.plot(frmTime3, medfilt(xtmag3[:,n],61), label='%i' %(n + 1))   ) #con filtrado
    xtmag3[:,n] = medfilt(xtmag3[:, n], 61) #aplicar filtro y guardar en la misma variable
#plt.legend(loc='upper right')

#plt.show()



"""
#Graficar las fases para ver su comportamiento
plt.title('Comportamiento de las fases')
linesFilt = []
for n in range(maxSines):
    linesFilt.append(   plt.plot(frmTime3, xtphase3[:,n], label='%i' %(n + 1))   )
    #plt.show()
#plt.legend(loc='upper right')
plt.show()
"""


#PArametrizacion de las armonicas
plt.subplot(212)
plt.xlabel("Tiempo")
plt.ylabel("dB")
matrixDBY = np.zeros(xtmag3.shape) #crear matriz para guardar las harmonicas parametrizadas
plt.title('Armonicas parametrizadas')
paramList = np.zeros((20,5))
linesFiltParam = []
for n in range(maxSines):
    harmNumber = np.array(10 ** ((xtmag3[:, n]) / 20))  # Convertir armonica a escala lineal
    harmSize = harmNumber.size # el numero de muestras
    harmMedia = np.median(harmNumber) # el promedio de la armonica

    #Calcular el valor maximo y su respectivo indice
    harmIndexMax = 0
    while True:
        if harmNumber[harmIndexMax + 1] < harmNumber[harmIndexMax]:
            break
        harmIndexMax = harmIndexMax + 1
    harmMax = harmNumber[harmIndexMax]

    if harmMedia > harmMax: #Para evitar la division por cero en el calculo de la armonica parametrizada
        harmMax = harmMedia*1.3

    if harmMax > 1.5*harmMedia: # si el valor pico es grande comparado con la media entonces hay que reducirlo
        harmMax = harmMedia * 1.3

    #calcular la cantidad de muestras para el delay
    harmDelaySamples = 0
    while True:
        if harmNumber[harmDelaySamples + 1] > harmNumber[harmDelaySamples]:
            break
        harmDelaySamples = harmDelaySamples + 1
    #plt.plot(frmTime3, xtmag3[:, 0])

    
    #print('armónico: ', n+1)
    #print('Media: ', harmMedia)
    #print('Valor pico: ', harmMax)
    #print('Numero de delay: ', harmDelaySamples)
    #print('Indice del valor pico: ', harmIndexMax)

    if n == 0:
        delayFund = harmDelaySamples

    paramList[n, 0] = harmMedia
    paramList[n, 1] = harmMax
    paramList[n, 2] = harmSize
    paramList[n, 3] = abs(harmDelaySamples - delayFund)
    paramList[n, 4] = abs(harmIndexMax - delayFund)

    x, y = sOM.soModel(valStable=harmMedia, valPeak=harmMax, numSamples=harmSize, delaySamples=harmDelaySamples,
                       peakValIndex=harmIndexMax, fs=44100.0, H=150, initCondition=harmNumber[0])
    dbY = 20 * np.log10(np.abs(y))
    #dbY[dbY < -60] = -60  #Comentar esta linea solo es para mostrar el graficado de mejor forma

    linesFiltParam.append(   plt.plot(x, dbY, label='%i' %(n + 1))   ) #Crear etiquetas para cada armonica
    matrixDBY[:,n] = dbY #guardar las armonicas parametrizadas en la matriz correspondiente

#plt.legend(loc='upper right')
plt.show()
np.savetxt('paramList.txt', paramList)
#a = np.loadtxt('datos.txt') #para cargar datos desde un archivo

'''
#Realizar la comparacion de la armonica generada con la original
nHarmonicPrint = 3
harmNumber = np.array(10**((xtmag3[:,nHarmonicPrint])/20)) #Convertir armonica a escala lineal
harmMedia = np.median(harmNumber)
harmMax = np.max(harmNumber)
harmIndexMax = np.argmax(harmNumber)
harmSize = harmNumber.size
harmDelaySamples = 0
while True:
    if harmNumber[harmDelaySamples + 1] > harmNumber[harmDelaySamples]:
        break
    harmDelaySamples = harmDelaySamples + 1

plt.plot(frmTime3, xtmag3[:,nHarmonicPrint])
x,y = sOM.soModel(valStable=harmMedia, valPeak=harmMax, numSamples=harmSize, delaySamples= harmDelaySamples, peakValIndex=harmIndexMax, fs=44100.0, H=150, initCondition=harmNumber[0])
y[y<np.finfo(float).eps] = np.finfo(float).eps    # if zeros add epsilon to handle log
dbY = 20 * np.log10(np.abs(y))
#dbY[dbY<-60] = -60
plt.plot(x,dbY)
plt.scatter(harmDelaySamples*H/fs,-60)
plt.show()



'''
#Sintesis de audio aplicando filtrado
y4 = SMM.sineModelSynth(xtfreq3, xtmag3, np.array([]), 600, 150, fs)  # el valor de N debe cuatro veces H
waves.write('audioFiltrado.wav', fs, y4)

#Sintesi de audio con las armonicas parametrizadas
zeroPhase=np.zeros(xtphase3.shape)
y5 = SMM.sineModelSynth(xtfreq3, matrixDBY, np.array([]), 600, 150, fs)  # el valor de N debe cuatro veces H
waves.write('harmonicParam.wav', fs, y5)
#np.savetxt('freq.txt', matrixDBY)

