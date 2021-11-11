import numpy as np
import matplotlib.pyplot as plt
import math
import secondOrdenModel as soM
import scipy.io.wavfile as waves #para guardar archivos de wav
import sys, os
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'software/models/'))
import utilFunctions as UF
import sineModelModified as SMM #modelo sinusoidal modificado
import harmonicModel as HM
import onsetDetect as ODF
from scipy.signal import get_window
import stft as STFT
from matplotlib.legend_handler import HandlerLine2D


def createTone(note, duration, freqSample, hopSize):
    totalSamples = int(math.floor(duration*freqSample/hopSize)) #Convertir la duracion a samples
    #totalSamples = duration  # duration esta dando en samples por lo que no se requiere convertir
    harmonicParam = np.loadtxt('paramList.txt')

    maxSines = harmonicParam.shape[0] #obtener el numero de armonicas del archivo

    #En esta seccion se crea la matriz con las frecuencias de las armonicas
    '''
    El rango del sax alto es de C#3-A5 que son los limites superiores minimos y maximos con respecto al piano
    los numeros de las teclas correspondientes al piano son 29-61, la primera tecla es C=27.5HZ
    formula = 27.5*math.pow(math.pow(2,1/12),note - 1) fuente: http://elclubdelautodidacta.es/wp/2012/08/calculo-de-la-frecuencia-de-nuestras-notas-musicales/
    https://es.wikipedia.org/wiki/Frecuencias_de_afinaci%C3%B3n_del_piano
    '''
    #freqFund = 27.5*math.pow(math.pow(2,1/12),note - 1) #si se da el numero de nota MIDI utilizar esta linea
    freqFund = note #utilizar esta nota si lo que brinda es una frecuencia
    freqMatrix = np.zeros((totalSamples ,maxSines))  # crear matriz para guardar las harmonicas
    for n in range(maxSines):
        freqMatrix[:, n] = freqFund*(n + 1)
    #np.savetxt('freqMatrix.txt', freqMatrix)

    #En la siguiente seccion se crean las magnitudes
    magnitudDB = np.zeros((totalSamples ,maxSines))  # crear matriz para guardar las harmonicas
    #linesFiltParam = []
    for n in range(maxSines):
        #if n == 3:
        #print(harmonicParam[n,0], harmonicParam[n,1], totalSamples, int(harmonicParam[n,3]), int(harmonicParam[n,4]))
        x, y = soM.soModel(valStable=harmonicParam[n,0], valPeak=harmonicParam[n,1], numSamples=totalSamples, delaySamples=int(harmonicParam[n,3]),
                           peakValIndex=int(harmonicParam[n,4]), fs=freqSample, H=hopSize, initCondition=0)
        y[y < np.finfo(float).eps] = np.finfo(float).eps #para evitar la division entre 0
        dbY = 20 * np.log10(np.abs(y))
        #linesFiltParam.append(plt.plot(x, dbY, label='%i' % (n + 1)))  # Crear etiquetas para cada armonica
        #print(len(dbY))
        #print(magnitudDB.shape)
        if len(dbY) != magnitudDB.shape[0]:
            magnitudDB = np.delete(magnitudDB, -1, 0)
            freqMatrix = np.delete(freqMatrix, -1, 0)
        #print(magnitudDB.shape)
        magnitudDB[:, n] = np.transpose(dbY)
        #magnitudDB[0:len(dbY), n] = np.transpose(dbY)  # guardar las armonicas parametrizadas en la matriz correspondiente

    #plt.legend(loc='upper right')
    #plt.show()
    zeroPhase = np.zeros(freqMatrix.shape) #creacion de la matriz con ceros para las fases
    '''
    #Realizar la sintesis con las frecuencias y magnitudes obtenidas
    y = SMM.sineModelSynth(freqMatrix, magnitudDB, zeroPhase, 600, 150, fs)  # el valor de N debe cuatro veces H
    #waves.write('reCreate.wav', fs, y)
    '''
    return freqMatrix, magnitudDB, zeroPhase

def shapeMatrix(listOnset, freqSample, hopSize, H=150):
    nRaws = 0
    for n in range(listOnset.size - 1):
        duration = (listOnset[n + 1] - listOnset[n]) * hopSize / fs
        totalSamples = int(math.floor(duration * freqSample / H))  # Convertir la duracion a samples
        nRaws = nRaws + totalSamples
    return nRaws

def graphOnset(ODF):
    # Seccion para graficar onset
    numFrames = int(ODF[:, 0].size)
    frmTime = H * np.arange(numFrames) / float(fs)
    binFreq = np.arange(N / 2 + 1) * float(fs) / N
    # line1, = plt.plot(frmTime, ODF[:,0], label='Onset Detection Low Frequency')
    line2, = plt.plot(frmTime, ODF[:, 1], label='Onset Detection High Frequency')
    plt.legend(handler_map={line2: HandlerLine2D(numpoints=4)})
    plt.title('Onset Detection Function, M=%d, N=%d, H=%d' % (M, N, H))
    plt.autoscale(tight=True)
    plt.tight_layout()
    #plt.show()
    return


#Parametros de entrada
"""
inputFile = '../../sounds/cello-phrase.wav'
fs = 44100;
M = 1501 # debe ser un numero impar
w = np.blackman(M)
N = 2048 #2048
H = 128 #128
"""
inputFile = 'sounds/piano.wav'
fs = 44100;
M = 1501 # debe ser un numero impar
w = np.blackman(M)
N = 2048 #2048
H = 128 #128

'''
window= 'blackman'
M = 512
N = 1024
H = 64
w = get_window(window,M,False)

'''


#ANALISIS ARMONICO
#Detectar la fundamental utilizando Two way mismatch
(fs, x) = UF.wavread(inputFile)
t = -65 #-90
minf0 = 100
maxf0 = 800 #300
f0et = 1
maxnpeaksTwm = 4
f0 = HM.f0Detection(x, fs, w, N, H, t, minf0, maxf0, f0et)
f0 = UF.cleaningTrack(f0, 5)
copyF0 = f0.copy() #crer copia de f0 para utilizarse despues
f0[f0 == 0] = np.nan
mX, pX = STFT.stftAnal(x, w, N, H)
numFrames = int(mX[:,0].size)
frmTime = H*np.arange(numFrames)/float(fs)
print("Two way mistach length: ", f0.shape[0])
line1, = plt.plot(frmTime, f0, label='Tracking Two-Way Mismatch')
plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)})
#plt.show()


#obtener ONSET
factorODF = 20 #Factor para que al graficar el ODF con el tracking se vean en la misma escala
window= 'blackman'
ODF = ODF.computeODF(inputFile, window, M, N, H)
print("Maximo valor ODF: ", max(ODF[1][:])) #Este valor se itilizara para determinar el thresholdODF
print("Onset Detection length: ", ODF.shape[0])

#Codigo para recrear notas
thresholdODF = 1 #threslhold para deteccion de peaks
mainAux = []
fault = False
for n in range(ODF.shape[0]): #iterar en los frames para encontrar los onsets
    if ODF[n,1]>=thresholdODF:
        print(n)
        dataAux = []
        dataAux.append(n)#guardar indice de peak detectado
        m = np.copy(n)
        while copyF0[m] == 0: #seguir hasta que la frecuencia sea diferente de 0
            m = m + 1
            if m == len(copyF0):
                fault = True
                break
        if fault:
            dataAux.pop(-1)
            break
        dataAux.append(copyF0[m])  # guardar la fracuencia
        dataAux.append(m)  # guardar el indice del inicio de la nota
        while copyF0[m] != 0:
            m = m + 1
        dataAux.append(m)  # guardar el indice del final de la nota
        mainAux.append(dataAux) #agregar los datos al array general
        if len(mainAux) >=2: #eliminar los onset que son continuos
            if mainAux[-2][0] + 1 == mainAux[-1][0] or mainAux[-2][2] == mainAux[-1][2]: #solo dejar el ultimo de los peaks y eliminar el penultimo si
                mainAux.pop(-2)

#Graficar los peaks encontrados y obtener los parametros necesarios
onsetPeaksAux = []
notesFoundAux = []
notesMIDI = []
for n in range(len(mainAux)):
    plt.scatter(mainAux[n][0]*H/float(fs), ODF[mainAux[n][0], 1]*factorODF, color='black') #incrementar y * 20 para mejor visualización
    onsetPeaksAux.append(mainAux[n][0])
    notesFoundAux.append(mainAux[n][1])
    notesMIDI.append( math.log(mainAux[n][1] / 27.5 , 2**(1/12)) + 1 )
onsetPeaksAux.append(ODF.shape[0]) #agregar el largo total al array
onsetPeaks = np.array(onsetPeaksAux)
notesFound = np.array(notesFoundAux)
print(notesMIDI)
print(mainAux)

#graficar onset
graphOnset(ODF*factorODF) #*20
plt.xlabel("Tiempo")
plt.ylabel("Frecuencia")
plt.show()
#notesFound = [44, 47, 45, 40, 52] #Son las notas encontradas

#nRaws = shapeMatrix(listOnset=onsetPeaks,freqSample=44100,hopSize=H)
#print(nRaws)

for n in range(onsetPeaks.size - 1): #Recrear las notas y su respectiva duracion
    #print(n)
    noteDuration = (onsetPeaks[n+1] - onsetPeaks[n])*H/fs
    freqMatrixAux, magnitudDBAux, zeroPhaseAux = createTone(note=notesFound[n], duration=noteDuration, freqSample=fs, hopSize=150)
    #rawAux = freqMatrixAux.shape[0]
    #print(freqMatrix.shape)
    if n == 0:
        #yFinal = SMM.sineModelSynth(freqMatrixAux, magnitudDBAux, zeroPhaseAux, 600, 150, 44100)  # el valor de N debe cuatro veces H
        freqMatrixFinal = freqMatrixAux
        magnitudDBFinal = magnitudDBAux
        zeroPhaseFinal = zeroPhaseAux
        continue
    #yFinal =np.concatenate((yFinal,SMM.sineModelSynth(freqMatrixAux, magnitudDBAux, zeroPhaseAux, 600, 150, 44100)))
    freqMatrixFinal = np.concatenate((freqMatrixFinal, freqMatrixAux))
    magnitudDBFinal = np.concatenate((magnitudDBFinal, magnitudDBAux))
    zeroPhaseFinal = np.concatenate((zeroPhaseFinal, zeroPhaseAux))


yFinal = SMM.sineModelSynth(freqMatrixFinal, magnitudDBFinal, np.array([]), 600, 150, 44100)  # el valor de N debe cuatro veces H

#np.savetxt('test.txt', freqMatrixFinal)
waves.write('reCreate.wav', 44100, yFinal)

#final del archivo