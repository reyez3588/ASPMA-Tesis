import numpy as np
import matplotlib.pyplot as plt
import math
import sympy as sp

def transferFunction(fs, ts, xi):
    '''
    :param fs: Frecuencia de muestreo
    :param ts: tiempo de asentamiento
    :param xi: coeficiente de amortiguamiento relativo
    :return:
    '''
    Wn = 4/(xi*ts) #Frecuencia natural no amortiguada
    Wd = Wn * math.sqrt(1-xi**2) #Frecuencia natural amortiguada
    Ws = fs #Frecuencia de muestreo
    T= 1/fs
    #Calculo del polo dominante
    magZ = math.exp(-1*T*xi*Wn) #magnitud del polo dominante
    angZ = T*Wd #angulo del polo dominante
    angZDeg= angZ * 180 / math.pi # angulo del polo dominante en grados
    rectZ = complex(magZ*math.cos(angZ), magZ*math.sin(angZ)) # polo dominante en coordenadas polares

    # Es la funcion de tranferencia del retenedor de orden cero y una planta de segundo orden
    z = sp.Symbol('z')
    systemTranfer = ((((-1*math.exp(-2*T))*(0.5*T + 0.25)) + 0.25) + ((0.25*math.exp(-2*T)) + 0.5*T -0.25)*z) / ((z-1)*(z-math.exp(-2*T)))

    #Calculo de polos de la funcion anterior
    pole1 = 1
    pole2 = math.exp(-2*T)
    zero = (-1*(((-1*math.exp(-2*T))*(0.5*T + 0.25)) + 0.25)) / ((0.25*math.exp(-2*T)) + 0.5*T -0.25)

    #angulos con respecto al polo dominante
    anglepole1 = np.rad2deg(np.arctan2(rectZ.imag - 0, rectZ.real-pole1))
    anglepole2 = np.rad2deg(np.arctan2(rectZ.imag - 0, rectZ.real - pole2))
    angleZero = np.rad2deg(np.arctan2(rectZ.imag - 0, rectZ.real - zero))
    #print(anglepole1,anglepole2, angleZero) para ver los valores de angulos con respecto al polo dominante

    angleOut = -1*(180 - (anglepole1 + anglepole2 - angleZero)) #angulo de salida del polo dominante
    angleAlpha = anglepole2 - angleOut #angulo que debe tener el nuevo polo
    #print(angleAlpha) #el valor del angulo que debe tener el nuevo polo

    #Determinar la posicion en el eje z del nuevo polo
    newPole = -1
    while True:
        angleAct = anglepole1 = np.rad2deg(np.arctan2(rectZ.imag - 0, rectZ.real-newPole))
        if angleAct > angleAlpha:
            break
        newPole = newPole + 0.0001
    #print(newPole) #el valor del polo nuevo

    #Diseño del controlador
    #K = sp.Symbol('K')
    controlFunction = (z - pole2) / (z - newPole)  # le falta la k pero se calculara posteriormente
    transferFunctionOL = controlFunction * systemTranfer

    symbZ = rectZ
    # sustituye el polo dominante en la funcion de transferencia en lazo abierto
    subsZ = (sp.simplify(transferFunctionOL.subs(z, symbZ), complex = True)).as_real_imag()
    K = 1 / math.sqrt((subsZ[0]**2)+(subsZ[1]**2))
    controlFunction = K*(z - pole2) / (z - newPole) #Se actualiza la funcion del controlador ya con el valor K
    transferFunctionOL = (((((-1*math.exp(-2*T))*(0.5*T + 0.25)) + 0.25) + ((0.25*math.exp(-2*T)) + 0.5*T -0.25)*z) / (z-1))*(K / (z - newPole)) #quitando manualmente el polo eliminaado
    #controlFunction * systemTranfer #Se actualiza la funcion de tranferencia en lazo abierto ya con el valor K
    print('La funcion de tranferencia del controlador es: ', controlFunction)
    print('La funcion de transferencia en lazo abierto es:', sp.simplify(transferFunctionOL))

    #Funcion de tranferencia en lazo cerrado
    transferFunctionCL = sp.simplify(transferFunctionOL/(1+transferFunctionOL))
    print('La funcion de transferencia en lazo cerrado es:', sp.simplify(transferFunctionCL))

    '''
    posNew = np.array([[0], [0], [1]])
    print(posNew)
    a = traslateRotate(rectZ.real, rectZ.imag, np.pi, 'rad')
    print(a)
    b = np.matmul(a,posNew)
    print(b)
    '''
    drawCircle(1,rectZ,pole1,pole2,zero, newPole)

    return magZ, angZDeg, rectZ, systemTranfer

def drawCircle(radio,P, pole1, pole2, zero, newPole):
    '''
    :param radio:
    :param P: polo dominante a graficar
    :return:
    '''
    #Graficar circulo
    theta = np.linspace(0, 2 * np.pi, 100)
    radius = radio
    a = radius * np.cos(theta)
    b = radius * np.sin(theta)
    figure, axes = plt.subplots(1)
    axes.plot(a, b)
    axes.set_aspect(1)
    plt.title('Lugar geométrico de las raíces')

    #Graficar puntos
    plt.plot(P.real, P.imag, marker="X", color="red")
    plt.plot(pole1, 0, marker="X", color="blue")
    plt.plot(pole2, 0, marker="X", color="blue")
    plt.plot(zero, 0, marker="o", color="blue")
    plt.plot(newPole, 0, marker="X", color="green")

    #Mostrar grafica
    plt.grid()
    plt.show()
    return

def rot2(ang):
    matrixR = np.zeros([2,2])
    matrixR[0, 0] = math.cos(ang)
    matrixR[0, 1] = -1 * math.sin(ang)
    matrixR[1, 0] = math.sin(ang)
    matrixR[1, 1] = math.cos(ang)
    return matrixR

def transl2(x,y):
    matrixT = np.zeros((3, 3))
    matrixT[0:2, 0:2]=rot2(0)
    matrixT[0, 2] = x
    matrixT[1, 2] = y
    matrixT[2, 2] = 1
    return matrixT

def trot2(ang, typeAng):
    if typeAng == 'deg':
        matrixTR = np.zeros((3, 3))
        matrixTR[0: 2, 0: 2] = rot2(ang * np.pi / 180)
        matrixTR[0, 2] = 0
        matrixTR[1, 2] = 0
        matrixTR[2, 2] = 1
    if typeAng == 'rad':
        matrixTR = np.zeros((3, 3))
        matrixTR[0: 2, 0: 2] = rot2(ang * np.pi / 180)
        matrixTR[0, 2] = 0
        matrixTR[1, 2] = 0
        matrixTR[2, 2] = 1
    return matrixTR

def traslateRotate(x, y, ang, typeAng):
    a = transl2(x, y)
    b = trot2(ang, typeAng)
    matrizTR = np.matmul(a, b)
    return matrizTR

#z =  transferFunction(1/0.2, 2, 0.5)
z =  transferFunction(44100, 2, 0.5)

#M = traslateRotate(1,2,0.3,'deg')
#print(M)