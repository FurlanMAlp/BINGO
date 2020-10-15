# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 13:47:29 2020

@author: Acer
"""
import numpy as np
import numpy.fft as fft
import matplotlib as mpl
import matplotlib.pyplot as plt


#Definição do tamanho dos vetores e de variáveis auxiliares
n=4096
nm1 = n-1
nd2=n//2
M=int(np.log2(n))
l=nd2

#Definição dos sinais de teste arbitrários
x=np.linspace(-np.pi/2,np.pi/2,n)
ReX= np.zeros(n)
ImX= np.zeros(n)
A = [127, 63, 31]
F = [8, 64, 256]
Phi = [0, np.pi/8, np.pi/6]
ReX2=np.sum([a*np.sin(f*x+phi) for a,f,phi in zip (A,F,Phi)],axis=0)
for i in range(n):
    ReX[i]=np.sum([a*np.sin(f*x[i]+phi) for a,f,phi in zip (A,F,Phi)],axis=0)#+rand[i]

#C
T_fft=fft.fft(ReX2)        

#Algoritmo para reversão de bits
for i in range(1,nm1):
    if i<l:
        tr=ReX[l]
        ti=ImX[l]
        ReX[l]=ReX[i]
        ImX[l]=ImX[i]
        ReX[i]=tr
        ImX[i]=ti
        k=nd2
    else:
        k=nd2
    while k <= l:
        l=l-k
        k=k//2
    l=l+k
    
#Algoritmo para FFT
#Loop dos estágios da FFT
for l in range(1,M):
    le=int(2**l)
    le2=le//2
    ur=1 #Primeiros valores para a multiplicação pela senoide complexa
    ui=0
    sr=np.cos(np.pi/le2)
    si=np.sin(np.pi/le2)
    
    #Loop para cada espectro individual dentro do estágio
    for j in range(1,le2+1):
        jm1=j-1
        
        #Loop para o cálculo da borboleta
        for i in range(jm1,nm1,le):
            ip=i+le2
            #print(str(ip))
            tr = ReX[ip]*ur-ReX[ip]*ui
            ti= ImX[ip]*ui+ImX[ip]*ur
            ReX[ip]=ReX[ip]=ReX[i]-tr
            ImX[ip]=ImX[ip]-ti
            ReX[i]=ReX[i]+tr
            ImX[i]=ImX[i]+ti
        tr=ur
        ur=tr*sr-ui*si
        ui=tr*si+ui*sr

#Plotagem dos gráficos
plt.subplot(2,2,1)
plt.plot(x,ReX); #plt.xlim(-np.pi/32, np.pi/32)
plt.subplot(2,2,2)
plt.plot(x,ImX)
plt.subplot(2,2,(3,4))
plt.plot(fft.fftfreq(n),T_fft.real)


          
    


        

