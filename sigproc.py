import circuit as cir
from sympy.abc import  s, t
from sympy import Heaviside, sin, cos, exp
from numpy import pi
from sympy.integrals import laplace_transform, inverse_laplace_transform
from sympy.parsing.sympy_parser import parse_expr
from scipy.fftpack import  fft, fftshift
from scipy.signal import ss2tf
import matplotlib.pyplot as plt
import copy
import numpy as np


class SignalCircuit(cir.Circuit):

    def __init__(self, node_array, el_array, sigtype, A, Tau, signal = 0):
        cir.Circuit.__init__(self, node_array, el_array)
        if sigtype == 1:
            self.signal = A*(sin(t*pi/Tau) - sin(t*pi/Tau)*Heaviside(t-Tau))
        if sigtype == 2:
            self.signal = 2*A/Tau*(t-t*Heaviside(t - Tau/2))+(-2*A/Tau*t +2*A)*(Heaviside(t- Tau/2)-Heaviside(t - Tau))
        if sigtype == 3:
            self.signal = A*(cos(t*pi/Tau) - cos(t*pi/Tau)*Heaviside(t-Tau))
        if sigtype == 4:
            self.signal = A*(1 - 2*Heaviside(t-Tau/2) + Heaviside(t-Tau))
        # Сигнал, заданный пользователем
        if sigtype == 5:
            self.signal = signal
        # Время сигнала
        self.Tau = Tau

    def H_t(self,num, type):
        ABCD = cir.Circuit.StateSpace(self, num, type)
        H_Sold=ss2tf(ABCD[0],ABCD[1],ABCD[2],ABCD[3])
        H_S = []
        # округляем значения до второй значащей цифры
        H_S.append(np.around(H_Sold[0], decimals = 4))
        H_S.append(np.around(H_Sold[1], decimals = 4))
        HS0=''
        HS1=''
        f=0
        j=H_S[0][0].size-1
        for i in H_S[0][0]:
            num=str(i)
            if i:
                if(j>1):
                    num=num+'*s**'+str(j)
                else:
                    if (j>0):
                        num=num+'*s'
                if f:
                    if (i<0):
                        HS0=HS0+'-'
                    else:
                        HS0=HS0+'+'
                else:
                    f=1
                HS0=HS0+num
            j=j-1
        j=H_S[1].size-1
        f=0
        for i in H_S[1]:
            num=str(i)
            if i:
                if(j>1):
                    num=num+'*s**'+str(j)
                else:
                    if (j>0):
                        num=num+'*s'
                if f:
                    if (i<0):
                        HS1=HS1+'-'
                    else:
                        HS1=HS1+'+'
                else:
                    f=1
                HS1=HS1+num
            j=j-1
        HS='('+HS0+')'+'/('+HS1+')'
        HS = parse_expr(HS, evaluate=False)
        HS1 = HS/s
        # Возвращает h(t), h1(t)
        return (inverse_laplace_transform(HS,s,t)), (inverse_laplace_transform(HS1,s,t))

    # N - Количество точек
    def Fourier(self, N = 200):
        # Задаём шаг
        T = self.Tau/N
        x = np.linspace(0.0, N * T, N)
        y = []
        f2 = copy.deepcopy(self.signal)
        # Заполняем массив значений сигнала
        for i in x:
            if int(i) == i:
                f2 = f2.subs(Heaviside(t - int(i)), 1)
            else:
                f2 = f2.subs(Heaviside(t - i), 1)
            y.append(f2.subs(t, i))
        # Заполняем значения амплитудного спектра
        yf = fftshift(fft(y))
        # Заполняем частоты
        xf = np.linspace(-1.0 / (2.0 * T), 1.0 / (2.0 * T), N)
        yf2 = copy.deepcopy(2.0 / N * yf)
        # Фильтруем помехи
        threshold = np.max(np.abs(yf2)) / 10000
        for i in range(len(yf2)):
            if abs(yf2[i]) < threshold:
                yf2[i] = 0
        xp = []
        yp = []
        # Заполняем новые массивы без помех
        for i in range(len(yf2)):
            if yf2[i] != 0:
                xp.append(xf[i])
                yp.append(yf2[i])
        yp = np.angle(yp)
        yf = np.abs(yf)
        # Первые 2 элемента х и у амплитудного спектра, вторые 2 фазового
        return xf, yf, xp, yp

    # Создаёт точки для графика сигнала
    def SignalGraph(self,x):
        y = []
        f2 = copy.deepcopy(self.signal)
        # Заполняем массив значений сигнала
        for i in x:
            if int(i) == i:
                f2 = f2.subs(Heaviside(t - int(i)), 1)
            else:
                f2 = f2.subs(Heaviside(t - i), 1)
            y.append(f2.subs(t, i))
        return x, y

    # Создаёт точки для графиков h(t), h1(t) f2(t) для элемента num и типа реакции type
    def hth1tGraph(self, num, type, N = 200, deltay = 0.0001):
        T = self.Tau / N
        x = np.linspace(0.0, N * T, N)
        hs = self.H_t(num,type)
        h = hs[0]
        h1 = hs[1]
        yh = []
        for i in x:
            if int(i) == i:
                h = h.subs(Heaviside(t - int(i)), 1)
            else:
                h = h.subs(Heaviside(t - i), 1)
            yh.append(h.subs(t, i))
        yh1 = []
        for i in x:
            if int(i) == i:
                h1 = h1.subs(Heaviside(t - int(i)), 1)
            else:
                h1 = h1.subs(Heaviside(t - i), 1)
            yh1.append(h1.subs(t, i))
        while abs(yh[len(yh) - 1] - yh[len(yh) - 2]) > deltay:
            x = np.append(x, x[len(x) - 1] + T)
            i = x[len(x) - 1]
            if int(i) == i:
                h = h.subs(Heaviside(t - int(i)), 1)
                h1 = h1.subs(Heaviside(t - int(i)), 1)
            else:
                h = h.subs(Heaviside(t - i), 1)
                h1 = h1.subs(Heaviside(t - i), 1)
            yh.append(h.subs(t, i))
            yh1.append(h1.subs(t, i))
        return x, yh, yh1




node_array = [cir.Node(0), cir.Node(1), cir.Node(2), cir.Node(3)]
node_elem = [cir.Element(0, "I", 1, 1, None, 0, 1), cir.Element(1, "R", 2, None, None, 1, 0),
             cir.Element(2, "L", 2, None, None, 1, 2), cir.Element(3, "R", 1, None, None, 2, 3),
             cir.Element(4, "C", 4, None, None, 3, 0), cir.Element(5, "R", 0.5, None, None, 3, 0)]
circ = SignalCircuit(node_array, node_elem, 3, 1, 2)
hxy = circ.hth1tGraph(5,'I', 100, 0.00001)
fftxy = circ.Fourier()
sigxy = circ.SignalGraph(hxy[0])
f2x = sigxy[0]
f2xy = np.convolve(sigxy[1], hxy[1])
f2y = f2xy[0:len(sigxy[0])]

# Сигнал
plt.figure(1)
plt.xlabel('t')
plt.ylabel('f1(t)')
plt.title('Signal f1(t)')
plt.grid()
plt.plot(sigxy[0],sigxy[1])

# Амплитуда
plt.figure(2)
plt.xlabel('Omega')
plt.ylabel('|A|')
plt.title('Ampletude spectre')
plt.grid()
plt.plot(fftxy[0],fftxy[1])

# Фаза
plt.figure(3)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase spectre')
plt.grid()
plt.plot(fftxy[2],fftxy[3])

# Выходной сигнал
plt.figure(4)
plt.xlabel('t')
plt.ylabel('f2(t)')
plt.title('Reaction f2(t)')
plt.grid()
plt.plot(f2x,f2y)

# h(t)
plt.figure(5)
plt.xlabel('t')
plt.ylabel('h(t)')
plt.title('h(t)')
plt.grid()
plt.plot(hxy[0],hxy[1])

# h1(t)
plt.figure(6)
plt.xlabel('t')
plt.ylabel('h1(t)')
plt.title('h1(t)')
plt.grid()
plt.plot(hxy[0],hxy[2])

plt.show()
pass