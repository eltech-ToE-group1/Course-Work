import circuit as cir
from sympy import Heaviside, sin, cos, symbols
from scipy.fftpack import fft, fftshift
from scipy.signal import ss2tf
import matplotlib.pyplot as plt
import copy
import numpy as np

# Расчёт значений обратного преобразования лапласа функции f_s, заданной матрично в точках t_z
def talbot_inverse(f_s, t_z, M = 64):
    k = np.arange(M)
    delta = np.zeros(M, dtype=complex)
    for i in k:
        if i != 0:
            delta[i] = 2*np.pi/5 * i * (np.tan(np.pi/M*i)**(-1)+1.j)
    delta[0] = 2*M/5
    gamma = np.zeros(M, dtype=complex)
    for i in k:
        if i != 0:
            gamma[i] = (1 + 1.j*np.pi/M*i*(1+np.tan(np.pi/M*i)**(-2))-1.j*np.tan(np.pi/M*i)**(-1))*np.exp(delta[i])
    gamma[0] = 0.5*np.exp(delta[0])
    delta_mesh, t_mesh = np.meshgrid(delta, t_z)
    gamma_mesh = np.meshgrid(gamma,t_z)[0]
    points = delta_mesh/t_mesh
    fun_res_ch = np.zeros(np.shape(points), dtype=complex)
    fun_res_zn = np.zeros(np.shape(points), dtype=complex)
    for i in range(len(f_s[0])):
        fun_res_ch = fun_res_ch + (f_s[0][i])*points**(len(f_s[0]) - i)
    for i in range(len(f_s[1])):
        fun_res_zn = fun_res_zn + (f_s[1][i])*points**(len(f_s[1]) - i)
    fun_res = fun_res_ch/fun_res_zn
    sum_ar = np.real(gamma_mesh*fun_res)
    sum_ar = np.sum(sum_ar, axis = 1)
    ilt = 0.4/t_z * sum_ar
    return ilt

class SignalCircuit(cir.Circuit):

    def __init__(self, node_array, el_array, sigtype, A, Tau, signal = 0):
        cir.Circuit.__init__(self, node_array, el_array)
        if sigtype == 1:
            self.signal = A*(sin(t*np.pi/Tau) - sin(t*np.pi/Tau)*Heaviside(t-Tau))
        if sigtype == 2:
            self.signal = 2*A/Tau*(t-t*Heaviside(t - Tau/2))+(-2*A/Tau*t +2*A)*(Heaviside(t- Tau/2)-Heaviside(t - Tau))
        if sigtype == 3:
            self.signal = A*(cos(t*np.pi/Tau) - cos(t*np.pi/Tau)*Heaviside(t-Tau))
        if sigtype == 4:
            self.signal = A*(1 - 2*Heaviside(t-Tau/2) + Heaviside(t-Tau))
        # Сигнал, заданный пользователем
        if sigtype == 5:
            self.signal = signal
        # Время сигнала
        self.Tau = Tau

    def H_s(self,num, type):
        ABCD = cir.Circuit.StateSpace(self, num, type)
        HStp=ss2tf(ABCD[0],ABCD[1],ABCD[2],ABCD[3])
        H1S_zn = np.zeros(len(HStp[1]) + 1)
        for i in range(len(HStp[1])):
            H1S_zn[i] = HStp[1][i]
        H1S = HStp[0][0], H1S_zn
        HS = HStp[0][0], HStp[1]
        # Возвращает H(s), H1(s)
        return HS, H1S

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
        x = np.linspace(deltay, N * T, N)
        hs = self.H_s(num,type)
        h = hs[0]
        h1 = hs[1]
        yh = talbot_inverse(h,x)
        yh1 = talbot_inverse(h1,x)
        while abs(yh[len(yh) - 1] - yh[len(yh) - 2]) > deltay:
            temp = x[len(x) - 1] + T
            x = np.append(x,temp)
            yh = np.append(yh, talbot_inverse(h,temp))
            yh1 = np.append(yh1, talbot_inverse(h1, temp))
        return x, yh, yh1


t = symbols("t", positive=True)
node_array = [cir.Node(0), cir.Node(1), cir.Node(2), cir.Node(3)]
node_elem = [cir.Element(0, "U", 1, None, 1, 0, 1), cir.Element(1, "C", 4, None, None, 1, 2),
             cir.Element(2, "R", 1, None, None, 2, 1), cir.Element(3, "C", 1, None, None, 2, 3),
             cir.Element(4, "R", 0.5, None, None, 3, 0), cir.Element(5, "R", 0.5, None, None, 3, 0)]
circ = SignalCircuit(node_array, node_elem, 3, 1, 5)
hxy = circ.hth1tGraph(5,'U', 100)
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