import circuit as cir
from sympy import Heaviside, sin, cos, symbols
from scipy.fftpack import fft, fftshift, fftfreq
from scipy.signal import ss2tf
import matplotlib.pyplot as plt
import copy
import numpy as np
from math import exp,isnan

# Расчёт значений обратного преобразования лапласа функции f_s, заданной матрично в точках t_z
def talbot_inverse(f_s, t_z, M = 64):
    k = np.arange(M)
    # Задаём вектор значений функции дельта
    delta = np.zeros(M, dtype=complex)
    for i in k:
        if i != 0:
            delta[i] = 2*np.pi/5 * i * (np.tan(np.pi/M*i)**(-1)+1.j)
    delta[0] = 2*M/5
    # Задаём вектор значений функции гамма
    gamma = np.zeros(M, dtype=complex)
    for i in k:
        if i != 0:
            gamma[i] = (1 + 1.j*np.pi/M*i*(1+np.tan(np.pi/M*i)**(-2))-1.j*np.tan(np.pi/M*i)**(-1))*np.exp(delta[i])
    gamma[0] = 0.5*np.exp(delta[0])
    # Создаём сетки, чтобы избежать использования циклов
    delta_mesh, t_mesh = np.meshgrid(delta, t_z)
    gamma_mesh = np.meshgrid(gamma,t_z)[0]
    points = delta_mesh/t_mesh
    # Ищем значене f(s) нужных точках
    fun_res_ch = np.zeros(np.shape(points), dtype=complex)
    fun_res_zn = np.zeros(np.shape(points), dtype=complex)
    # Обходим числитель
    for i in range(len(f_s[0])):
        fun_res_ch = fun_res_ch + (f_s[0][i])*points**(len(f_s[0]) - i)
    # Обходим знаменатель
    for i in range(len(f_s[1])):
        fun_res_zn = fun_res_zn + (f_s[1][i])*points**(len(f_s[1]) - i)
    fun_res = fun_res_ch/fun_res_zn
    # Выделяем вещественную часть поэлементного произведения матриц
    sum_ar = np.real(gamma_mesh*fun_res)
    # Суммируем столбцы
    sum_ar = np.sum(sum_ar, axis = 1)
    # Получаем значение f(t) в заданных точках
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
        # Заполняем нули
        y = np.append(y, np.zeros(N))
        N = N*2
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
    
    def FourierPeriod(self, Stype, H_S, A, Tau,Period = 0, N = 7):
        HS0=0
        HS1=0
        j=H_S[0].size-1;
        for i in H_S[0]:
            num=i
            if i:
                if(j>1):
                    num=num*(t**j)
                else:
                    if (j>0):
                        num=num*(t)
                HS0=HS0+num
            j=j-1
        j=H_S[1].size-1
        for i in H_S[1]:
            num=i
            if i:
                if(j>1):
                    num=num*(t**j)
                else:
                    if (j>0):
                        num=num*(t)
                HS1=HS1+num
            j=j-1
        HS=HS0/(HS1)
        if (Period == 0):
            Period=self.Tau*2
        OMEGA = 2*np.pi/Period
        x = np.linspace(0.0,Period,200)
        xf = np.linspace(0.0,OMEGA*N,N+1)
        #1-Синус
        #2-треугольный импульс
        #3-Косинус
        #4-прямоугольный импульс
        if (Stype==1):
            A_jw=(2*A*np.pi/(Tau*(-t*t+np.pi*np.pi/(Tau*Tau))))*sin((Tau/2)*t)
            F_jw=-(Tau/2)*t
        if (Stype==2):
            A_jw=(-8*A/(Tau*t*t))*sin((Tau/4)*t)*sin((Tau/4)*t)
            F_jw=np.pi-(Tau/2)*t
        if (Stype==3):
            A_jw=(-2*A*t/(-t*t+np.pi*np.pi/(Tau*Tau)))*sin((Tau/2)*t)
            F_jw=np.pi/2+(Tau/2)*t
        if (Stype==4):
            A_jw=(4*A/t)*sin((Tau/4)*t)*sin((Tau/4)*t)
            F_jw=np.pi/2-(Tau/2)*t
        xs=np.linspace(0.0,3,1201)
        As=[]
        Fs=[]
        for i in range(len(xs)):
            As.append(A_jw.subs(t,xs[i]))
            if (isnan(abs(complex(As[i]))) or isinf(abs(complex(As[i])))):
                As[i]=0 
            if (abs(As[i])<0.01):
                As[i]=0
            #print(xs[i],As[i])
            if (As[i]==0):
                j=0;
            Fs.append(F_jw.subs(t,xs[j]))
            j=j+1
        a1 = []
        f1=[]
        a2=[]
        f2=[]
        for i in range(len(xf)):
            a1.append(A_jw.subs(t,xf[i])/(Period/2))
            if (isnan(abs(complex(a1[i]))) or isinf(abs(complex(a1[i])))):
                a1[i]=0 
            if (abs(a1[i])<0.00001):
                a1[i]=0
            if (a1[i]==0):
                j=0;
            f1.append(F_jw.subs(t,xf[j]))
            a2.append(a1[i]*np.abs(HS.subs(t,xf[i]*1.j)))
            f2.append(f1[i]+np.angle(complex(HS.subs(t,xf[i]*1.j))))
            j=j+1
            # Заполняем значения амплитудного спектра
        #print(a1)
        F1=a1[0]/2
        F2=a2[0]/2
        for i in range(len(f1)-1):
            F1=F1+a1[i+1]*cos(OMEGA*(i+1)*t+f1[i+1])
            F2=F2+a2[i+1]*cos(OMEGA*(i+1)*t+f2[i+1])
        Y1=[]
        Y2=[]
        for i in x:
            Y1.append(F1.subs(t,i))
            Y2.append(F2.subs(t,i))            
        return x,Y1,Y2,xs,As,Fs
        

            

    # Создаёт точки для графика сигнала
    def SignalGraph(self, x):
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
        hs = self.H_s(num, type)
        h = hs[0]
        h1 = hs[1]
        yh = talbot_inverse(h, x)
        yh1 = talbot_inverse(h1, x)
        while abs(yh[len(yh) - 1] - yh[len(yh) - 2]) > deltay:
            temp = x[len(x) - 1] + T
            x = np.append(x,temp)
            yh = np.append(yh, talbot_inverse(h, temp))
            yh1 = np.append(yh1, talbot_inverse(h1, temp))
        return x, yh, yh1, h

    def frequency_analysis(self, h_s, N = 100):
        # Задаём шаг
        T = self.Tau/N
        xf = np.linspace( 0 , 1.0 / (2.0 * T), N, dtype=complex)
        fun_res_ch = np.zeros(np.shape(xf), dtype=complex)
        fun_res_zn = np.zeros(np.shape(xf), dtype=complex)
        # Заполняем массив значений сигнала
        for i in range(len(h_s[0])):
            fun_res_ch = fun_res_ch + (h_s[0][i]) * (xf*1.j) ** (len(h_s[0]) - i)
        for i in range(len(h_s[1])):
            fun_res_zn = fun_res_zn + (h_s[1][i]) * (xf*1.j) ** (len(h_s[1]) - i)
        # Заполняем значения амплитудного спектра
        y = np.ndarray(shape = np.shape(fun_res_ch), dtype=complex, buffer = fun_res_ch/fun_res_zn)
        # Заполняем частоты
        yf = np.abs(y)
        yf2 = copy.deepcopy(2.0 / N * y)
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

A=10
Stype=4
Tau=20
Period=Tau*2
t = symbols("t", positive=True)
node_array = [cir.Node(0), cir.Node(1), cir.Node(2), cir.Node(3)]
node_elem = [cir.Element(0, "I", 1, 1, None, 0, 1), cir.Element(1, "R", 2, None, None, 1, 0),
             cir.Element(2, "L", 2, None, None, 1, 2), cir.Element(3, "R", 1, None, None, 2, 3),
             cir.Element(4, "C", 4, None, None, 3, 0), cir.Element(5, "R", 0.5, None, None, 3, 0)]
circ = SignalCircuit(node_array, node_elem, Stype, A, Tau)
hxy = circ.hth1tGraph(5, 'I', 100, 0.00005)
fftxy = circ.Fourier()
sigxy = circ.SignalGraph(hxy[0])
f2x = sigxy[0]
f2xy = np.convolve(sigxy[1], hxy[1])
f2y = f2xy[0:len(sigxy[0])]
frq_a = circ.frequency_analysis(hxy[3])
F1F2=circ.FourierPeriod(Stype,hxy[3],A,Tau,Period)

# Сигнал
plt.figure(1)
plt.xlabel('t')
plt.ylabel('f1(t)')
plt.title('Signal f1(t)')
plt.grid()
plt.plot(sigxy[0], sigxy[1])
plt.plot(F1F2[0], F1F2[1])

# Амплитуда
plt.figure(2)
plt.xlabel('Omega')
plt.ylabel('|A|')
plt.title('Ampletude spectre')
plt.grid()
plt.plot(F1F2[3], F1F2[4])

# Фаза
plt.figure(3)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase spectre')
plt.grid()
plt.plot(F1F2[3], F1F2[5])

# Выходной сигнал
plt.figure(4)
plt.xlabel('t')
plt.ylabel('f2(t)')
plt.title('Reaction f2(t)')
plt.grid()
plt.plot(f2x, f2y)
plt.plot(F1F2[0], F1F2[2])

# h(t)
plt.figure(5)
plt.xlabel('t')
plt.ylabel('h(t)')
plt.title('h(t)')
plt.grid()
plt.plot(hxy[0], hxy[1])

# h1(t)
plt.figure(6)
plt.xlabel('t')
plt.ylabel('h1(t)')
plt.title('h1(t)')
plt.grid()
plt.plot(hxy[0], hxy[2])

# АЧХ
plt.figure(7)
plt.xlabel('Omega')
plt.ylabel('|H(jw)|')
plt.title('АЧХ')
plt.grid()
plt.plot(frq_a[0], frq_a[1])

# ФЧХ
plt.figure(8)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('ФЧХ')
plt.grid()
plt.plot(frq_a[2], frq_a[3])
plt.show()

pass
