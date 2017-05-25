import circuit as cir
from sympy import Heaviside, sin, cos, symbols
from scipy.signal import ss2tf
import matplotlib.pyplot as plt
import copy
import numpy as np
from math import isnan, isinf


def convolve_integral(x, y1, y2):
    delta = x[1] - x[0]
    y_ret = np.zeros(np.shape(x))
    for i in range(len(x) - 1):
        temp_y1 = y1[0:i+1]
        temp_y2 = y2[0:i+1]
        temp_y2 = np.flipud(temp_y2)
        y_int = temp_y1*temp_y2
        y_ret[i] = np.trapz(y_int, dx=delta)
    return y_ret


def delta_omega(x1,ACH,x2,Spectre):
    H707=ACH[1]*0.707
    i=1;
    while(ACH[i]>=H707):
        i=i+1
    dx=x1[i-1]
    dwx=[0,dx,dx+0.00001]
    dwy=[ACH[i-1],ACH[i-1],0]
    A01=Spectre[0]*0.1
    maxi=0
    for i in range(len(Spectre)):
        if (Spectre[i]>=A01):
            maxi=i
    dx=x1[maxi]
    dw1x=[0,x2[len(x2)-1]]
    dw1y=[Spectre[maxi],Spectre[maxi]]
    return dwx,dwy,dw1x,dw1y

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
        # Время сигнала
        self.Tau = Tau
        # Амплитуда сигнала
        self.A = A
        # Тип сигнала
        self.sigtype = sigtype

    def H_s(self,num, type):
        ABCD = cir.Circuit.StateSpace(self, num, type)
        HStp=ss2tf(ABCD[0],ABCD[1],ABCD[2],ABCD[3])
        H1S_zn = np.zeros(len(HStp[1]) + 1)
        for i in range(len(HStp[1])):
            H1S_zn[i] = HStp[1][i]
        H1S = HStp[0][0], H1S_zn
        HS = HStp[0][0], HStp[1]
        # Возвращает H(s), H1(s)
        return HS, H1S,ABCD

    # N - Количество точек для преобразования фурье
    def Fourier(self, stype, h_s, A, Tau, Period,  K = 7, N = 200):
        # Удаляем плохие значения хардкодом(чтобы не считать пределы)
        if (Period / Tau) % 2 == 0:
            Period = Period + 0.0001 * abs(Period - Tau)
        # Задаём шаг
        T = self.Tau/N
        xf = np.linspace( 0.001 , 1.0 / (2.0 * T), N)
        # аргумент x для значений f1(t) и f2(t)
        x = np.linspace(0, Period, N)
        # Значения k*w0 для разложения входного и выходного сигнала по частотам
        omega = np.zeros(K+1)
        for i in range(K):
            omega[i+1] = 2*np.pi/Period * (i+1)
        omega[0] =0.00001
        # Фазовый и амплитудный спектры
        yp = np.zeros(np.shape(xf))
        ya = np.zeros(np.shape(xf))
        # Фазы и амплитуды для разложения
        a1 = np.zeros(np.shape(omega))
        fi1 = np.zeros(np.shape(omega))
        # Считаем амплитудный и фазовый спектры входного сигнала
        if stype == 1:
            temp = A*Tau*np.pi*(np.exp(-Tau*xf*1.j) + 1)/((Tau*xf*1.j)**2 + np.pi**2)
            ya = np.abs(temp)
            yp = np.angle(temp)
            temp = A * Tau * np.pi*(np.exp(-Tau* 1.j * omega) + 1) / ((Tau * 1.j * omega) ** 2 + np.pi ** 2)
            a1 = np.abs(2*temp/Period)
            fi1 = np.angle(temp)
        if stype == 2:
            temp = A*np.exp(-Tau*xf*1.j)*(np.exp(Tau/2*xf*1.j) - 1)**2/(Tau/2*(xf*1.j)**2)
            ya = np.abs(temp)
            yp = np.angle(temp)
            temp = A*np.exp(-Tau*omega*1.j)*(np.exp(Tau/2*omega*1.j) - 1)**2/(Tau/2*(omega*1.j)**2)
            a1 = np.abs(2*temp/Period)
            fi1 = np.angle(temp)
        if stype == 3:
            temp = A*(Tau**2)*(np.exp(-Tau*xf*1.j) + 1)*xf*1.j/((Tau*xf*1.j)**2 + np.pi**2)
            ya = np.abs(temp)
            yp = np.angle(temp)
            temp = A*(Tau**2)*(np.exp(-Tau*omega*1.j) + 1)*omega*1.j/((Tau*omega*1.j)**2 + np.pi**2)
            a1 = np.abs(2*temp/Period)
            fi1 = np.angle(temp)
        if stype == 4:
            temp = A*np.exp(-Tau*xf*1.j)*(np.exp((Tau/2)*xf*1.j) - 1)**2/(xf*1.j)
            ya = np.abs(temp)
            yp = np.angle(temp)
            temp = A*np.exp(-Tau*omega*1.j)*(np.exp((Tau/2)*omega*1.j) - 1)**2/(omega*1.j)
            a1 = np.abs(2*temp/Period)
            fi1 = np.angle(temp)
        # Считаем H(jw)
        fun_res_ch = np.zeros(np.shape(omega), dtype=complex)
        fun_res_zn = np.zeros(np.shape(omega), dtype=complex)
        for i in range(len(h_s[0])):
            fun_res_ch = fun_res_ch + (h_s[0][i]) * (omega*1.j) ** (len(h_s[0]) - i)
        for i in range(len(h_s[1])):
            fun_res_zn = fun_res_zn + (h_s[1][i]) * (omega*1.j) ** (len(h_s[1]) - i)
        y = np.ndarray(shape = np.shape(fun_res_ch), dtype=complex, buffer = fun_res_ch/fun_res_zn)
        # Расчитываем A2 и Ф2
        a2 = np.abs(y)
        fi2 = np.angle(y)
        a2 = a1*a2
        fi2 = fi1+fi2
        yf1 = np.zeros(np.shape(x))
        yf2 = np.zeros(np.shape(x))
        yf1 = yf1+a1[0]/2
        yf2 = yf2+a2[0]/2
        # Расчитываем f1(t) и f2(t) через частоты и амплитуды
        for i in range(K):
            yf1 = yf1 + a1[i+1]*np.cos(omega[i+1]*x + fi1[i+1])
            yf2 = yf2 + a2[i + 1] * np.cos(omega[i + 1] * x + fi2[i + 1])
        # Чиним фазовый спектр
        f = 0
        for i in range(len(yp) - 1):
            if (f == 1):
                if (yp[i + 1] < 0):
                    yp[i + 1] = (abs(yp[i + 1]) / yp[i + 1]) * ((abs(yp[i + 1]) + 2 * np.pi) % (2 * np.pi))
                else:
                    yp[i + 1] = yp[i + 1] - 2 * np.pi
            else:
                if (yp[i + 1] > yp[i]):
                    f = 1
                    if (yp[i + 1] < 0):
                        yp[i + 1] = (abs(yp[i + 1]) / yp[i + 1]) * ((abs(yp[i + 1]) + 2 * np.pi) % (2 * np.pi))
                    else:
                        yp[i + 1] = yp[i + 1] - 2 * np.pi
        return xf, ya, yp, x, yf1, yf2, omega, a1, fi1, a2, fi2      

            

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
        x = np.linspace(deltay, 2*N * T, 2*N)
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
        return x, yh, yh1, h, hs[2]

    def frequency_analysis(self, h_s, N = 200):
        # Задаём шаг
        T = self.Tau/N
        xf = np.linspace( 0 , 1.0 / (2.0 * T), N, dtype=complex)
        fun_res_ch = np.zeros(np.shape(xf), dtype=complex)
        fun_res_zn = np.zeros(np.shape(xf), dtype=complex)
        # Заполняем массив значений сигнала
        for i in range(len(h_s[0])):
            fun_res_ch = fun_res_ch + (h_s[0][i]) * (xf*1j) ** (len(h_s[0]) - i)
        for i in range(len(h_s[1])):
            fun_res_zn = fun_res_zn + (h_s[1][i]) * (xf*1j) ** (len(h_s[1]) - i)
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
Stype=1
Tau=20
Period=40
t = symbols("t", positive=True)
node_array = [cir.Node(0), cir.Node(1), cir.Node(2), cir.Node(3)]
node_elem = [cir.Element(0, "U", 1, None, 1, 0, 1), cir.Element(1, "R", 1, None, None, 1, 3),
             cir.Element(2, "C", 1, None, None, 3, 0), cir.Element(3, "R", 1, None, None, 3, 2),
             cir.Element(4, "C", 0.25, None, None, 2, 0), cir.Element(5, "R", 1, None, None, 2, 0)]
circ = SignalCircuit(node_array, node_elem, Stype, A, Tau)
hxy = circ.hth1tGraph(5, 'U', 100, 0.00005)
sigxy = circ.SignalGraph(hxy[0])
f2x = sigxy[0]
f2y = convolve_integral(f2x, sigxy[1], hxy[1])
frq_a = circ.frequency_analysis(hxy[3])
f1f2 = circ.Fourier(Stype, hxy[3], A,  Tau, Period)
delta=delta_omega(frq_a[0], frq_a[1],f1f2[0], f1f2[1])
# f1f2=circ.FourierPeriod(Stype,hxy[3],A,Tau,Period)
print('H(S): num ',hxy[3][0],'denum',hxy[3][1])
print('A:',hxy[4][0])
print('B:',hxy[4][1])
print('C:',hxy[4][2])
print('D:',hxy[4][3])

# Сигнал
plt.figure(1)
plt.xlabel('t')
plt.ylabel('f1(t)')
plt.title('Signal f1(t)')
plt.grid()
plt.plot(sigxy[0], sigxy[1])
plt.plot(f1f2[3],f1f2[4])

# Амплитуда входного сигнала
plt.figure(2)
plt.xlabel('Omega')
plt.ylabel('|A|')
plt.title('Ampletude spectre IN')
plt.grid()
plt.plot(f1f2[0], f1f2[1])
plt.plot(delta[2],delta[3],'r--')
plt.show()


print('dw1=',delta[2][1])


# Фаза входного сигнала
plt.figure(3)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase spectre IN')
plt.grid()
plt.plot(f1f2[0], f1f2[2])

# Выходной сигнал
plt.figure(4)
plt.xlabel('t')
plt.ylabel('f2(t)')
plt.title('Reaction f2(t)')
plt.grid()
plt.plot(f2x, f2y)
plt.plot(f1f2[3], f1f2[5])

# Амплитуда выходного сигнала
plt.figure(5)
plt.xlabel('Omega')
plt.ylabel('|A|')
plt.title('Ampletude spectre OUT')
plt.grid()
plt.plot(f1f2[0], f1f2[1]*frq_a[1])

# Фаза выходного сигнала
plt.figure(6)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase spectre OUT')
plt.grid()
plt.plot(f1f2[0], f1f2[2]+frq_a[3])

# h(t)
plt.figure(7)
plt.xlabel('t')
plt.ylabel('h(t)')
plt.title('h(t)')
plt.grid()
plt.plot(hxy[0], hxy[1])

# h1(t)
plt.figure(8)
plt.xlabel('t')
plt.ylabel('h1(t)')
plt.title('h1(t)')
plt.grid()
plt.plot(hxy[0], hxy[2])

# АЧХ
plt.figure(9)
plt.xlabel('Omega')
plt.ylabel('|H(jw)|')
plt.title('АЧХ')
plt.grid()
plt.plot(frq_a[0], frq_a[1])
plt.plot(delta[0],delta[1],'r--')
plt.show()

print('dw1=',delta[0][1])

# ФЧХ
plt.figure(10)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('ФЧХ')
plt.grid()
plt.plot(frq_a[2], frq_a[3])


# Амплитуда Фурье  IN
plt.figure(11)
plt.xlabel('Omega')
plt.ylabel('A')
plt.title('Ampletude Fourier IN')
plt.grid()
plt.stem(f1f2[6], f1f2[7])


# Фаза Фурье IN
plt.figure(12)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase Fourier IN')
plt.grid()
plt.stem(f1f2[6], f1f2[8])


# Амплитуда Фурье Out
plt.figure(13)
plt.xlabel('Omega')
plt.ylabel('A')
plt.title('Ampletude Fourier OUT')
plt.grid()
plt.stem(f1f2[6], f1f2[9])


# Фаза Фурье OUT
plt.figure(14)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase Fourier OUT')
plt.grid()
plt.stem(f1f2[6], f1f2[10])
plt.show()

pass
