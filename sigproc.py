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
        return HS, H1S

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

    """def FourierPeriod(self, Stype, H_S, A, Tau,Period = 0, N = 7):
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
        As=[]
        Fs=[]
        xs=np.linspace(0.0,3,1201)
        a1 = []
        f1=[]
        a2=[]
        f2=[]
        if (Stype==1):
            A_jw=(2*A*np.pi/(Tau*(-t*t+np.pi*np.pi/(Tau*Tau))))*sin((Tau/2)*t)
            F_jw=-(Tau/2)*t
            As.append(0)
            Fs.append(np.pi/2)
            a1.append(0)
            f1.append(np.pi/2)
            a2.append(0)
            f2.append(f1[0]+np.angle(complex(HS.subs(t,xf[0]*1.j))))
        if (Stype==2):
            A_jw=(8*A/(Tau*t*t))*sin((Tau/4)*t)*sin((Tau/4)*t)
            F_jw=-(Tau/2)*t
            As.append(A*Tau/2)
            Fs.append(0)
            a1.append(A*Tau/Period)
            f1.append(0)
            a2.append(a1[0]*np.abs(HS.subs(t,xf[0]*1.j)))
            f2.append(f1[0]+np.angle(complex(HS.subs(t,xf[0]*1.j))))
        if (Stype==3):
            A_jw=(-2*A*t/(-t*t+np.pi*np.pi/(Tau*Tau)))*sin((Tau/2)*t)
            F_jw=np.pi/2-(Tau/2)*t
            As.append(0)
            Fs.append(0)
            a1.append(0)
            f1.append(0)
            a2.append(0)
            f2.append(f1[0]+np.angle(complex(HS.subs(t,xf[0]*1.j))))
        if (Stype==4):
            A_jw=(4*A/t)*sin((Tau/4)*t)*sin((Tau/4)*t)
            F_jw=np.pi/2-(Tau/2)*t
            As.append(0)
            Fs.append(np.pi/2)
            a1.append(0)
            f1.append(np.pi/2)
            a2.append(0)
            f2.append(f1[0]+np.angle(complex(HS.subs(t,xf[0]*1.j))))
        temp=A_jw.subs(t,xs[1])
        sflag=0
        fflag=0
        if  (Stype==1 or Stype==3):
            for i in range(len(xs)):
                if (xs[i]==np.pi/Tau):
                    sflag=1
                    scount=i
        if  (Stype==1 or Stype==3):
            for i in range(len(xf)):
                if (xf[i]==np.pi/Tau):
                    fflag=1
                    fcount=i
        for k in range(len(xs)-1):
            i=k+1
            As.append(temp)
            if (i<len(xs)-1):
                temp=A_jw.subs(t,xs[i+1])
            else:
                temp=0
            if (isnan(abs(complex(As[i]))) or isinf(abs(complex(As[i])))):
                As[i]=0 
            if (abs(As[i])<abs(temp) and abs(As[i])<abs(As[i-1]) ):
                As[i]=0
            Fs.append(F_jw.subs(t,xs[i]))
            if As[i]<0:
                As[i]=As[i]*(-1)
                As[i]=As[i]+np.pi
            if (abs(Fs[i])>2*np.pi):
                Fs[i]=(abs(Fs[i])/Fs[i])*(abs(Fs[i])%(2*np.pi))
        temp=A_jw.subs(t,xf[1])/(Period/2)
        for k in range(len(xf)-1):
            i=k+1
            a1.append(temp)
            if (i<len(xf)-1):
                temp=A_jw.subs(t,xf[i+1])/(Period/2)
            else:
                temp=0
            if (isnan(abs(complex(a1[i]))) or isinf(abs(complex(a1[i])))):
                a1[i]=0 
            if (abs(a1[i])<abs(temp) and abs(a1[i])<abs(a1[i-1]) ):
                a1[i]=0
            f1.append(F_jw.subs(t,xf[i]))
            if a1[i]<0:
                a1[i]=a1[i]*(-1)
            if (abs(f1[i])>2*np.pi):
                f1[i]=(abs(f1[i])/f1[i])*(abs(f1[i])%(2*np.pi))
            a2.append(a1[i]*np.abs(HS.subs(t,xf[i]*1.j)))
            f2.append(f1[i]+np.angle(complex(HS.subs(t,xf[i]*1.j))))
            j=j+1
            # Заполняем значения амплитудного спектра
        if (sflag==1):
            if(Stype==1):
                As[scount]=A*Tau/4
                Fs[scount]=-np.pi/2
            if(Stype==3):
                As[scount]=A*np.pi*np.pi/4
                Fs[scount]=0
        if (fflag==1):
            if(Stype==1):
                a1[fcount]=A*Tau/(2*Period)
                f1[fcount]=-np.pi/2
                a2[fcount]=a1[fcount]*np.abs(HS.subs(t,xf[fcount]*1.j))
                f2[fcount]=f1[fcount]+np.angle(complex(HS.subs(t,xf[fcount]*1.j)))
            if(Stype==3):
                a1[fcount]=A*np.pi*np.pi/(2*Period)
                f1[fcount]=0
                a2[fcount]=a1[fcount]*np.abs(HS.subs(t,xf[fcount]*1.j))
                f2[fcount]=f1[fcount]+np.angle(complex(HS.subs(t,xf[fcount]*1.j)))
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
        return xs,As,Fs,x,Y1,Y2"""
        

            

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
# f1f2=circ.FourierPeriod(Stype,hxy[3],A,Tau,Period)

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

# ФЧХ
plt.figure(10)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('ФЧХ')
plt.grid()
plt.plot(frq_a[2], frq_a[3])
plt.show()

# Амплитуда Фурье  IN
plt.figure(10)
plt.xlabel('Omega')
plt.ylabel('A')
plt.title('Ampletude Fourier IN')
plt.grid()
plt.stem(f1f2[6], f1f2[7])
plt.show()

# Фаза Фурье IN
plt.figure(10)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase Fourier IN')
plt.grid()
plt.stem(f1f2[6], f1f2[8])
plt.show()

# Амплитуда Фурье Out
plt.figure(10)
plt.xlabel('Omega')
plt.ylabel('A')
plt.title('Ampletude Fourier OUT')
plt.grid()
plt.stem(f1f2[6], f1f2[9])
plt.show()

# Фаза Фурье OUT
plt.figure(10)
plt.xlabel('Omega')
plt.ylabel('arg(A)')
plt.title('Phase Fourier OUT')
plt.grid()
plt.stem(f1f2[6], f1f2[10])
plt.show()

pass
