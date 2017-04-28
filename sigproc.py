import circuit as cir
from sympy.abc import  s, t
from sympy import Heaviside, sin, cos, exp
from numpy import pi
from sympy.integrals import inverse_laplace_transform
from sympy.parsing.sympy_parser import parse_expr
from scipy.fftpack import  fft, fftshift
from scipy.signal import ss2tf
import matplotlib.pyplot as plt
import numpy as np


class SignalCircuit(cir.Circuit):

    def __init__(self, node_array, el_array, sigtype, A, Tau, signal = 0):
        cir.Circuit.__init__(self, node_array, el_array)
        if sigtype == 1:
            self.signal = A*(sin(t*pi/Tau) - sin(t*pi/Tau)*Heaviside(t-Tau,1))
        if sigtype == 2:
            self.signal = 2*A/Tau*(t-t*Heaviside(t - Tau/2, 1))+(-2*A/Tau*t +2*A)*(Heaviside(t- Tau/2, 1)-Heaviside(t - Tau, 1))
        if sigtype == 3:
            self.signal = A*(cos(t*pi/Tau) - cos(t*pi/Tau)*Heaviside(t-Tau,1))
        if sigtype == 4:
            self.signal = A*(1 - 2*Heaviside(t-Tau/2),1)
        # Сигнал, заданный пользователем
        if sigtype == 5:
            self.signal = signal
        # Время сигнала
        self.Tau = Tau

    def H_t(self):
        ABCD = cir.Circuit.StateSpace(self, 5, 'U')
        H_S=ss2tf(ABCD[0],ABCD[1],ABCD[2],ABCD[3])
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
        return (inverse_laplace_transform(HS,s,t))

    def H1_t(self):
        ABCD = cir.Circuit.StateSpace(self, 5, 'U')
        H_S=ss2tf(ABCD[0],ABCD[1],ABCD[2],ABCD[3])
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
        HS = parse_expr(HS, evaluate=False)/s
        return (inverse_laplace_transform(HS,s,t))

node_array = [cir.Node(0), cir.Node(1), cir.Node(2), cir.Node(3)]
node_elem = [cir.Element(0, "I", 1, 1, None, 0, 1), cir.Element(1, "R", 2, None, None, 1, 0),
             cir.Element(2, "L", 2, None, None, 1, 2), cir.Element(3, "R", 1, None, None, 2, 3),
             cir.Element(4, "C", 4, None, None, 3, 0), cir.Element(5, "R", 0.5, None, None, 3, 0)]
circ = SignalCircuit(node_array, node_elem, 1, 1, 1)
a = circ.H_t()
b = circ.H1_t()
x = np.arange(0.01, pi, 0.01)
y = []
for i in x:
    y.append(b.subs(t, i))
plt.plot(x,y)
plt.grid()
plt.show()
pass