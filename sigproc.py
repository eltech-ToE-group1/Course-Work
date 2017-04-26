import circuit as cir
from sympy.abc import  s, t
from sympy import Heaviside, sin, cos, exp
from numpy import pi
from sympy.integrals import inverse_laplace_transform
from scipy.fftpack import  fft, fftshift
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