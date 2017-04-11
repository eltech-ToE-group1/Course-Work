# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 17:13:37 2017

@author: Дудко Кирилл
"""

from numpy.linalg import inv
import numpy as np

class Voltage:
    def __init__(self, function, plus, minus):
        # Пока-что будем считать сигнал постоянным
        self.function = function
        # Номер узла, который соответствует плюсу
        self.plus = plus
        # Номер узла, который соответствует минусу
        self.minus = minus

        @property
        def function(self):
            return self._function

        @function.setter
        def function(self, function):
            self._function = function

        @property
        def plus(self):
            return self._plus

        @plus.setter
        def plus(self, plus):
            self._plus = plus

        @property
        def minus(self):
            return self._minus

        @minus.setter
        def minus(self, minus):
            self._minus = minus


class Amperage:
    def __init__(self, function, from_this, to_this):
        # Пока-что будем считать сигнал постоянным
        self.function = function
        # Номер узла, откуда течёт ток
        self.from_this= from_this
        # Номер узла, куда течёт ток
        self.to_this = to_this

        @property
        def function(self):
            return self._function

        @function.setter
        def function(self, function):
            self._function = function

        @property
        def from_this(self):
            return self._from_this

        @from_this.setter
        def from_this(self, from_this):
            self._from_this = from_this

        @property
        def to_this(self):
            return self._to_this

        @to_this.setter
        def to_this(self, to_this):
            self._to_this = to_this

class Element:
    def __init__(self, key , el_type, value, amperage, voltage, from_, to_ ):
        # R - резистор
        # L - катушка
        # C - конденсатор
        # U - ИН
        # I - ИТ
        # NL - ХХ
        # SC - КЗ
        # el_type имеет тип строки

        self.el_type = el_type
        # Хранит в себе номер элемента
        self.key = key
        # Хранит в себе соответствующие значение для L, R, C элементов
        self.value = value
        # Обязательно определять для ИТ и ИН силу тока и напряжение соответственно
        # Для остальных типов приравнимать function к None и хранить только направления
        if(amperage):
            self.amperage = Amperage(amperage,from_,to_)
        else:
            self.amperage=None
        if(voltage):
            self.voltage = Voltage(voltage,to_,from_)
        else:
            self.voltage = None
        # Хранить направление тока, если сила тока неизвестна
        self.from_=from_
        self.to_=to_

        @property
        def el_type(self):
            return self._el_type

        @el_type.setter
        def el_type(self, el_type):
            self._el_type = el_type

        @property
        def value (self):
            return self._value

        @value .setter
        def value (self, value ):
            self._value  = value

        @property
        def amperage(self):
            return self._amperage

        @amperage.setter
        def amperage(self, amperage):
            self._amperage = amperage

        @property
        def voltage(self):
            return self._voltage

        @voltage.setter
        def voltage(self, voltage):
            self._voltage = voltage


class Node:
    def __init__(self,key):
        # Массивы содержащии ключи элементов 
        self.To = None
        self.From = None
        self.key = key

        @property
        def left(self):
            return self._left

        @left.setter
        def left(self, left):
            self._left = left

        @property
        def right(self):
            return self._right

        @right .setter
        def right(self, right):
            self._right  = right

        @property
        def mid(self):
            return self._mid

        @mid.setter
        def mid(self, mid):
            self._mid = mid

        @property
        def right_array(self):
            return self._right_array

        @right_array.setter
        def right_array(self, right_array):
            self._right_array = right_array

        @property
        def mid_array(self):
            return self._mid_array

        @mid_array.setter
        def mid_array(self, mid_array):
            self._mid_array = mid_array


class Circuit:
    def __init__(self, node_array, el_array):
        # Массив узлов (номер в массиве является ключем)
        self.node_array = node_array
        # Массив элементов (номер в массиве является ключем)
        self.el_array = el_array
        for i in el_array:
            if(node_array[i.from_].To):
                node_array[i.from_].To.append(i.key)
            else:
                node_array[i.from_].To=[i.key]
            if(node_array[i.to_].From):
                node_array[i.to_].From.append(i.key)
            else:
                node_array[i.to_].From=[i.key]
            

    def NORMALIZE(self):
        for i in self.el_array:
            if i.el_type.find('SC'):
                self.node_array[i.from_].To.extend(self.node_array[i.to_].To)
                self.node_array[i.from_].From.extend(self.node_array[i.to_].From)
                self.node_array[i.from_].el_To.extend(self.node_array[i.to_].el_To)
                self.node_array[i.from_].el_From.extend(self.node_array[i.to_].el_From)
                self.node_array.pop([i.to_])
                i.el_type='N'
                self.el_array.remove(i)
            #if i.el_type.find('NL'):
                
               


    def MUN(self,elem): # На вход подаётся элемент 
        #self.NORMILIZE #-когда будет готов ХХ убрать коммент
        N=0
        for i in self.node_array:
            N=N+1
        # Задаем матрицы проводимойстей и токов
        print(N)
        Basic=0
        G=np.zeros((N,N))
        I=np.zeros((N))
        # Находим базовый узел
        for i in self.el_array:
            if (i.el_type.find('U')!=-1):
                Basic=i.from_
        # Заполняем собственные проводимости узлов
        for i in range(N):
            if(i!=Basic):
                if(self.node_array[i].To):
                    for j in self.node_array[i].To:
                        G[i][i]=G[i][i]+1/self.el_array[j].value
                if(self.node_array[i].From):
                    for j in self.node_array[i].From:
                        G[i][i]=G[i][i]+1/self.el_array[j].value
        # Заполняем матрицу проводимойстей до конца
        for i in range(N):
            if(i!=Basic):
                if(self.node_array[i].To):
                    for j in self.node_array[i].To:
                        if (self.el_array[j].el_type.find('R')!=-1):
                            G[i][self.el_array[j].to_]=G[i][self.el_array[j].to_]-1/self.el_array[j].value
                if(self.node_array[i].From):
                    for j in self.node_array[i].From:
                        if (self.el_array[j].el_type.find('R')!=-1):
                            G[i][self.el_array[j].from_]=G[i][self.el_array[j].from_]-1/self.el_array[j].value
        # Заполняем матрицу токов
        for i in self.el_array:
            if (i.el_type.find('I')!=-1):
                if(i.amperage.to_this!=Basic):
                    I[i.amperage.to_this]=I[i.amperage.to_this]+i.amperage.function
                if(i.amperage.from_this!=Basic):
                    I[i.amperage.from_this]=I[i.amperage.from_this]-i.amperage.function
            if (i.el_type.find('U')!=-1): # Если у нас есть ИН, то соответствующию строку в матрице проводимостей переводим в единичную, а соответсвующее значение матрицы токов приравниваем напряжению
                for j in range(N-1):
                    G[i.to_][j]=0
                    G[i.from_][j]=0
                if (i.voltage.plus!=Basic):
                    G[i.to_][i.voltage.plus]=1
                    I[i.to_]=i.voltage.function
                if (i.voltage.minus!=Basic):
                    G[i.from_][i.voltage.minus]=1
                    I[i.from_]=-1*i.voltage.function
        #print(I)
        N=N-1;
        F=np.zeros((N,N))
        k=0;
        n=0;
        for i in range(N):
            if (i!=Basic):
                if(i<Basic):
                    k=i;
                else:
                    k=i-1
            for j in range(N):
                if (j!=Basic):
                    if(j<Basic):
                        n=j;
                    else:
                        n=j-1
                    F[k][n]=G[i][j]
        II=np.zeros((N))
        for i in range(N):
            if (i!=Basic):
                if(i<Basic):
                    k=i;
                else:
                    k=i+1
                II[k]=I[i]
        V=inv(F)*II #получаем вектор узловых напряжений.
        V1=np.zeros(N+1)
        for i in range(N+1):
            if (i!=Basic):
                if(i<Basic):
                    k=i;
                else:
                    k=i+1
                print(k,i,V[i])
                V1[k]=V[i][0]
        V1[Basic]=0;
        return (V1[elem.from_]-V1[elem.to_]) #возвращаем разницу узловых напряжений на элементе, то есть его напряжение.


#Зададим цепь
node_array=[Node(0),Node(1),Node(2)]
node_elem=[Element(0,"U",5,None,5,2,0),Element(1,"R",1,None,None,0,1),Element(2,"R",2,None,None,1,2),Element(3,"R",3,None,None,1,2)]
circ=Circuit(node_array,node_elem)
print(circ.MUN(node_elem[1]))
