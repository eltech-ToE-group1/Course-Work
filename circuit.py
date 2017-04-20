# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 17:13:37 2017

@author: Дудко Кирилл
"""

from numpy.linalg import inv
import numpy as np
import copy

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
        # SC - КЗ
        # ХХ представляеться как ИТ с нулевым током
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
        def key(self):
            return self.key
        
        @key.setter
        def key(self, key):
            self.key  = key
            
            
        @property
        def el_type(self):
            return self._el_type

        @el_type.setter
        def el_type(self, el_type):
            self._el_type = el_type

        @property
        def value (self):
            return self._value

        @value.setter
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
        @property    
        def to_(self):
            return self.to_

        @to_.setter
        def to_(self, To):
            self.to_ = To

        @property
        def from_(self):
            return self.from_

        @from_.setter
        def from_(self, from_):
            self.from_  = from_      
    
class Node:
    def __init__(self,key):
        # Массивы содержащии ключи элементов 
        self.To = None
        self.From = None
        self.key = key

        @property
        def To(self):
            return self.To

        @To.setter
        def To(self, To):
            self.To = To

        @property
        def From(self):
            return self.From

        @From.setter
        def From(self, From):
            self.From  = From

        @property
        def key(self):
            return self.key
        
        @key.setter
        def key(self, key):
            self.key  = key
    

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
    
    def FindReact(self):
        f=0
        for i in self.el_array:
            if (i.el_type.find('R')!=-1 or i.el_type.find('L')!=-1 or i.el_type.find('C')!=-1):
                if(f==0):
                    React_arr=[i]
                else:
                    React_arr.append(i)
        return React_arr
    
    def FindSource(self):
        f=0
        for i in self.el_array:
            if (i.el_type.find('U')!=-1 or i.el_type.find('I')!=-1):
                if(f==0):
                    Source_arr=[i]
                else:
                    Source_arr.append(i)
        return Source_arr
    
    # ищет элемент заданного типа, если он не замкнут на одном узле
    def FindEl (self, eltype):
        for i in range(len(self.el_array)):
            if self.el_array[i].to_ != self.el_array[i].from_:
                if self.el_array[i].el_type == eltype:
                    return i
        return -1
    
    # Меняет направления в напряжении и токах элементов в соответствии с from_ и to_, а также заменяет ключи узлов на их номера в массиве и соответственно меняет from_ и to_
    def Refresh(self):
        for i in range(len(self.node_array)):
            for j in self.node_array[i].From:
                self.el_array[j].to_ = i
            for j in self.node_array[i].To:
                self.el_array[j].from_ = i
            self.node_array[i].key = i
        for i in self.el_array:
            if (i.amperage):
                i.amperage.from_this = i.from_
                i.amperage.to_this = i.to_
            if (i.voltage):
                i.voltage.plus = i.to_
                i.voltage.minus = i.from_

    # На вход подаёться цепь и номера узлов, если между узлами КЗ, они объеденяються в один, функция возвращает цепь.
    def NodeMerge(self):
        new_circ = copy.deepcopy(self)
        num = new_circ.FindEl('SC')
        # Если КЗ нету
        if not (num+1):
            return new_circ
        node1 = new_circ.el_array[num].from_
        node2 = new_circ.el_array[num].to_
        # Объеденяем массивы элементов
        new_circ.node_array[node2].From.extend(self.node_array[node1].From)
        new_circ.node_array[node2].To.extend(self.node_array[node1].To)
        # Меняем направления в элементах
        for i in new_circ.node_array[node2].From:
            new_circ.el_array[i].from_ = node2
        for i in new_circ.node_array[node2].To:
            new_circ.el_array[i].to_ = node2
        # Удаляем лишний узел
        del new_circ.node_array[node1]
        # Обновляем цепь
        new_circ.Refresh()
        return new_circ


     def MUN(self): 
        N=0
        for i in self.node_array:
            N=N+1
        # Задаем матрицы проводимойстей и токов
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
                        if (self.el_array[j].el_type.find('R')!=-1):
                            G[i][i]=G[i][i]+1/self.el_array[j].value
                if(self.node_array[i].From):
                    for j in self.node_array[i].From:
                        if (self.el_array[j].el_type.find('R')!=-1):
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
        N1=N-1;
        F=np.zeros((N1,N1))
        k=0
        for i in range(N1):
            if (k==Basic):
                k=k+1;  
            n=0;
            for j in range(N1):
                if(n==Basic):
                    n=n+1
                F[i][j]=G[k][n]
                n=n+1;
            k=k+1;
        II=np.zeros((N1))
        k=0;
        for i in range(N1):
            if(k==Basic):
                k=k+1;
            II[i]=I[k]
            k=k+1
        V=np.matmul(inv(F),II) #получаем вектор узловых напряжений.
        print(inv(F))
        print(II)
        print(V)
        V1=np.zeros(N)
        k=0;
        for i in range(N1):
            while(k==Basic):
                k=k+1;
            V1[k]=V[i]
            k=k+1;
        V1[Basic]=0
        print(V1)
        for i in self.el_array:
            if i.el_type.find('U')==-1:
                i.voltage = Voltage(V1[i.from_]-V1[i.to_], i.to_, i.from_)
                if i.el_type == 'R':
                    i.amperage = Amperage((V1[i.from_]-V1[i.to_])/i.value, i.from_, i.to_)


#Зададим цепь
node_array=[Node(0),Node(1),Node(2),Node(3)]
node_elem=[Element(0,"I",1,1,None,3,0),Element(1,"R",3,None,None,0,1),Element(2,"R",2,None,None,1,2),Element(3,"SC",None,None,None,1,3),Element(4,"R",4,None,None,2,3),Element(5,"R",5,None,None,2,3)]
circ=Circuit(node_array,node_elem)
#circ.MUN()
n_circ = circ.NodeMerge()
#n_circ.MUN()
pass

