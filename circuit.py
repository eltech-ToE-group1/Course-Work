from numpy.linalg import inv
import numpy as np
import copy
from scipy.signal import ss2tf

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
        self.from_this = from_this
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
    def __init__(self, key, el_type, value, amperage, voltage, from_, to_):
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
        if el_type == 'I' or amperage:
            self.amperage = Amperage(amperage, from_, to_)
        else:
            self.amperage = None
        if el_type == 'U' or voltage:
            self.voltage = Voltage(voltage, to_, from_)
        else:
            self.voltage = None
        # Хранить направление тока, если сила тока неизвестна
        self.from_ = from_
        self.to_ = to_

        @property
        def key(self):
            return self.key

        @key.setter
        def key(self, key):
            self.key = key

        @property
        def el_type(self):
            return self._el_type

        @el_type.setter
        def el_type(self, el_type):
            self._el_type = el_type

        @property
        def value(self):
            return self._value

        @value.setter
        def value(self, value):
            self._value = value

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
            self.from_ = from_


class Node:
    def __init__(self, key):
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
            self.From = From

        @property
        def key(self):
            return self.key

        @key.setter
        def key(self, key):
            self.key = key


class Circuit:
    def __init__(self, node_array, el_array):
        # Массив узлов (номер в массиве является ключем)
        self.node_array = node_array
        # Массив элементов (номер в массиве является ключем)
        self.el_array = el_array
        for i in el_array:
            if (node_array[i.from_].To):
                node_array[i.from_].To.append(i.key)
            else:
                node_array[i.from_].To = [i.key]
            if (node_array[i.to_].From):
                node_array[i.to_].From.append(i.key)
            else:
                node_array[i.to_].From = [i.key]

    def FindReact(self):
        Source_arr = []
        for i in self.el_array:
            if (i.el_type == 'L' or i.el_type == 'C'):
                Source_arr.append(i.key)
        return Source_arr

    def FindSource(self):
        Source_arr = []
        for i in self.el_array:
            if (i.el_type == 'U' or i.el_type == 'I'):
                Source_arr.append(i.key)
        return Source_arr

    # ищет элемент заданного типа, если он не замкнут на одном узле
    def FindEl(self, eltype):
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
        if not (num + 1):
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
        N = 0
        for i in self.node_array:
            N = N + 1
        # Задаем матрицы проводимойстей и токов
        Basic = 0
        G = np.zeros((N, N))
        I = np.zeros((N))
        # Находим базовый узел
        for i in self.el_array:
            if (i.el_type.find('U') != -1):
                Basic = i.from_
        # Заполняем собственные проводимости узлов
        for i in range(N):
            if (i != Basic):
                if (self.node_array[i].To):
                    for j in self.node_array[i].To:
                        if (self.el_array[j].el_type.find('R') != -1):
                            G[i][i] = G[i][i] + 1 / self.el_array[j].value
                if (self.node_array[i].From):
                    for j in self.node_array[i].From:
                        if (self.el_array[j].el_type.find('R') != -1):
                            G[i][i] = G[i][i] + 1 / self.el_array[j].value
        # Заполняем матрицу проводимойстей до конца
        for i in range(N):
            if (i != Basic):
                if (self.node_array[i].To):
                    for j in self.node_array[i].To:
                        if (self.el_array[j].el_type.find('R') != -1):
                            G[i][self.el_array[j].to_] = G[i][self.el_array[j].to_] - 1 / self.el_array[j].value
                if (self.node_array[i].From):
                    for j in self.node_array[i].From:
                        if (self.el_array[j].el_type.find('R') != -1):
                            G[i][self.el_array[j].from_] = G[i][self.el_array[j].from_] - 1 / self.el_array[j].value
        # Заполняем матрицу токов
        for i in self.el_array:
            if (i.el_type.find('I') != -1):
                if (i.amperage.to_this != Basic):
                    I[i.amperage.to_this] = I[i.amperage.to_this] + i.amperage.function
                if (i.amperage.from_this != Basic):
                    I[i.amperage.from_this] = I[i.amperage.from_this] - i.amperage.function
            if (i.el_type.find(
                    'U') != -1):  # Если у нас есть ИН, то соответствующию строку в матрице проводимостей переводим в единичную, а соответсвующее значение матрицы токов приравниваем напряжению
                for j in range(N - 1):
                    G[i.to_][j] = 0
                    G[i.from_][j] = 0
                if (i.voltage.plus != Basic):
                    G[i.to_][i.voltage.plus] = 1
                    I[i.to_] = i.voltage.function
                if (i.voltage.minus != Basic):
                    G[i.from_][i.voltage.minus] = 1
                    I[i.from_] = -1 * i.voltage.function
        N1 = N - 1
        F = np.zeros((N1, N1))
        k = 0
        for i in range(N1):
            if (k == Basic):
                k = k + 1
            n = 0
            for j in range(N1):
                if (n == Basic):
                    n = n + 1
                F[i][j] = G[k][n]
                n = n + 1
            k = k + 1
        II = np.zeros((N1))
        k = 0
        for i in range(N1):
            if (k == Basic):
                k = k + 1
            II[i] = I[k]
            k = k + 1
        V = np.matmul(inv(F), II)  # получаем вектор узловых напряжений.
        # print(inv(F))
        # print(II)
        # print(V)
        V1 = np.zeros(N)
        k = 0
        for i in range(N1):
            while (k == Basic):
                k = k + 1
            V1[k] = V[i]
            k = k + 1
        V1[Basic] = 0
        # print(V1)
        for i in self.el_array:
            if i.el_type.find('U') == -1:
                i.voltage = Voltage(V1[i.from_] - V1[i.to_], i.to_, i.from_)
                if i.el_type == 'R':
                    i.amperage = Amperage((V1[i.from_] - V1[i.to_]) / i.value, i.from_, i.to_)

    # Копирует значения сил токов и напряжений в self из n_cir
    def va_cpy(self, n_cir):
        for i in range(len(self.el_array)):
            if n_cir.el_array[i].amperage is not None:
                self.el_array[i].amperage = Amperage(copy.deepcopy(n_cir.el_array[i].amperage.function),
                                                         None, None)
            if n_cir.el_array[i].voltage is not None:
                self.el_array[i].voltage = Voltage(copy.deepcopy(n_cir.el_array[i].voltage.function),
                                                       None, None)
            # Проставляем напрявления
            self.Refresh()

    # Перед вызовом этой функции нужно расчитать цепь в МУН  и "развернуть" её если там были свёрнутые КЗ
    # Данная функция расчитывает неизвестные токи, используя закон токов Кирхгофа
    def KCL(self):
        # Все токи расчитаны
        all_amp = True
        # Проверка, расчитаны ли все токи
        for i in self.el_array:
            if i.amperage is None:
                all_amp = False
        # Пока все токи не расчитаны
        while not all_amp:
            # Проходим по всем узлам
            for i in self.node_array:
                # Число неизвестных в узле
                num = 0
                # Номер неизвестного элемента
                un = 0
                for j in i.From:
                    if self.el_array[j].amperage is None:
                        num = num + 1
                        un = j
                for j in i.To:
                    if self.el_array[j].amperage is None:
                        num = num + 1
                        un = j
                # Если неизвестный только 1
                if num == 1:
                    # Входящий ток
                    current_in = 0
                    # Выходящий ток
                    current_out = 0
                    for j in i.From:
                        if self.el_array[j].amperage is not None:
                            current_in = current_in + self.el_array[j].amperage.function
                    for j in i.To:
                        if self.el_array[j].amperage is not None:
                            current_out = current_out + self.el_array[j].amperage.function
                    cur = abs(current_in - current_out)
                    self.el_array[un].amperage = Amperage(cur, self.el_array[un].from_,
                                                              self.el_array[un].to_)
                all_amp = True
                for j in self.el_array:
                    if j.amperage is None:
                        all_amp = False

    # Функция, решающая цепь, если в ней отсутствуют реактивные элементы и не боле одного ИН
    def solve(self):
        if len(self.FindReact()) > 0:
            return
        source_arr = self.FindSource()
        # Количество ИН
        num = 0
        for i in source_arr:
            if self.el_array[i].el_type == 'U':
                num = num + 1
        if num > 1:
            return
        n_cir = copy.deepcopy(self)
        # Пока есть не замкнутые КЗ
        while n_cir.FindEl('SC') + 1:
            n_cir = n_cir.NodeMerge()
        n_cir.MUN()
        # Копируем данные из n_cir
        self.va_cpy(n_cir)
        # Окончательно решаем цепь
        self.KCL()
        del n_cir

    # Исключает все источники, за исключением одного, чей номер подан в качестве аргумента
    def change_sources(self, num):
        for i in range(len(self.el_array)):
            if i != num:
                if self.el_array[i].el_type == 'I':
                    self.el_array[i].value = 0
                    if self.el_array[i].amperage is not None:
                        self.el_array[i].amperage.function = 0
                if self.el_array[i].el_type == 'U':
                    self.el_array[i].el_type = 'SC'
                    self.el_array[i].amperage = None
                    self.el_array[i].voltage = None

    # Изменяет полярность элемента с выборанным номером
    def change_polarity(self, num):
        nodeto = self.el_array[num].to_
        self.node_array[nodeto].From.remove(num)
        self.node_array[nodeto].To.append(num)
        nodefrom = self.el_array[num].from_
        self.node_array[nodefrom].To.remove(num)
        self.node_array[nodefrom].From.append(num)
        self.Refresh()

    # На вход подаёться цепь и номер элемента, также тип реакции на элементе, относительно которого считаються матрицы C и D
    def StateSpace(self, el_num, f2_type):
        # Список реактивных элементов
        react_list = self.FindReact()
        # Список источников сигнала
        source_list = self.FindSource()
        # Инициируем матрицы уравнения состояния
        Num1 = len(react_list)
        Num2 = len(source_list)
        A = np.zeros((Num1, Num1))
        B = np.zeros((Num1, Num2))
        C = np.zeros((1, Num1))
        D = np.zeros((1, Num2))
        n_cir = copy.deepcopy(self)
        # Проходимся по реактивным элементам
        for i in react_list:
            # Заменяем индуктивность на ИТ с током 1
            if n_cir.el_array[i].el_type == 'L':
                n_cir.el_array[i].el_type = 'I'
                n_cir.el_array[i].value = 1
                n_cir.el_array[i].amperage = Amperage(1, n_cir.el_array[i].from_, n_cir.el_array[i].to_)
            # Заменяем конденсаторы на ИН с напряжением 1
            if n_cir.el_array[i].el_type == 'C':
                n_cir.el_array[i].el_type = 'U'
                n_cir.el_array[i].value = 1
                n_cir.el_array[i].voltage = Voltage(1, n_cir.el_array[i].to_, n_cir.el_array[i].from_)
        # Проходимся по источникам
        for i in source_list:
            # Заменяем силу тока в ИТ
            if n_cir.el_array[i].el_type == 'I':
                n_cir.el_array[i].value = 1
                n_cir.el_array[i].amperage.function = 1
            # Заменяем напряжение в ИН
            if n_cir.el_array[i].el_type == 'U':
                n_cir.el_array[i].value = 1
                n_cir.el_array[i].voltage.function = 1
        # Заполняем матрицы B и D
        for i in range(len(source_list)):
            n2_cir = copy.deepcopy(n_cir)
            n2_cir.change_sources(source_list[i])
            n2_cir.solve()
            for j in range(len(react_list)):
                if n2_cir.el_array[react_list[j]].el_type == 'I':
                    B[j][i] = n2_cir.el_array[react_list[j]].voltage.function / self.el_array[react_list[j]].value
                if n2_cir.el_array[react_list[j]].el_type == 'SC':
                    B[j][i] = n2_cir.el_array[react_list[j]].amperage.function / self.el_array[
                        react_list[j]].value
            if f2_type == 'U':
                D[0][i] = n2_cir.el_array[el_num].voltage.function
            if f2_type == 'I':
                D[0][i] = n2_cir.el_array[el_num].amperage.function
            del n2_cir
        # Заполняем матрицы A и C
        for i in range(len(react_list)):
            n2_cir = copy.deepcopy(n_cir)
            n2_cir.change_sources(react_list[i])
            if n2_cir.el_array[react_list[i]].el_type == 'U':
                # Меняем полярность для правильного расчёта C элемента
                n2_cir.change_polarity(react_list[i])
            n2_cir.solve()
            for j in range(len(react_list)):
                if n2_cir.el_array[react_list[j]].el_type == 'U':
                    A[j][i] = -(
                        n2_cir.el_array[react_list[j]].amperage.function / self.el_array[react_list[j]].value)
                if n2_cir.el_array[react_list[j]].el_type == 'I':
                    A[j][i] = n2_cir.el_array[react_list[j]].voltage.function / self.el_array[react_list[j]].value
                if n2_cir.el_array[react_list[j]].el_type == 'SC':
                    A[j][i] = n2_cir.el_array[react_list[j]].amperage.function / self.el_array[
                        react_list[j]].value
            if f2_type == 'U':
                C[0][i] = n2_cir.el_array[el_num].voltage.function
            if f2_type == 'I':
                C[0][i] = n2_cir.el_array[el_num].amperage.function
            del n2_cir
        return A, B, C, D

# Зададим цепь
"""node_array = [Node(0), Node(1), Node(2), Node(3)]
node_elem = [Element(0, "I", 1, 1, None, 0, 1), Element(1, "R", 2, None, None, 1, 0),
             Element(2, "L", 2, None, None, 1, 2), Element(3, "R", 1, None, None, 2, 3),
             Element(4, "C", 4, None, None, 3, 0), Element(5, "R", 0.5, None, None, 3, 0)]
circ = Circuit(node_array, node_elem)
a = circ.StateSpace(5, 'U')
HS = ss2tf(a[0], a[1], a[2], a[3])"""
pass
