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
    def __init__(self, el_type, value, amperage, voltage ):
       # R - резистор
       # L - катушка
       # C - конденсатор
       # U - ИН
       # I - ИТ
       # NL - ХХ
       # SC - КЗ
       # el_type имеет тип строки

        self.el_type = el_type
       # Хранит в себе соответствующие значение для L, R, C элементов
        self.value = value
       # Обязательно определять для ИТ и ИН силу тока и напряжение соответственно
       # Для остальных типов желательно приравнимать к None
        self.amperage = amperage
        self.voltage = voltage

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
    def __init__(self, key):
        # Связанные с данным узлы, если левый или правый и средний совпадают, это обозначает край цепи
        self.left = None
        self.right = None
        self.mid = None
        # Массивы содержащии ключи элементов между данным узлом и правым, левым, средним соответственно
        self.right_array = None
        self.left_array = None
        self.mid_array = None
        # Ключ (номер) узла
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

        @property
        def key(self):
            return self._key

        @key.setter
        def key(self, key):
            self._key = key


class Circuit:
    def __init__(self, node_array, el_array):
        # Массив узлов (порядок не имеет значения)
        self.node_array = node_array
        # Массив элементов (номер в массиве являеться ключем)
        self.el_array = el_array
