import circuit as cir
from scipy.signal import ss2tf
import numpy as np
import copy

# На вход подаёться цепь и номер элемента, также тип реакции на элементе, относительно которого считаються матрицы C и D
def StateSpace(source_cir, el_num, f2_type):
    # Список реактивных элементов
    react_list = source_cir.FindReact()
    # Список источников сигнала
    source_list = source_cir.FindSource()
    # Инициируем матрицы уравнения состояния
    A = np.zeros(len(react_list), len(react_list))
    B = np.zeros(len(react_list), len(source_list))
    C = np.zeros(1, len(react_list))
    D = np.zeros(1, len(source_list))
    n_cir = copy.deepcopy(source_cir)
    # Проходимся по реактивным элементам
    for i in react_list:
        # Заменяем индуктивность на ИТ с током 1
        if n_cir.el_array[i].el_type == 'L':
            n_cir.el_array[i].el_type = 'I'
            n_cir.el_array[i].value = 1
            n_cir.el_array[i].amperage = cir.Amperage(1, n_cir.el_array[i].from_, n_cir.el_array[i].to_ )
        # Заменяем конденсаторы на ИН с напряжением 1
        if n_cir.el_array[i].el_type == 'C':
            n_cir.el_array[i].el_type = 'U'
            n_cir.el_array[i].value = 1
            n_cir.el_array[i].voltage = cir.Voltage(1, n_cir.el_array[i].to_, n_cir.el_array[i].from_)
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
                B[j][i] = n2_cir.el_array[react_list[j]].voltage.function/source_cir.el_array[react_list[j]].value
            if n2_cir.el_array[react_list[j]].el_type == 'SC':
                B[j][i] = n2_cir.el_array[react_list[j]].amperage.function/source_cir.el_array[react_list[j]].value
        if f2_type == 'U':
            D[1][i] = n2_cir.el_array[el_num].voltage.function
        if f2_type == 'I':
            D[1][i] = n2_cir.el_array[el_num].amperage.function
        del n2_cir

    # Заполняем матрицы A и C
    for i in range(len(react_list)):
        n2_cir = copy.deepcopy(n_cir)
        n2_cir.change_sources(react_list[i])
        if n2_cir.el_array[react_list[i]].el_type == 'U':
            n2_cir.el_array[react_list[i]].from_, n2_cir.el_array[react_list[i]].to_ = \
                n2_cir.el_array[react_list[i]].to_, n2_cir.el_array[react_list[i]].from_
            n2_cir.Refresh()
        n2_cir.solve()
        for j in range(len(react_list)):
            if j == i:
                if n2_cir.el_array[react_list[j]].el_type == 'U':
                    A[j][i] = -n2_cir.el_array[react_list[j]].voltage.function/source_cir.el_array[react_list[j]].value
            if n2_cir.el_array[react_list[j]].el_type == 'I':
                A[j][i] = n2_cir.el_array[react_list[j]].voltage.function/source_cir.el_array[react_list[j]].value
            if n2_cir.el_array[react_list[j]].el_type == 'SC':
                A[j][i] = n2_cir.el_array[react_list[j]].amperage.function/source_cir.el_array[react_list[j]].value
        if f2_type == 'U':
            C[1][i] = n2_cir.el_array[el_num].voltage.function
        if f2_type == 'I':
            C[1][i] = n2_cir.el_array[el_num].amperage.function
        del n2_cir
    return A, B, C, D


