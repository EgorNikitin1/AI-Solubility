from random import choices
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def check():  # Проверка на растворимость
    weight = ()
    flag = True
    if solubility == "Р":
        print("Вещество растворимо")
        weight = (50, 50, 0)
    elif solubility == "М":
        print("Вещество мало растворимо")
        weight = (10, 10, 80)
    elif solubility == "Н":
        print("Вещество не растворимо")
        weight = (2, 2, 96)
    elif solubility == "-":
        print("Вещество не растворяется в водной среде")
        flag = False
    elif solubility == "?":
        print("Нет достоверных сведений о существовании соединения")
        flag = False
    return weight, flag


def randomizer(p):  # Рандомайзер
    ch = choices((True, False, "Stop"), weights=p, k=1)
    return ch[0]


def rotation(temp):  # Поворот
    r = randomizer(check_tuple[0])
    if r is True:  # По часовой
        rotated = list(zip(*temp[::-1]))
    elif r is False:  # Против часовой
        rotated = list(zip(*temp))[::-1]
    else:
        rotated = temp
    return rotated


def algorithm(k):  # Алгоритм окрестности Марголуса
    for i in range(1, n, 2):
        for j in range(1, m, 2):
            temp = []
            if (not i == n - 1 and not j == m - 1) and k == 1:
                temp.append([data[i][j], data[i][j + k]])
                temp.append([data[i + k][j], data[i + k][j + k]])
                rotated = rotation(temp)
                data[i][j], data[i][j + k], data[i + k][j], data[i + k][j + k] = rotated[0][0], rotated[0][1], rotated[1][0], rotated[1][1]
            elif k == -1:
                temp.append([data[i][j], data[i][j + k]])
                temp.append([data[i + k][j], data[i + k][j + k]])
                rotated = rotation(temp)
                data[i][j], data[i][j + k], data[i + k][j], data[i + k][j + k] = rotated[-1][-1], rotated[-1][0], rotated[0][-1], rotated[0][0]


def space(n, m):  # Создание поля
    data = [[0 for _ in range(n)] for _ in range(m)]
    for i in range(n):
        for j in range(m):
            if n // 2 - k < i < n // 2 + k and n // 2 - k < j < n // 2 + k:
                data[i][j] = 1
    return data


# Main
n, m = 200, 200
k = 20
data = space(n, m)
chart = pd.read_excel('Solubility Chart.xlsx', index_col=0)
kation, anion = input("Введите катион вещества: "), input("Введите анион вещества: ")
solubility = chart.loc[anion, kation]
check_tuple = check()
epoch = 100

plt.ion()
fig, ax = plt.subplots()
axim = ax.imshow(np.array(data), vmin=0, vmax=1, cmap='cool')

if check_tuple[1]:
    for _ in range(epoch):
        algorithm(1)
        algorithm(-1)
        axim.set_data(np.array(data))
        fig.canvas.flush_events()
        # plt.pause(0.1)

    plt.savefig("Solubility.png")
    plt.show()
