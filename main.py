from random import choices
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

n, l, P = 200, 20, 0
data, weight = [], []
solubility, kation, anion, mes = "", "", "", ""


def check():  # Проверка на растворимость
    global mes
    p = 0
    flag = True
    if solubility == "Р":
        mes = "Вещество растворимо"
        p = 0
    elif solubility == "М":
        mes = "Вещество мало растворимо"
        p = 50
    elif solubility == "Н":
        mes = "Вещество не растворимо"
        p = 90
    elif solubility == "-":
        mes = "Вещество не растворяется в водной среде"
        flag = False
    elif solubility == "?":
        mes = "Нет достоверных сведений о существовании соединения"
        flag = False
    else:
        flag = False
    return p, flag


def rotation(temp):  # Поворот
    global weight

    def randomizer(p):  # Рандомайзер
        ch = choices((True, False, "Stop"), weights=p, k=1)
        return ch[0]

    r = randomizer(weight)
    if r is True:  # По часовой
        rotated = list(zip(*temp[::-1]))
    elif r is False:  # Против часовой
        rotated = list(zip(*temp))[::-1]
    else:
        rotated = temp
    return rotated


def algorithm_margolus(k):  # Алгоритм окрестности Марголуса
    global weight
    for i in range(1, n, 2):
        for j in range(1, n, 2):
            temp = []
            if (not i == n - 1 and not j == n - 1) and k == 1:  # Четное
                counted = algorithm_moore_even(i, j)
                weight = set_weight(P, counted)
                temp.append([data[i][j], data[i][j + k]])
                temp.append([data[i + k][j], data[i + k][j + k]])
                rotated = rotation(temp)
                data[i][j], data[i][j + k], data[i + k][j], data[i + k][j + k] = rotated[0][0], rotated[0][1], rotated[1][0], rotated[1][1]
            elif k == -1:  # Нечетное
                counted = algorithm_moore_odd(i, j)
                weight = set_weight(P, counted)
                temp.append([data[i][j], data[i][j + k]])
                temp.append([data[i + k][j], data[i + k][j + k]])
                rotated = rotation(temp)
                data[i][j], data[i][j + k], data[i + k][j], data[i + k][j + k] = rotated[-1][-1], rotated[-1][0], rotated[0][-1], rotated[0][0]


def algorithm_moore_even(row, col):  # Алгоритм окрестности Мура
    count = 0
    for x, y in (
    (row - 2, col - 2), (row - 2, col - 1), (row - 2, col), (row - 2, col + 1), (row - 2, col + 2), (row - 2, col + 3),
    (row - 1, col - 2), (row - 1, col - 1), (row - 1, col), (row - 1, col + 1), (row - 1, col + 2), (row - 1, col + 3),
    (row, col - 2), (row, col - 1), (row, col), (row, col + 1), (row, col + 2), (row, col + 3),
    (row + 1, col - 2), (row + 1, col - 1), (row + 1, col), (row + 1, col + 1), (row + 1, col + 2), (row + 1, col + 3),
    (row + 2, col - 2), (row + 2, col - 1), (row + 2, col), (row + 2, col + 1), (row + 2, col + 2), (row + 2, col + 3),
    (row + 3, col - 2), (row + 3, col - 1), (row + 3, col), (row + 3, col + 1), (row + 3, col + 2), (row + 3, col + 3)):
        if not (0 <= x < len(data) and 0 <= y < len(data[x])):
            continue  # Вне границ
        if data[x][y] == 1:
            count += 1
    return count


def algorithm_moore_odd(row, col):  # Алгоритм окрестности Мура
    count = 0
    for x, y in (
    (row - 3, col - 3), (row - 3, col - 2), (row - 3, col - 1), (row - 3, col), (row - 3, col + 1), (row - 3, col + 2),
    (row - 2, col - 3), (row - 2, col - 2), (row - 2, col - 1), (row - 2, col), (row - 2, col + 1), (row - 2, col + 2),
    (row - 1, col - 3), (row - 1, col - 2), (row - 1, col - 1), (row - 1, col), (row - 1, col + 1), (row - 1, col + 2),
    (row, col - 3), (row, col - 2), (row, col - 1), (row, col), (row, col + 1), (row, col + 2),
    (row + 1, col - 3), (row + 1, col - 2), (row + 1, col - 1), (row + 1, col), (row + 1, col + 1), (row + 1, col + 2),
    (row + 2, col - 3), (row + 2, col - 2), (row + 2, col - 1), (row + 2, col), (row + 2, col + 1), (row + 2, col + 2)):
        if not (0 <= x < len(data) and 0 <= y < len(data[x])):
            continue  # Вне границ
        if data[x][y] == 1:
            count += 1
    return count


def space(n):  # Создание поля
    data = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if n // 2 - l < i < n // 2 + l and n // 2 - l < j < n // 2 + l:
                data[i][j] = 1
    return data


def set_weight(p, count):
    if count != 0:
        if p == 0:
            p += count * 2.5
            if p >= 100:
                p = 100
            w = [(100 - p) / 2, (100 - p) / 2, p]
        elif p == 50:
            p += count * 1.5
            if p >= 100:
                p = 100
            w = [(100 - p) / 2, (100 - p) / 2, p]
        elif p == 90:
            p += count * 0.4
            if p >= 100:
                p = 100
            w = [(100 - p) / 2, (100 - p) / 2, p]
    else:
        w = [(100 - p) / 2, (100 - p) / 2, p]
    return w


def main():
    global n, l, data, solubility, weight, P
    data = space(n)
    chart = pd.read_excel('Solubility Chart.xlsx', index_col=0)
    solubility = chart.loc[anion, kation]
    check_tuple = check()
    Label(up_frame, text=mes).grid(row=3, column=0, columnspan=2, sticky="we")
    P = int(check_tuple[0])
    COUNT = 0
    weight = set_weight(P, COUNT)
    epoch = 100

    plt.ion()
    fig, ax = plt.subplots()
    canvas = FigureCanvasTkAgg(fig, master=down_frame)
    canvas.get_tk_widget().pack()
    axim = ax.imshow(np.array(data), vmin=0, vmax=1, cmap='cool')

    if check_tuple[1]:
        for _ in range(epoch):
            algorithm_margolus(1)
            algorithm_margolus(-1)
            axim.set_data(np.array(data))
            canvas.draw()
            fig.canvas.flush_events()

        # plt.savefig("Solubility.png")


def btn_click():
    global kation, anion
    kation = ktxt.get()
    anion = atxt.get()
    main()


root = Tk()
root.title("Модель процесса диффузии вещества в водной среде")
root.geometry("525x725")
root.resizable(False, False)

up_frame = Frame(root)
up_frame.pack(fill="both", expand=True)
down_frame = Frame(root)
down_frame.pack(fill="both", expand=True)

klbl = Label(up_frame, text="Катион вещества:", anchor="w", width=20, height=1)
albl = Label(up_frame, text="Анион вещества:", anchor="w", width=20, height=1)
klbl.grid(row=0, column=0)
albl.grid(row=1, column=0)

ktxt = Entry(up_frame)
atxt = Entry(up_frame)
ktxt.grid(row=0, column=1, sticky="we")
atxt.grid(row=1, column=1, sticky="we")

Button(up_frame, text="Ввести", command=btn_click).grid(row=2, column=0, columnspan=2, sticky="we")

img = ImageTk.PhotoImage(Image.open("Solubility Chart.png").resize((250, 250)))
panel = Label(up_frame, image=img, width=250)
panel.grid(row=0, column=2, rowspan=4)

root.mainloop()
