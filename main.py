from random import choices
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

###
# КОНФИГУРАЦИЯ
n, l, epoch = 100, 10, 100
# ПЕРЕМЕННЫЕ
P, click, coef = 100, 0, 0
data, weight = [], []
kation, anion, mes, canvas, toolbar = "", "", "", "", ""
flag = True
###


def create_data():  # Создание поля
    arr = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if n // 2 - l < i < n // 2 + l and n // 2 - l < j < n // 2 + l:
                arr[i][j] = 1
    return arr


def check(df):  # Проверка на растворимость
    global mes, flag, coef
    p = 0

    try:
        solubility = str(df.loc[anion, kation])
    except KeyError:
        solubility = "?"

    if solubility[0].isdigit() or solubility[0] == "-":
        solubility = float(solubility)

        def fun_coef(x):
            return 0.0019 * x ** 2 + 0.1094 * x + 1.7417

        coef = fun_coef(solubility)
        if float(solubility) >= 0:
            mes = "Вещество растворимо"
            p = 80 + solubility*3.193
        elif -2.3 <= float(solubility) < 0:
            mes = "Вещество мало растворимо"
            p = 60 + solubility*8.695
        elif float(solubility) < -2.3:
            mes = "Вещество не растворимо"
            p = 20 + solubility*1.145
    else:
        flag = False
        if solubility == "x":
            mes = "Вещество не растворяется в водной среде"
        else:
            mes = "Нет сведений о существовании соединения"
    return p


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


def set_weight(p, count):
    w = []
    if count != 0:
        p -= count * coef
        if p >= 100:
            p = 100
        elif p <= 0:
            p = 1
        w = [p / 2, p / 2, 100 - p]
    else:
        w = [p / 2, p / 2, 100 - p]
    return w


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


def main():
    global l, weight, data, P, kation, anion, canvas, toolbar, click
    click += 1
    kation = ktxt.get()
    anion = atxt.get()
    data = create_data()
    df = pd.read_excel('Solubility Chart.xlsx', sheet_name="Values log", index_col=0)
    P = check(df)
    Label(up_frame, text=mes).grid(row=3, column=0, columnspan=2, sticky="we")
    weight = set_weight(P, 0)
    if flag:
        if click > 1:
            canvas.get_tk_widget().destroy()
            toolbar.destroy()
        fig = plt.figure(1)
        ax = fig.add_subplot()
        canvas = FigureCanvasTkAgg(fig, master=down_frame)
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas, down_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack()

        plt.ion()
        axim = ax.imshow(np.array(data), vmin=0, vmax=1, cmap='cool')
        counter = 0
        Label(up_frame, text=f"Количество итераций: {counter}").grid(row=4, column=0, columnspan=2, sticky="we")
        for _ in range(epoch):
            algorithm_margolus(1)
            algorithm_margolus(-1)
            axim.set_data(np.array(data))
            fig.canvas.flush_events()
            counter += 1
            Label(up_frame, text=f"Количество итераций: {counter}").grid(row=4, column=0, columnspan=2, sticky="we")
        # plt.savefig("Solubility new H.png")


root = Tk()
root.title("Моделирование процесса диффузии вещества в водной среде")
root.geometry("735x875")

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

Button(up_frame, text="Ввести", command=main).grid(row=2, column=0, columnspan=2, sticky="we")
root.bind_all('<Return>', lambda event: main())

img = ImageTk.PhotoImage(Image.open("Solubility Chart.png").resize((460, 350)))
panel = Label(up_frame, image=img, width=460)
panel.grid(row=0, column=2, rowspan=5)

root.mainloop()
