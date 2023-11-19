from random import choices
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

###
# КОНФИГУРАЦИЯ
n, l, epoch, q = 100, 10, 100, 0
# ПЕРЕМЕННЫЕ
P, click, click2, coef = 100, 0, 0, 0
data, data_conc, weight = [], [], []
kation, anion, mes, canvas, toolbar, canvas2, toolbar2 = "", "", "", "", "", "", ""
flag = True
###


def create_data():  # Создание поля
    global q
    arr = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if n // 2 - l < i < n // 2 + l and n // 2 - l < j < n // 2 + l:
                arr[i][j] = 1
                q += 1
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


def get_conc(i, j):
    count = 0
    for x, y in (
            (i - 5, j - 5), (i - 5, j - 4), (i - 5, j - 3), (i - 5, j - 2), (i - 5, j - 1), (i - 5, j), (i - 5, j + 1),
            (i - 5, j + 2), (i - 5, j + 3), (i - 5, j + 4), (i - 5, j + 5),
            (i - 4, j - 5), (i - 4, j - 4), (i - 4, j - 3), (i - 4, j - 2), (i - 4, j - 1), (i - 4, j), (i - 4, j + 1),
            (i - 4, j + 2), (i - 4, j + 3), (i - 4, j + 4), (i - 4, j + 5),
            (i - 3, j - 5), (i - 3, j - 4), (i - 3, j - 3), (i - 3, j - 2), (i - 3, j - 1), (i - 3, j), (i - 3, j + 1),
            (i - 3, j + 2), (i - 3, j + 3), (i - 3, j + 4), (i - 3, j + 5),
            (i - 2, j - 5), (i - 2, j - 4), (i - 2, j - 3), (i - 2, j - 2), (i - 2, j - 1), (i - 2, j), (i - 2, j + 1),
            (i - 2, j + 2), (i - 2, j + 3), (i - 2, j + 4), (i - 2, j + 5),
            (i - 1, j - 5), (i - 1, j - 4), (i - 1, j - 3), (i - 1, j - 2), (i - 1, j - 1), (i - 1, j), (i - 1, j + 1),
            (i - 1, j + 2), (i - 1, j + 3), (i - 1, j + 4), (i - 1, j + 5),
            (i, j - 5), (i, j - 4), (i, j - 3), (i, j - 2), (i, j - 1), (i, j), (i, j + 1), (i, j + 2), (i, j + 3),
            (i, j + 4), (i, j + 5),
            (i + 1, j - 5), (i + 1, j - 4), (i + 1, j - 3), (i + 1, j - 2), (i + 1, j - 1), (i + 1, j), (i + 1, j + 1),
            (i + 1, j + 2), (i + 1, j + 3), (i + 1, j + 4), (i + 1, j + 5),
            (i + 2, j - 5), (i + 2, j - 4), (i + 2, j - 3), (i + 2, j - 2), (i + 2, j - 1), (i + 2, j), (i + 2, j + 1),
            (i + 2, j + 2), (i + 2, j + 3), (i + 2, j + 4), (i + 2, j + 5),
            (i + 3, j - 5), (i + 3, j - 4), (i + 3, j - 3), (i + 3, j - 2), (i + 3, j - 1), (i + 3, j), (i + 3, j + 1),
            (i + 3, j + 2), (i + 3, j + 3), (i + 3, j + 4), (i + 3, j + 5),
            (i + 4, j - 5), (i + 4, j - 4), (i + 4, j - 3), (i + 4, j - 2), (i + 4, j - 1), (i + 4, j), (i + 4, j + 1),
            (i + 4, j + 2), (i + 4, j + 3), (i + 4, j + 4), (i + 4, j + 5),
            (i + 5, j - 5), (i + 5, j - 4), (i + 5, j - 3), (i + 5, j - 2), (i + 5, j - 1), (i + 5, j), (i + 5, j + 1),
            (i + 5, j + 2), (i + 5, j + 3), (i + 5, j + 4), (i + 5, j + 5)):
        if not (0 <= x < len(data_conc) and 0 <= y < len(data[x])):
            continue  # Вне границ
        if data[x][y] == 1:
            count += 1
    return count


def onclick(event):
    global slbl
    ix, iy = int(event.xdata), int(event.ydata)
    slbl = Label(head_frame, text=f"Концентрация в точке {ix, iy} равна {round(data_conc[iy][ix]/q, 3)}")
    slbl.grid(row=6, column=0, columnspan=2, sticky="we", pady=20)


def btn_click():
    global canvas2, toolbar2, click2
    click2 += 1
    if click2 > 1:
        canvas2.get_tk_widget().destroy()
        toolbar2.destroy()
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    canvas2 = FigureCanvasTkAgg(fig2, master=fig_frame2)
    canvas2.draw()
    canvas2.get_tk_widget().pack()
    toolbar2 = NavigationToolbar2Tk(canvas2, fig_frame2, pack_toolbar=False)
    toolbar2.update()
    toolbar2.pack()
    #axim2 = ax2.imshow(data_conc, vmin=0, vmax=120, cmap='cool')
    #axim2.set_data(data_conc)
    axim2 = ax2.pcolormesh(data_conc, vmin=0, vmax=120, cmap='cool')
    fig2.colorbar(axim2, ax=ax2)
    axim2.set_mouseover(True)
    plt.close()

    canvas2.mpl_connect("button_press_event", onclick)


def main():
    global l, weight, data, data_conc, P, kation, anion, canvas, toolbar, click, slbl, ilbl
    click += 1
    kation = ktxt.get()
    anion = atxt.get()
    data = create_data()
    df = pd.read_excel('Solubility Chart.xlsx', sheet_name="Values log", index_col=0)
    P = check(df)
    slbl = Label(head_frame, text=mes)
    slbl.grid(row=4, column=0, columnspan=2, sticky="we", pady=20)
    weight = set_weight(P, 0)
    if flag:
        if click > 1:
            canvas.get_tk_widget().destroy()
            toolbar.destroy()

        fig = plt.figure(1)
        ax = fig.add_subplot()
        canvas = FigureCanvasTkAgg(fig, master=fig_frame1)
        canvas.draw()
        canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(canvas, fig_frame1, pack_toolbar=False)
        toolbar.update()
        toolbar.pack()

        plt.ion()
        axim = ax.imshow(np.array(data), cmap='cool')
        fig.colorbar(axim, ax=ax)
        counter = 0
        ilbl = Label(head_frame, text=f"Количество итераций: {counter}")
        ilbl.grid(row=5, column=0, columnspan=2, sticky="we", pady=20)
        for _ in range(epoch):
            algorithm_margolus(1)
            algorithm_margolus(-1)
            axim.set_data(np.array(data))
            fig.canvas.flush_events()
            counter += 1
            ilbl = Label(head_frame, text=f"Количество итераций: {counter}")
            ilbl.grid(row=5, column=0, columnspan=2, sticky="we", pady=20)

        data_conc = [[0 for _ in range(n)] for _ in range(n)]
        for i in range(len(data_conc)):
            for j in range(len(data_conc[i])):
                data_conc[i][j] = get_conc(i, j)
        data_conc = np.array(data_conc)


root = Tk()
root.title("Моделирование процесса диффузии вещества в водной среде")
root.geometry("1110x1045")

head_frame = Frame(root, width=640, height=522)
head_frame.grid(row=0, column=0)
table_frame = Frame(root, width=640, height=522)
table_frame.grid(row=1, column=0)
fig_frame1 = Frame(root, width=640, height=522)
fig_frame1.grid(row=0, column=1)
fig_frame2 = Frame(root, width=640, height=522)
fig_frame2.grid(row=1, column=1)

klbl = Label(head_frame, text="Катион вещества:", anchor="w", width=20, height=1)
albl = Label(head_frame, text="Анион вещества:", anchor="w", width=20, height=1)
klbl.grid(row=0, column=0, pady=20)
albl.grid(row=1, column=0, pady=20)

ktxt = Entry(head_frame)
atxt = Entry(head_frame)
ktxt.insert(0, "Na")
atxt.insert(0, "OH")
ktxt.grid(row=0, column=1, sticky="we", pady=20)
atxt.grid(row=1, column=1, sticky="we", pady=20)

Button(head_frame, text="Ввести", command=main).grid(row=2, column=0, columnspan=2, sticky="we", pady=20)
root.bind_all('<Return>', lambda event: main())

Button(head_frame, text="Построить градиент", command=btn_click).grid(row=3, column=0, columnspan=2, sticky="we", pady=20)

slbl = Label(head_frame, text="Введите вещество")
slbl.grid(row=4, column=0, columnspan=2, sticky="we", pady=20)

ilbl = Label(head_frame, text="Количество итераций: 0")
ilbl.grid(row=5, column=0, columnspan=2, sticky="we", pady=20)

slbl = Label(head_frame, text="Выберите точку")
slbl.grid(row=6, column=0, columnspan=2, sticky="we", pady=20)

img = ImageTk.PhotoImage(Image.open("Solubility Chart.png").resize((460, 350)))
table = Label(table_frame, image=img)
table.pack()

root.mainloop()
