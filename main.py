from random import choices
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tkinter import *
from PIL import ImageTk, Image
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

###
# КОНСТАНТЫ
SIZE, EPOCH = 150, 150
# ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ
FIELD = []
FIELD_OF_CONCENTRATION = []
START = True
AMOUNT_OF_SUBSTANCE = 0
PROBABILITY = 100
COEFFICIENT = 0
click1, click2 = 0, 0
mes, canvas1, toolbar1, canvas2, toolbar2 = "", "", "", "", ""
###


def create_field():
    """Создание водной среды"""
    global AMOUNT_OF_SUBSTANCE, FIELD, FIELD_OF_CONCENTRATION
    r = SIZE // 10
    FIELD = [[0 for _ in range(SIZE)] for _ in range(SIZE)]
    FIELD_OF_CONCENTRATION = np.zeros((SIZE, SIZE))
    for i in range(SIZE):
        for j in range(SIZE):
            p2 = (i - SIZE // 2) ** 2 + (j - SIZE // 2) ** 2
            if p2 < r**2:
                FIELD[i][j] = 1
                AMOUNT_OF_SUBSTANCE += 1


def check(df, kation, anion):
    """Проверка на растворимость"""
    global mes, START, PROBABILITY, COEFFICIENT
    try:
        solubility = str(df.loc[anion, kation])
    except KeyError:
        solubility = "?"

    if solubility[0].isdigit() or solubility[0] == "-":
        solubility = float(solubility)
        COEFFICIENT = (lambda x: 0.0019 * x ** 2 + 0.1094 * x + 1.7417)(solubility)
        if solubility >= 0:
            mes = "Вещество растворимо"
            PROBABILITY = 80 + solubility*3.193
        elif -2.3 <= solubility < 0:
            mes = "Вещество мало растворимо"
            PROBABILITY = 60 + solubility*8.695
        elif solubility < -2.3:
            mes = "Вещество не растворимо"
            PROBABILITY = 20 + solubility*1.145
    else:
        START = False
        if solubility == "x":
            mes = "Вещество не растворяется в водной среде"
        else:
            mes = "Нет сведений о существовании соединения"


def rotation(temp, w):  # Поворот
    """Поворот блока по и против часовой стрелки"""
    def randomizer(p):  # Рандомайзер
        ch = choices((True, False, "Stop"), weights=p, k=1)
        return ch[0]

    r = randomizer(w)
    if r is True:  # По часовой
        rotated = list(zip(*temp[::-1]))
    elif r is False:  # Против часовой
        rotated = list(zip(*temp))[::-1]
    else:
        rotated = temp
    return rotated


def set_weight(count, p):
    """Установить весовые параметры для блока"""
    if count != 0:
        p -= count * COEFFICIENT
        if p >= 100:
            p = 100
        elif p <= 0:
            p = 1
        w = [p / 2, p / 2, 100 - p]
    else:
        w = [p / 2, p / 2, 100 - p]
    return w


def algorithm_margolus(k):
    """Алгоритм окрестности Марголуса"""
    for row in range(1, SIZE, 2):
        for col in range(1, SIZE, 2):
            temp = []
            counted = algorithm_moore(row, col, k)
            weight = set_weight(counted, PROBABILITY)
            if (not row == SIZE - 1 and not col == SIZE - 1) and k == 1:
                temp.append([FIELD[row][col], FIELD[row][col + k]])
                temp.append([FIELD[row + k][col], FIELD[row + k][col + k]])
                rotated = rotation(temp, weight)
                FIELD[row][col], FIELD[row][col + k], FIELD[row + k][col], FIELD[row + k][col + k] = rotated[0][0], rotated[0][1], rotated[1][0], rotated[1][1]
            elif k == -1:
                temp.append([FIELD[row][col], FIELD[row][col + k]])
                temp.append([FIELD[row + k][col], FIELD[row + k][col + k]])
                rotated = rotation(temp, weight)
                FIELD[row][col], FIELD[row][col + k], FIELD[row + k][col], FIELD[row + k][col + k] = rotated[-1][-1], rotated[-1][0], rotated[0][-1], rotated[0][0]


def algorithm_moore(i, j, k):
    """Алгоритм окрестности Мура"""
    count = 0
    if k == 1:
        for x, y in (
                (i - 2, j - 2), (i - 2, j - 1), (i - 2, j), (i - 2, j + 1), (i - 2, j + 2), (i - 2, j + 3),
                (i - 1, j - 2), (i - 1, j - 1), (i - 1, j), (i - 1, j + 1), (i - 1, j + 2), (i - 1, j + 3),
                (i, j - 2), (i, j - 1), (i, j), (i, j + 1), (i, j + 2), (i, j + 3),
                (i + 1, j - 2), (i + 1, j - 1), (i + 1, j), (i + 1, j + 1), (i + 1, j + 2), (i + 1, j + 3),
                (i + 2, j - 2), (i + 2, j - 1), (i + 2, j), (i + 2, j + 1), (i + 2, j + 2), (i + 2, j + 3),
                (i + 3, j - 2), (i + 3, j - 1), (i + 3, j), (i + 3, j + 1), (i + 3, j + 2), (i + 3, j + 3)):
            if not (0 <= x < SIZE and 0 <= y < SIZE):
                continue  # Вне границ
            if FIELD[x][y] == 1:
                count += 1
    elif k == -1:
        for x, y in (
                (i - 3, j - 3), (i - 3, j - 2), (i - 3, j - 1), (i - 3, j), (i - 3, j + 1), (i - 3, j + 2),
                (i - 2, j - 3), (i - 2, j - 2), (i - 2, j - 1), (i - 2, j), (i - 2, j + 1), (i - 2, j + 2),
                (i - 1, j - 3), (i - 1, j - 2), (i - 1, j - 1), (i - 1, j), (i - 1, j + 1), (i - 1, j + 2),
                (i, j - 3), (i, j - 2), (i, j - 1), (i, j), (i, j + 1), (i, j + 2),
                (i + 1, j - 3), (i + 1, j - 2), (i + 1, j - 1), (i + 1, j), (i + 1, j + 1), (i + 1, j + 2),
                (i + 2, j - 3), (i + 2, j - 2), (i + 2, j - 1), (i + 2, j), (i + 2, j + 1), (i + 2, j + 2)):
            if not (0 <= x < SIZE and 0 <= y < SIZE):
                continue  # Вне границ
            if FIELD[x][y] == 1:
                count += 1
    return count


def get_conc(i, j):
    """Рассчитать количество вещества в точке, окрестностью 11х11"""
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
        if not (0 <= x < SIZE and 0 <= y < SIZE):
            continue  # Вне границ
        if FIELD[x][y] == 1:
            count += 1
    return count


def onclick(event):
    """Ивент клика мышкой по графику"""
    global slbl
    ix, iy = int(event.xdata), int(event.ydata)
    slbl = Label(head_frame, text=f"Концентрация в точке {ix, iy} равна {round(FIELD_OF_CONCENTRATION[iy][ix] / AMOUNT_OF_SUBSTANCE, 3)}")
    slbl.grid(row=6, column=0, columnspan=2, sticky="we", pady=20)


def btn_click():
    """Ивент клика по кнопке запуска"""
    global canvas2, toolbar2, click2
    click2 += 1
    if click2 > 1:
        canvas2.get_tk_widget().destroy()
        toolbar2.destroy()

    ###
    # График
    fig2 = plt.figure()
    ax2 = fig2.add_subplot()
    fig2.patch.set_facecolor('#F0F0F0')
    plt.xlim([0, SIZE])
    plt.ylim([0, SIZE])
    canvas2 = FigureCanvasTkAgg(fig2, master=fig_frame2)
    canvas2.draw()
    canvas2.get_tk_widget().pack()
    toolbar2 = NavigationToolbar2Tk(canvas2, fig_frame2, pack_toolbar=False)
    toolbar2.update()
    toolbar2.pack()
    axim2 = ax2.pcolormesh(FIELD_OF_CONCENTRATION, vmin=0, vmax=120, cmap='cool')
    fig2.colorbar(axim2, ax=ax2)
    axim2.set_mouseover(True)
    plt.close()
    ###

    canvas2.mpl_connect("button_press_event", onclick)


def main():
    """Приложение"""
    global canvas1, toolbar1, click1, slbl, ilbl
    click1 += 1

    create_field()
    df = pd.read_excel('Solubility Chart.xlsx', sheet_name="Values log", index_col=0)
    check(df, ktxt.get(), atxt.get())

    slbl = Label(head_frame, text=mes)
    slbl.grid(row=4, column=0, columnspan=2, sticky="we", pady=20)

    if START:  # Если вещество растворимо, начинается процесс
        if click1 > 1:
            canvas1.get_tk_widget().destroy()
            toolbar1.destroy()

        ###
        # График
        fig = plt.figure(1)
        ax = fig.add_subplot()
        fig.patch.set_facecolor('#F0F0F0')
        plt.xlim([0, SIZE])
        plt.ylim([0, SIZE])
        canvas1 = FigureCanvasTkAgg(fig, master=fig_frame1)
        canvas1.draw()
        canvas1.get_tk_widget().pack()
        toolbar1 = NavigationToolbar2Tk(canvas1, fig_frame1, pack_toolbar=False)
        toolbar1.update()
        toolbar1.pack()
        plt.ion()
        axim = ax.imshow(np.array(FIELD), cmap='cool')
        fig.colorbar(axim, ax=ax)
        counter = 0
        ilbl = Label(head_frame, text=f"Количество итераций: {counter}")
        ilbl.grid(row=5, column=0, columnspan=2, sticky="we", pady=20)
        ###

        for _ in range(EPOCH):
            algorithm_margolus(1)
            algorithm_margolus(-1)
            axim.set_data(np.array(FIELD))
            fig.canvas.flush_events()
            counter += 1
            ilbl = Label(head_frame, text=f"Количество итераций: {counter}")
            ilbl.grid(row=5, column=0, columnspan=2, sticky="we", pady=20)

        for i in range(len(FIELD_OF_CONCENTRATION)):
            for j in range(len(FIELD_OF_CONCENTRATION[i])):
                FIELD_OF_CONCENTRATION[i][j] = get_conc(i, j)


# Запуск приложения
root = Tk()
root.title("Моделирование процесса диффузии вещества в водной среде")
root.geometry("1110x1045")
root.iconbitmap("Icon.ico")

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
