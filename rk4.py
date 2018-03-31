import matplotlib.pyplot as plt
import numpy as np

hp = dict()
hp = {
    "ht": 20.0,
    "p": 50.0,
    "d": [-0.7, -0.5],
    "g": [1.0, 0.3],
    "h": 0.1,
    "t": 0.0,
    "tf": 0.7,
    "plotH": [],
    "plotP": []
}


# função para desenhar os gráficos dos hospedeiros e parasitas em função do tempo
def draw_graph(x, ylabel):
    t = np.arange(hp["t"], hp["tf"] + hp["h"], hp["h"])
    plt.title("Time(s) x " + ylabel)
    plt.plot(t, x, c='#FFCC00', lw=1, marker='o', ms=2, mec='b', mew=1)
    plt.ylabel(ylabel)
    plt.xlabel("Time(s)")
    plt.show()


# equações diferencias Lotka-Volterra
def lv_eqt(x, y, a, b):
    return x * (a + b * y)


# procedimento runge kutta 4ª ordem
def rk4_host_paras(hp):
    k = [0, 0, 0, 0]
    l = [0, 0, 0, 0]
    print(" \t\t HOST's \t\t PARASIT's")
    for i in np.arange(hp["t"], hp["tf"] + hp["h"], hp["h"]):
        print("t= %.1f \t %.5f \t\t %.5f" % (i, hp["ht"], hp["p"]))
        k[0] = hp["h"] * lv_eqt(hp["ht"], hp["p"], hp["g"][0], hp["d"][0])
        l[0] = hp["h"] * lv_eqt(hp["ht"], hp["p"], hp["d"][1], hp["g"][1])
        k[1] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[0] * hp["h"], hp["p"] + 0.5 * l[0] * hp["h"], hp["g"][0], hp["d"][0])
        l[1] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[0] * hp["h"], hp["p"] + 0.5 * l[0] * hp["h"], hp["d"][1], hp["g"][1])
        k[2] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[1] * hp["h"], hp["p"] + 0.5 * l[1] * hp["h"], hp["g"][0], hp["d"][0])
        l[2] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[1] * hp["h"], hp["p"] + 0.5 * l[1] * hp["h"], hp["d"][1], hp["g"][1])
        k[3] = hp["h"] * lv_eqt(hp["ht"] + k[2] * hp["h"], hp["p"] + l[2] * hp["h"], hp["g"][0], hp["d"][0])
        l[3] = hp["h"] * lv_eqt(hp["ht"] + k[2] * hp["h"], hp["p"] + l[2] * hp["h"], hp["d"][1], hp["g"][1])
        hp["plotH"].append(hp["ht"])
        hp["plotP"].append(hp["p"])
        hp["ht"] += (k[0] + 2.0 * (k[1] + k[2]) + k[3]) / 6.0
        hp["p"] += (l[0] + 2.0 * (l[1] + l[2]) + l[3]) / 6.0


rk4_host_paras(hp)
draw_graph(hp["plotH"], "Host's")
draw_graph(hp["plotP"], "Parasit's")