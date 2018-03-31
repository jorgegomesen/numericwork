import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
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
    "plotH":[],
    "plotP":[]
}


def draw_graph(x, ylabel):
    ax = host_subplot(111, axes_class=AA.Axes)
    t = np.arange(hp["t"], hp["tf"] + hp["h"], hp["h"])
    ax.plot(t, x)

    ax2 = ax.twin()  # ax2 is responsible for "top" axis and "right" axis
    ax2.set_xticks([0., .5 * np.pi, np.pi, 1.5 * np.pi, 2 * np.pi])
    ax2.set_xticklabels(["$0$", r"$\frac{1}{2}\pi$",
                         r"$\pi$", r"$\frac{3}{2}\pi$", r"$2\pi$"])

    ax2.axis["right"].major_ticklabels.set_visible(False)
    ax2.axis["top"].major_ticklabels.set_visible(False)

    plt.xlabel("Time(s)")
    plt.ylabel(ylabel)
    plt.draw()
    plt.show()


# equações diferencias Lotka-Volterra
def lv_eqt(x, y, a, b):
    return x * (a + b * y)


# procedimento runge kutta 4ª ordem
def rk4_host_paras(hp):
    k = [0, 0, 0, 0]
    l = [0, 0, 0, 0]
    print(" \t\t HOST's \t\t PARASIT's")
    for i in np.arange(hp["t"], hp["tf"]+hp["h"], hp["h"]):
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
print(hp["plotH"])
print(hp["plotP"])
draw_graph(hp["plotH"], "Host's")
draw_graph(hp["plotP"], "Parasit's")
