import matplotlib.pyplot as plt
import numpy as np

# dicionário contendo os parâmetros utilizados nas funções diferenciais
hp = {
    "ht": 20.,
    "p": 50.,
    "d": [-0.7, -0.5],
    "g": [1., 0.3],
    "h": 0.1,
    "t": 0.,
    "tf": 35.,
    "plotH": [],
    "plotP": []
}


# função para desenhar os gráficos dos hospedeiros e parasitas em função do tempo
def draw_graph(x, y, xlabel, ylabel):
    plt.title(xlabel + " x " + ylabel)
    plt.plot(x, y, c='#FFCC00', lw=1, marker='o', ms=2, mec='b', mew=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    # plt.plot()
    plt.savefig("%s x %s.png" % (xlabel, ylabel), dpi=(300))


# equações diferencias Lotka-Volterra
def lv_eqt(x, y, a, b):
    return x * (a + b * y)


# procedimento runge kutta 4ª ordem
def rk4_host_paras(h_p):
    f = [0, 0, 0, 0]
    g = [0, 0, 0, 0]
    print(" \t\t HOST's \t\t\t\t PARASIT's")
    for i in np.arange(h_p["t"], h_p["tf"] + h_p["h"], h_p["h"]):
        print("t= %.2f \t %.5f \t\t\t\t %.5f" % (i, h_p["ht"], h_p["p"]))

        # kn1 = f(tn, xn, yn)
        # ln1 = g(tn, yn, xn)
        # kn2 = f(tn + h / 2, xn + hkn1/2, yn + hln1/2)
        # ln2 = g(tn + h / 2, yn + hln1/2, xn + hkn1/2)
        # kn3 = f(tn + h / 2, xn + hkn2/2, yn + hln2/2)
        # ln3 = g(tn + h / 2, yn + hln2/2, xn + hkn2/2)
        # kn4 = f(tn + h, xn + hkn3, yn + hln3)
        # ln4 = g(tn + h, yn + hln3, xn + hkn3)
        # xn + 1 = xn + h(kn1 + 2kn2 + 2kn3 + kn4) / 6.0
        # yn + 1 = yn + h(ln1 + 2ln2 + 2ln3 + ln4) / 6.0

        f[0] = lv_eqt(h_p["ht"], h_p["p"], h_p["g"][0], h_p["d"][0])
        g[0] = lv_eqt(h_p["p"], h_p["ht"], h_p["d"][1], h_p["g"][1])

        f[1] = lv_eqt(h_p["ht"] + 0.5 * f[0] * h_p["h"], h_p["p"] + 0.5 * g[0] * h_p["h"], h_p["g"][0], h_p["d"][0])
        g[1] = lv_eqt(h_p["p"] + 0.5 * g[0] * h_p["h"], h_p["ht"] + 0.5 * f[0] * h_p["h"], h_p["d"][1], h_p["g"][1])

        f[2] = lv_eqt(h_p["ht"] + 0.5 * f[1] * h_p["h"], h_p["p"] + 0.5 * g[1] * h_p["h"], h_p["g"][0], h_p["d"][0])
        g[2] = lv_eqt(h_p["p"] + 0.5 * g[1] * h_p["h"], h_p["ht"] + 0.5 * f[1] * h_p["h"], h_p["d"][1], h_p["g"][1])

        f[3] = lv_eqt(h_p["ht"] + f[2] * h_p["h"], h_p["p"] + g[2] * h_p["h"], h_p["g"][0], h_p["d"][0])
        g[3] = lv_eqt(h_p["p"] + g[2] * h_p["h"], h_p["ht"] + f[2] * h_p["h"], h_p["d"][1], h_p["g"][1])

        h_p["plotH"].append(h_p["ht"])
        h_p["plotP"].append(h_p["p"])

        h_p["ht"] += (f[0] + 2.0 * (f[1] + f[2]) + f[3]) * h_p["h"] / 6.
        h_p["p"] += (g[0] + 2.0 * (g[1] + g[2]) + g[3]) * h_p["h"] / 6.


rk4_host_paras(hp)
t = np.arange(hp["t"], hp["tf"] + hp["h"], hp["h"])
draw_graph(t, hp["plotH"], "Tempo(s)", "Host's")
draw_graph(t, hp["plotP"], "Tempo(s)", "Parasit's")
draw_graph(hp["plotH"], hp["plotP"], "Host's", "Parasit's")
