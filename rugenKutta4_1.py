import matplotlib.pyplot as plt
import numpy as np


# função para desenhar os gráficos dos hospedeiros e parasitas em função do tempo
def draw_graph(x, y, xlabel, ylabel):
    plt.title(xlabel + " x " + ylabel)
    plt.plot(x, y, c='#FFCC00', lw=1, marker='+', ms=2, mec='b', mew=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale('symlog')
    # plt.yscale('symlog')
    plt.savefig("%s x %s.png" % (xlabel, ylabel), dpi=(500))
    plt.show()


# equações diferencias Lotka-Volterra
# OBS: Foi utilizada uma única função para representação das duas funções diferencias
# tanto hospedeiros como parasitas podem ser calculados alternando-se os parâmetros
def lv_eqt(x, y, a, b):
    return x * (a + b * y)


# procedimento runge kutta 4ª ordem
def rk4_host_paras(h_p):
    f = [0, 0, 0, 0]
    g = [0, 0, 0, 0]
    print("\t\t\t HOST's \t\t PARASIT's")
    for i in np.arange(h_p["t"], h_p["tf"] + h_p["h"], h_p["h"]):
        print("t= %.1f \t\t %.5f \t\t %.5f" % (i, h_p["ht"], h_p["p"]))

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


def main():
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
    print(
        "\n\nPrograma para cálculo do número de Hospedeiros e Parasitas através da correlação existente entre os mesmos"
        "a partir das equações diferenciais que seguem: ")
    print("dH/dt = g1H - d1HP\ndP/dt = -d2P + g2HP\n")
    print("Para tal solução foi utilizado o método numérico Rugen-Kutta de ordem 4.")
    print("\nParâmetros padrões:\n")
    print("Hospedeiros: 20.0\nParasitas: 50.0\nd1: 0.7\nd2: 0.5\ng1: 1.0\ng2: 0.3\nh: 0.1\nTi: 0.0\nTf: 35.0")
    print("\n\n")
    ch = (int)(input("Gostaria de alterar os parâmetros? ('1' para sim, ou '0' para não)"))
    if ch:
        hp["ht"] = (float)(input("Entre com o número de Hospedeiros."))
        hp["p"] = (float)(input("Entre com o número de Parasitas."))
        hp["d"][0] = -(float)(input("Entre com o parâmetro 'd1' (insira seu valor absoluto)."))
        hp["d"][1] = -(float)(input("Entre com o parâmetro 'd2' (insira seu valor absoluto)."))
        hp["g"][0] = (float)(input("Entre com o parâmetro 'g1'."))
        hp["g"][1] = (float)(input("Entre com o parâmetro 'g2'."))
        hp["h"] = (float)(input("Entre com o parâmetro 'h'."))
        hp["t"] = (float)(input("Entre com o parâmetro 'ti'."))
        hp["tf"] = (float)(input("Entre com o parâmetro 'tf'."))
    print("\n\n")
    rk4_host_paras(hp)
    t = np.arange(hp["t"], hp["tf"] + hp["h"], hp["h"])
    draw_graph(t, hp["plotH"], "Tempo(s)", "Host's")
    draw_graph(t, hp["plotP"], "Tempo(s)", "Parasit's")
    draw_graph(hp["plotH"], hp["plotP"], "Host's", "Parasit's")


main()
