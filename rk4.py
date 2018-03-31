hp = dict()
hp = {"ht": 2.0, "p": 1.0, "d": [-1.0, -1.0], "g": [1.0, 1.0], "h": 0.1, "t": 0.0, "tf": 3.3}


# equações diferencias Lotka-Volterra
def lv_eqt(x, y, a, b):
    return x * (a + b * y)


# procedimento runge kutta 4ª ordem
def rk4_host_paras(hp):
    k = [0, 0, 0, 0]
    l = [0, 0, 0, 0]
    print("HOST's \t\t PARASIT's")
    for i in range(0, 50, 1):
        print("%.5f \t %.5f" % (hp["ht"], hp["p"]))
        k[0] = hp["h"] * lv_eqt(hp["ht"], hp["p"], hp["g"][0], hp["d"][0])
        l[0] = hp["h"] * lv_eqt(hp["ht"], hp["p"], hp["d"][1], hp["g"][1])

        k[1] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[0] * hp["h"], hp["p"] + 0.5 * l[0] * hp["h"], hp["g"][0], hp["d"][0])
        l[1] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[0] * hp["h"], hp["p"] + 0.5 * l[0] * hp["h"], hp["d"][1], hp["g"][1])

        k[2] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[1] * hp["h"], hp["p"] + 0.5 * l[1] * hp["h"], hp["g"][0], hp["d"][0])
        l[2] = hp["h"] * lv_eqt(hp["ht"] + 0.5 * k[1] * hp["h"], hp["p"] + 0.5 * l[1] * hp["h"], hp["d"][1], hp["g"][1])

        k[3] = hp["h"] * lv_eqt(hp["ht"] + k[2] * hp["h"], hp["p"] + l[2] * hp["h"], hp["g"][0], hp["d"][0])
        l[3] = hp["h"] * lv_eqt(hp["ht"] + k[2] * hp["h"], hp["p"] + l[2] * hp["h"], hp["d"][1], hp["g"][1])

        hp["ht"] += (k[0] + 2.0 * (k[1] + k[2]) + k[3]) / 6.0
        hp["p"] += (l[0] + 2.0 * (l[1] + l[2]) + l[3]) / 6.0


rk4_host_paras(hp)
