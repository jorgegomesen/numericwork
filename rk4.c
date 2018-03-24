#include <stdio.h>

typedef struct hosp_Paras{
  float g[2];
  float d[2];
  float h;
  float t;
  float tf;
  double Par;
  double Hosp;
  double k[4];
  double l[4];
}HOSP_PARAS;

double rungeFunct(double, double, float, float);
void HospParas(HOSP_PARAS *);

double rungeFunct(double x, double y, float a, float b){
    return ((a * x) + (b * x * y));
}

void HospParasRungeKutta(HOSP_PARAS *hp){
    double Hosp = hp->Hosp, Paras = hp->Par;
    float t = hp->t, tf = hp->tf, h=hp->h;
    for(; t <= tf; t += h){
        printf("\nt-> %.1f -- Hospedeiros -> %.8lf ---- Pararasitas -> %.8lf ", t, Hosp, Paras);
        hp->k[0] = h * rungeFunct(Hosp, Paras, hp->g[0], -(hp->d[0]));
        hp->l[0] = h * rungeFunct(Paras, Hosp, -(hp->d[1]), hp->g[1]);
        hp->k[1] = h * rungeFunct((Hosp + 0.5 * h * hp->k[0]), (Paras + 0.5 * h * hp->l[0]), hp->g[0], -(hp->d[0]));
        hp->l[1] = h * rungeFunct((Paras + 0.5 * h * hp->l[0]), (Hosp + 0.5 * h * hp->k[0]), -(hp->d[1]), hp->g[1]);
        hp->k[2] = h * rungeFunct((Hosp + 0.5 * h * hp->k[1]), (Paras + 0.5 * h * hp->l[1]), hp->g[0], -(hp->d[0]));
        hp->l[2] = h * rungeFunct((Paras + 0.5 * h * hp->l[1]), (Hosp + 0.5 * h * hp->k[1]), -(hp->d[1]), hp->g[1]);
        hp->k[3] = h * rungeFunct((Hosp + hp->k[2]), (Paras + hp->l[2]), hp->g[0], -(hp->d[0]));
        hp->l[3] = h * rungeFunct((Paras + hp->l[2]), (Hosp + hp->k[2]), -(hp->d[1]), hp->g[1]);
        Hosp += ((1 / 6.0) * (hp->k[0] + (2.0 * hp->k[1]) + (2.0 * hp->k[2]) + hp->k[3]));
        Paras += ((1 / 6.0) * (hp->l[0] + (2.0 * hp->l[1]) + (2.0 * hp->l[2]) + hp->l[3]));
    }
    hp->Hosp = Hosp;
    hp->Par = Paras;
    printf("\n\n");
}

int main(int argc, char** argv){
    HOSP_PARAS hp;
    hp.g[0] = 1.0;
    hp.g[1] = 0.3;
    hp.d[0] = 0.7;
    hp.d[1] = 0.5;
    hp.h = 0.1;
    hp.t = 0.0;
    hp.tf = 35.0;
    hp.Par = 50.0;
    hp.Hosp = 20.0;
    HospParasRungeKutta(&hp);
    return 0;
}