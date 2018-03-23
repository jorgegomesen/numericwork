#include <stdio.h>
#include <math.h>

// float rungeFunct(float, float);
// float rungeKutta4(float, float, float, float);

float rungeFunct(float y, float d, float g, float Par, int param){
    if (!param){
        // printf("\ny -> %.4f -- (%.4f * %.4f) - (%.4f * %.4f * %.4f) = %.4f", y, g, y, d, Par, y, (g * y) - (d * Par * y));
        return ((g * y) - (d * Par * y));
    }else{
        // printf("\ny -> %.4f -- -(%.4f * %.4f) + (%.4f * %.4f * %.4f) = %.4f", y, d, y, g, Par, y, -(d * y) + (g * Par * y));
        return (-(d * y) + (g * Par * y));
    }
}

float rungeKutta4(float y, float h, float d, float g, float Par, int param){
    float k1, k2, k3, k4;
    k1 = rungeFunct(y, d, g, Par, param);
    k2 = rungeFunct(y + (h/2) * k1, d, g, Par, param);
    k3 = rungeFunct(y + (h/2) * k2, d, g, Par, param);
    k4 = rungeFunct(y + h * k3, d, g, Par, param);
    y += (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
    return y;
}

float HospParas(float d1, float d2, float g1, float g2, float Hosp, float Par, float t, float tf, float h){
    float aux;
    for(; t <= 1; t+=h){
        aux = rungeKutta4(Hosp, h, d1, g1, Par, 0);
        Par = rungeKutta4(Par, h, d2, g2, Hosp, 1);
        printf("\nHospedeiros -> %.4f ---- Pararasitas -> %.4f ", aux, Par);
        Hosp = aux;
    }
}

int main(int argc, char** argv){
    float Par = 50, Hosp = 20, g1 = 1, g2 = 0.3, d1 = 0.7, d2 = 0.5, h = 0.1, t = 0, tf = 35;
    // printf("\n\nEntre com os Parasitas. ");
    // scanf("%f", &Par);
    // printf("\nEntre com os Hospedeiros. ");
    // scanf("%f", &Hosp);
    // printf("\nEntre com o g1. ");
    // scanf("%f", &g1);
    // printf("\nEntre com o g2. ");
    // scanf("%f", &g2);
    // printf("\nEntre com o d1. ");
    // scanf("%f", &d1);
    // printf("\nEntre com o d2. ");
    // scanf("%f", &d2);
    // printf("\nEntre com o h. ");
    // scanf("%f", &h);
    // printf("\nEntre com o Tempo inicial. ");
    // scanf("%f", &t);
    // printf("\nEntre com o Tempo final. ");
    // scanf("%f", &tf);
    printf("\n\nResultado: %.4f\n\n", HospParas(d1, d2, g1, g2, Hosp, Par, t, tf, h));
    return 0;
}