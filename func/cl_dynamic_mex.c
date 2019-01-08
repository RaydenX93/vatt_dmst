/* ========================================================
 * Funzione "cl_dynamic" pronta per uso MEX in MATLAB
 * Versione funzione C: 25/09/2018 (con leggere modifiche)
 * Versione funzione MEX: 04/10/2018
 * ========================================================
 * - Sintassi per compilare:
 *   mex cl_dynamic_mex.c
 * - Sintassi per chiamare binario da MATLAB:
 *   cl_dyn = cl_dynamic_mex(chord, omega, Uinf, modW, Re_c, alpha_in, alpha_in_0, alpha_max, alpha_min)
 *   (Tutte le variabili sono scalari) */

#include <math.h>
#include "mex.h"

/* Definizioni per controllo input */
#define IS_double_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_double_SCALAR(P) (IS_double_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)

/* ============== DICHIARAZIONE VARIABILI GLOBALI ============== */
double attack[181]={0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00, 19.00, 20.00, 21.00, 22.00, 23.00, 24.00, 25.00, 26.00, 27.00, 28.00, 29.00, 30.00, 31.00, 32.00, 33.00, 34.000, 35.00, 36.00, 37.00, 38.00, 39.00, 40.00, 41.00, 42.00, 43.00, 44.00, 45.00, 46.00, 47.00, 48.00, 49.00, 50.00, 51.00, 52.00, 53.00, 54.00, 55.00, 56.00, 57.00, 58.00, 59.00, 60.00, 61.00, 62.00, 63.00, 64.00, 65.00, 66.00, 67.00, 68.00, 69.00, 70.00, 71.00, 72.00, 73.00, 74.00, 75.00, 76.00, 77.00, 78.00, 79.00, 80.00, 81.00, 82.00, 83.00, 84.00, 85.00, 86.00, 87.00, 88.00, 89.00, 90.00, 91.00, 92.00, 93.00, 94.00, 95.00, 96.00, 97.00, 98.00, 99.00, 100.00, 101.00, 102.00, 103.00, 104.00, 105.00, 106.00, 107.00, 108.00, 109.00, 110.00, 111.00, 112.00, 113.00, 114.00, 115.00, 116.00, 117.00, 118.00, 119.00, 120.00, 121.00, 122.00, 123.00, 124.00, 125.00, 126.00, 127.00, 128.00, 129.00, 130.00, 131.00, 132.00, 133.00, 134.000, 135.00, 136.00, 137.00, 138.00, 139.00, 140.00, 141.00, 142.00, 143.00, 144.00, 145.00, 146.00, 147.00, 148.00, 149.00, 150.00, 151.00, 152.00, 153.00, 154.00, 155.00, 156.00, 157.00, 158.00, 159.00, 160.00, 161.00, 162.00, 163.00, 164.00, 165.00, 166.00, 167.00, 168.00, 169.00, 170.00, 171.00, 172.00, 173.00, 174.00, 175.00, 176.00, 177.00, 178.00, 179.00, 180.00};
double Re_levels[5]={7e5, 1e6, 2e6, 5e6, 1e7};
double cl_s[5][181]={{0,0.11,0.22,0.33,0.44,0.55,0.66,0.7483,0.8442,0.926,0.9937,1.0363,1.0508,1.0302,0.9801,0.9119,0.8401,0.7799,0.7305,0.7041,0.699,0.7097,0.7298,0.7593,0.7961,0.8353,0.8838,0.9473,0.92337,0.87893,0.855,0.86556,0.89167,0.925,0.95722,0.98,0.99474,1.0081,1.0196,1.0287,1.035,1.0396,1.0437,1.047,1.0492,1.05,1.0482,1.0434,1.0365,1.0284,1.02,1.0103,0.99809,0.98427,0.96963,0.955,0.94052,0.92557,0.90986,0.8931,0.875,0.85487,0.83265,0.809,0.78456,0.76,0.73502,0.70915,0.68276,0.65625,0.63,0.60408,0.57824,0.55235,0.52631,0.5,0.47333,0.44637,0.41924,0.39208,0.365,0.33808,0.31124,0.28435,0.25731,0.23,0.20233,0.17437,0.14624,0.11808,0.09,0.061919,0.033756,0.0056335,-0.022326,-0.05,-0.077314,-0.10435,-0.13124,-0.15808,-0.185,-0.21208,-0.23924,-0.26637,-0.29333,-0.32,-0.3464,-0.3726,-0.3986,-0.4244,-0.45,-0.47586,-0.50199,-0.52769,-0.55226,-0.575,-0.59574,-0.61511,-0.63361,-0.65174,-0.67,-0.68831,-0.70635,-0.72423,-0.74208,-0.76,-0.77817,-0.79651,-0.81476,-0.83268,-0.85,-0.86719,-0.88445,-0.90111,-0.91651,-0.93,-0.94308,-0.95646,-0.96831,-0.97677,-0.98,-0.97485,-0.96135,-0.94242,-0.921,-0.9,-0.87742,-0.85083,-0.82253,-0.79482,-0.77,-0.74679,-0.7235,-0.70181,-0.68342,-0.67,-0.65972,-0.65021,-0.64234,-0.63698,-0.635,-0.6374,-0.64401,-0.65391,-0.66621,-0.68,-0.70679,-0.75009,-0.79699,-0.8346,-0.85,-0.83968,-0.81145,-0.76937,-0.71753,-0.66,-0.58223,-0.47111,-0.33287,-0.17376,-1.1102e-16}, {0,0.11,0.22,0.33,0.44,0.55,0.66,0.77,0.8504,0.9387,1.0141,1.0686,1.0971,1.0957,1.0656,1.0145,0.9567,0.8996,0.8566,0.8226,0.8089,0.8063,0.8189,0.8408,0.8668,0.9023,0.9406,0.9912,0.95589,0.89031,0.855,0.86556,0.89167,0.925,0.95722,0.98,0.99474,1.0081,1.0196,1.0287,1.035,1.0396,1.0437,1.047,1.0492,1.05,1.0482,1.0434,1.0365,1.0284,1.02,1.0103,0.99809,0.98427,0.96963,0.955,0.94052,0.92557,0.90986,0.8931,0.875,0.85487,0.83265,0.809,0.78456,0.76,0.73502,0.70915,0.68276,0.65625,0.63,0.60408,0.57824,0.55235,0.52631,0.5,0.47333,0.44637,0.41924,0.39208,0.365,0.33808,0.31124,0.28435,0.25731,0.23,0.20233,0.17437,0.14624,0.11808,0.09,0.061919,0.033756,0.0056335,-0.022326,-0.05,-0.077314,-0.10435,-0.13124,-0.15808,-0.185,-0.21208,-0.23924,-0.26637,-0.29333,-0.32,-0.3464,-0.3726,-0.3986,-0.4244,-0.45,-0.47586,-0.50199,-0.52769,-0.55226,-0.575,-0.59574,-0.61511,-0.63361,-0.65174,-0.67,-0.68831,-0.70635,-0.72423,-0.74208,-0.76,-0.77817,-0.79651,-0.81476,-0.83268,-0.85,-0.86719,-0.88445,-0.90111,-0.91651,-0.93,-0.94308,-0.95646,-0.96831,-0.97677,-0.98,-0.97485,-0.96135,-0.94242,-0.921,-0.9,-0.87742,-0.85083,-0.82253,-0.79482,-0.77,-0.74679,-0.7235,-0.70181,-0.68342,-0.67,-0.65972,-0.65021,-0.64234,-0.63698,-0.635,-0.6374,-0.64401,-0.65391,-0.66621,-0.68,-0.70679,-0.75009,-0.79699,-0.8346,-0.85,-0.83968,-0.81145,-0.76937,-0.71753,-0.66,-0.58223,-0.47111,-0.33287,-0.17376,-1.1102e-16}, {0,0.11,0.22,0.33,0.44,0.55,0.66,0.77,0.88,0.9574,1.0433,1.1138,1.1667,1.1948,1.1962,1.1744,1.1356,1.0921,1.051,1.0173,0.9954,0.9837,0.9827,0.991,1.0078,1.0317,1.0591,1.081,1.0224,0.91359,0.855,0.86556,0.89167,0.925,0.95722,0.98,0.99474,1.0081,1.0196,1.0287,1.035,1.0396,1.0437,1.047,1.0492,1.05,1.0482,1.0434,1.0365,1.0284,1.02,1.0103,0.99809,0.98427,0.96963,0.955,0.94052,0.92557,0.90986,0.8931,0.875,0.85487,0.83265,0.809,0.78456,0.76,0.73502,0.70915,0.68276,0.65625,0.63,0.60408,0.57824,0.55235,0.52631,0.5,0.47333,0.44637,0.41924,0.39208,0.365,0.33808,0.31124,0.28435,0.25731,0.23,0.20233,0.17437,0.14624,0.11808,0.09,0.061919,0.033756,0.0056335,-0.022326,-0.05,-0.077314,-0.10435,-0.13124,-0.15808,-0.185,-0.21208,-0.23924,-0.26637,-0.29333,-0.32,-0.3464,-0.3726,-0.3986,-0.4244,-0.45,-0.47586,-0.50199,-0.52769,-0.55226,-0.575,-0.59574,-0.61511,-0.63361,-0.65174,-0.67,-0.68831,-0.70635,-0.72423,-0.74208,-0.76,-0.77817,-0.79651,-0.81476,-0.83268,-0.85,-0.86719,-0.88445,-0.90111,-0.91651,-0.93,-0.94308,-0.95646,-0.96831,-0.97677,-0.98,-0.97485,-0.96135,-0.94242,-0.921,-0.9,-0.87742,-0.85083,-0.82253,-0.79482,-0.77,-0.74679,-0.7235,-0.70181,-0.68342,-0.67,-0.65972,-0.65021,-0.64234,-0.63698,-0.635,-0.6374,-0.64401,-0.65391,-0.66621,-0.68,-0.70679,-0.75009,-0.79699,-0.8346,-0.85,-0.83968,-0.81145,-0.76937,-0.71753,-0.66,-0.58223,-0.47111,-0.33287,-0.17376,-1.1102e-16}, {0,0.11,0.22,0.33,0.44,0.55,0.66,0.77,0.88,0.99,1.0685,1.1553,1.229,1.2847,1.3187,1.3298,1.3186,1.2917,1.2576,1.2242,1.1965,1.1771,1.1647,1.1611,1.1563,1.1322,1.1268,1.1397,1.0659,0.92881,0.855,0.86556,0.89167,0.925,0.95722,0.98,0.99474,1.0081,1.0196,1.0287,1.035,1.0396,1.0437,1.047,1.0492,1.05,1.0482,1.0434,1.0365,1.0284,1.02,1.0103,0.99809,0.98427,0.96963,0.955,0.94052,0.92557,0.90986,0.8931,0.875,0.85487,0.83265,0.809,0.78456,0.76,0.73502,0.70915,0.68276,0.65625,0.63,0.60408,0.57824,0.55235,0.52631,0.5,0.47333,0.44637,0.41924,0.39208,0.365,0.33808,0.31124,0.28435,0.25731,0.23,0.20233,0.17437,0.14624,0.11808,0.09,0.061919,0.033756,0.0056335,-0.022326,-0.05,-0.077314,-0.10435,-0.13124,-0.15808,-0.185,-0.21208,-0.23924,-0.26637,-0.29333,-0.32,-0.3464,-0.3726,-0.3986,-0.4244,-0.45,-0.47586,-0.50199,-0.52769,-0.55226,-0.575,-0.59574,-0.61511,-0.63361,-0.65174,-0.67,-0.68831,-0.70635,-0.72423,-0.74208,-0.76,-0.77817,-0.79651,-0.81476,-0.83268,-0.85,-0.86719,-0.88445,-0.90111,-0.91651,-0.93,-0.94308,-0.95646,-0.96831,-0.97677,-0.98,-0.97485,-0.96135,-0.94242,-0.921,-0.9,-0.87742,-0.85083,-0.82253,-0.79482,-0.77,-0.74679,-0.7235,-0.70181,-0.68342,-0.67,-0.65972,-0.65021,-0.64234,-0.63698,-0.635,-0.6374,-0.64401,-0.65391,-0.66621,-0.68,-0.70679,-0.75009,-0.79699,-0.8346,-0.85,-0.83968,-0.81145,-0.76937,-0.71753,-0.66,-0.58223,-0.47111,-0.33287,-0.17376,-1.1102e-16}, {0,0.11,0.22,0.33,0.44,0.55,0.66,0.77,0.88,0.99,1.1,1.1749,1.2591,1.33,1.3825,1.4136,1.4233,1.4136,1.3897,1.3608,1.3325,1.3077,1.2767,1.1981,1.1538,1.138,1.1374,1.1519,1.0749,0.93197,0.855,0.86556,0.89167,0.925,0.95722,0.98,0.99474,1.0081,1.0196,1.0287,1.035,1.0396,1.0437,1.047,1.0492,1.05,1.0482,1.0434,1.0365,1.0284,1.02,1.0103,0.99809,0.98427,0.96963,0.955,0.94052,0.92557,0.90986,0.8931,0.875,0.85487,0.83265,0.809,0.78456,0.76,0.73502,0.70915,0.68276,0.65625,0.63,0.60408,0.57824,0.55235,0.52631,0.5,0.47333,0.44637,0.41924,0.39208,0.365,0.33808,0.31124,0.28435,0.25731,0.23,0.20233,0.17437,0.14624,0.11808,0.09,0.061919,0.033756,0.0056335,-0.022326,-0.05,-0.077314,-0.10435,-0.13124,-0.15808,-0.185,-0.21208,-0.23924,-0.26637,-0.29333,-0.32,-0.3464,-0.3726,-0.3986,-0.4244,-0.45,-0.47586,-0.50199,-0.52769,-0.55226,-0.575,-0.59574,-0.61511,-0.63361,-0.65174,-0.67,-0.68831,-0.70635,-0.72423,-0.74208,-0.76,-0.77817,-0.79651,-0.81476,-0.83268,-0.85,-0.86719,-0.88445,-0.90111,-0.91651,-0.93,-0.94308,-0.95646,-0.96831,-0.97677,-0.98,-0.97485,-0.96135,-0.94242,-0.921,-0.9,-0.87742,-0.85083,-0.82253,-0.79482,-0.77,-0.74679,-0.7235,-0.70181,-0.68342,-0.67,-0.65972,-0.65021,-0.64234,-0.63698,-0.635,-0.6374,-0.64401,-0.65391,-0.66621,-0.68,-0.70679,-0.75009,-0.79699,-0.8346,-0.85,-0.83968,-0.81145,-0.76937,-0.71753,-0.66,-0.58223,-0.47111,-0.33287,-0.17376,-1.1102e-16}};
int flag_debug = 0;
int flag_shed = 0; /* 0 disattiva vortex shedding */

/* ============== FUNZIONI DI SUPPORTO ============== */
/* Funzione segno */
int signum(double x)
{
    if (x>0) return 1;
    if (x<0) return -1;
    return 0;
}

/* Funzione minimo fra due numeri */
double minimum(double x, double y)
{
    return (x<y) ? x : y;
}

/* Funzione massimo fra due numeri */
double maximum(double x, double y)
{
    return (x>y) ? x : y;
}

/* Funzione che calcola Cl statico */
double cl_static(double angolo, double Reynolds)
{
    double f_Re, f_attack, Cl, segno_angolo;
    int i, j;
    FILE *error;
    
    /* Controllo che il Reynolds sia compreso nel range. Se non lo e' viene impostato il valore consentito più vicino e viene prodotto un file di errore */
    if (Reynolds<Re_levels[0])
    {
        Reynolds=Re_levels[0];
    }
    if (Reynolds>Re_levels[4])
    {
        Reynolds=Re_levels[4];
    }
    
    /* Controllo che l'angolo sia compreso nel range. Se non lo e' viene limitato a +o-35 e viene prodotto un file di errore. Inoltre memorizzo il segno dell'angolo, di cui faccio il valore assoluto */
    segno_angolo=signum(angolo);
    angolo=fabs(angolo);
    
    if (angolo>180)
    {
        angolo=180;
    }
    
    /* Trovo f_Re, la posizione adimensionale rispetto all'intervallo di Reynolds in cui si trova il valore di input della funzione. Il primo intervallo e' 1 l'ultimo e' il 4*/
    for (j=1;j<5;j++)
    {
        if (Reynolds>=Re_levels[j-1] && Reynolds<=Re_levels[j])
        {
            f_Re=(Reynolds-Re_levels[j-1])/(Re_levels[j]-	Re_levels[j-1]);
            break;
        }
    }
    
    /* Trovo f_attack, la posizione adimensionale rispetto all'intervallo di Reynolds in cui si trova il valore di input della funzione */
    
    for (i=1;i<182;i++)
    {
        if (angolo>=attack[i-1] && angolo<=attack[i])
        {
            f_attack=(angolo-attack[i-1])/(attack[i]-attack[i-1]);
            break;
        }
    }
    
    /* Calcolo il Cl statico interpolando su angolo di attacco e Reynolds */
    
    Cl = cl_s[j-1][i-1]*(1-f_Re)*(1-f_attack) + cl_s[j][i-1]*f_Re*(1-f_attack) + cl_s[j-1][i]*(1-f_Re)*f_attack + cl_s[j][i]*f_Re*f_attack;
    
    /* Restituisco il valore del Cl statico con segno */
    return segno_angolo*Cl;
}
/* ============== END FUNZIONI DI SUPPORTO ============== */

/* Questa routine calcola l'intera curva del CL-alpha in condizioni dinamiche
 * ogni time-step. si entra con angolo massimo e minimo calcolati dal campo
 * di velocita', con omega, con il Reynolds e con il modulo della velocita'.
 * cosi' facendo siamo in grado di trovare la curva esatta che serve per il
 * caso in esame. In output abbiamo la curva intera, pertanto entrando con
 * alpha siamo in grado di calcolare il cl. */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Leggi da mxArray in input */
    double chord = mxGetScalar(prhs[0]);
    double omega = mxGetScalar(prhs[1]);
    double Uinf = mxGetScalar(prhs[2]);
    double modW = mxGetScalar(prhs[3]);
    double Re = mxGetScalar(prhs[4]);
    double alpha_in = mxGetScalar(prhs[5]);
    double alpha_in_0 = mxGetScalar(prhs[6]);
    double alpha_max = mxGetScalar(prhs[7]);
    double alpha_min = mxGetScalar(prhs[8]);
    double cosTh = mxGetScalar(prhs[9]);
    
    /* Dichiarazione delle variabili del flusso */
    double tau=chord/(2*modW);
    double k=(omega*chord)/(2*modW);
    double dt=(2*M_PI)/(1000*fabs(omega));
    double t_max=4*2*M_PI/fabs(omega);
    double U_lev=modW/(3*chord);
    
    /* Parametri del modello */
    double omega_3=0.08;
    double omega_5=0.12;
    double omega_4;
    double omega_6;
    double alpha_ds;
    double alpha_m;
    double cl_min;
    double puls=2*M_PI*(0.235*modW)/(chord);
    
    /* variabili locali */
    double cl_par0, cl_par, cl_stallo, aa0, aa, alpha_stallo;
    int Amp_init;
    int i, t, nt_in, count;
    int nt=t_max/dt+1;
    double clmax,clmin;
    double *al, *alpha_dot, *cls, *cl0s, *clns;
    double *x_revers, *x_lev;
    int *phase;
    double *alpha_d, *clsd, *cl0sd, *f, *fd, *cl0d;
    double *cld, *der_cl, *cl_vort, *cl_shed, *A, *cl_fin;
    double segno, cl_alpha, ka, kf;
    int imem, imem2;
    double foo, A_amp, A_med;
    double out;
    double rapp_incr,delta_alpha;
    double tau_omega_6;
    int i_vortex_start,vortex_start, reattached;
    char dyn[30];
    FILE *debug;
    char dyn1[30];
    FILE *debug1;
    char dyn2[30];
    FILE *debug2;
    
    /* Allocazioni ulteriori */
    int myid = 0;
    int N_TIME = 0;
    
    /* Alloco la memoria: */
    /* nt e' gia' dichiarato come intero, pertanto dovrebbe fare la divisione intera */
    al=(double * )malloc(nt*sizeof(double));
    alpha_dot=(double * )malloc(nt*sizeof(double));
    cls=(double * )malloc(nt*sizeof(double));
    cl0s=(double * )malloc(nt*sizeof(double));
    clns=(double * )malloc(nt*sizeof(double));
    phase=(int * )malloc(nt*sizeof(int));
    x_revers=(double * )malloc(nt*sizeof(double));
    x_lev=(double * )malloc(nt*sizeof(double));
    alpha_d=(double * )malloc(nt*sizeof(double));
    clsd=(double * )malloc(nt*sizeof(double));
    cl0sd=(double * )malloc(nt*sizeof(double));
    f=(double * )malloc(nt*sizeof(double));
    fd=(double * )malloc(nt*sizeof(double));
    cl0d=(double * )malloc(nt*sizeof(double));
    cld=(double * )malloc(nt*sizeof(double));
    der_cl=(double * )malloc(nt*sizeof(double));
    cl_vort=(double * )malloc(nt*sizeof(double));
    cl_shed=(double * )malloc(nt*sizeof(double));
    cl_fin=(double * )malloc(nt*sizeof(double));
    A=(double * )malloc(nt*sizeof(double));
    sprintf(dyn,"dynamic_%d.txt",myid);
    debug=fopen(dyn,"a");
    sprintf(dyn1,"stallo_%d.txt",myid);
    debug1=fopen(dyn1,"a");
    sprintf(dyn2,"interp_%d.txt",myid);
    debug2=fopen(dyn2,"a");
    
    /* ============== INPUT AND OUTPUT ARGS VALIDATION ============== */
    /* Check for proper number of input and output arguments */
    if (nrhs != 10)
    {
        mexErrMsgIdAndTxt( "cldynmex:invalidNumInputs",
                "9 input arguments are required.");
    }
    if (nlhs != 1)
    {
        mexErrMsgIdAndTxt( "cldynmex:invalidNumOutputs",
                "1 output arguments is required.");
    }
    
    /* Check for proper type of input arguments */
    for (i = 0; i < 9; i++)
    {
        if(!IS_double_SCALAR(prhs[i]))
        {
            mexErrMsgIdAndTxt( "cldynmex:invalidTypeInputs",
                    "One or more inputs arguments are not scalar.");
        }
    }
    
    /* ============== END INPUT AND OUTPUT ARGS VALIDATION ============== */
    if(flag_debug!=0)
        printf("chord=%f, omega=%f, Uinf=%f, modW=%f, Re=%f, alpha_in=%f, alpha_in_0=%f, alpha_max=%f, alpha_min=%f\n", chord, omega, Uinf, modW, Re, alpha_in, alpha_in_0, alpha_max, alpha_min);
    
    /* if(N_TIME==1 && flag_debug!=0)
	fprintf(debug," i  phase  al  alpha_d  cls  clns  cld  der_cl  A  cl_vort  cl_shed  cl_fin\n"); */
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Angolo di stallo statico, dalla curva cl-alpha in ingresso */
    cl_par0=0;
    aa=0;
    aa0=0;
    for(i=0;i<25;i++)
    {
        aa=aa+1;
        cl_par=cl_static(aa,Re);
        if(cl_par>cl_par0)
        {
            cl_par0=cl_par;
            aa0=aa;
        }
    }
    alpha_stallo=aa0;
    cl_stallo=cl_par0;
    
    /*  Ampiezza e valor medio del moto di pitching */
    A_amp=0.5*(alpha_max-alpha_min);
    A_med=alpha_max-A_amp;
    
    /*if(N_TIME<=96)
     * {
     *	A_amp=36.0;
     *	A_med=0.0;
     * }
     * else
     * {
     * 	A_amp=0.5*(alpha_max-alpha_min);
     * 	A_med=alpha_max-A_amp;
     * } */
    
    /* coefficiente di portanza lineare */
    cl_alpha=2*M_PI*M_PI/180;
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* calcolo alpha di stallo dinamico: */
    if(Re<490000)
        alpha_ds=320*(k*A_amp*M_PI/180)+12;
    else if (Re>=490000 && Re<980000)
        alpha_ds=320*(k*A_amp*M_PI/180)+14.5;
    else if (Re>=980000 && Re<2500000)
        alpha_ds=320*(k*A_amp*M_PI/180)+18;
    else
        alpha_ds=320*(k*A_amp*M_PI/180)+19;
    /* calcolo di alpha minimo: */
    if(Re<490000)
        alpha_m=-372*(k*A_amp*M_PI/180)+18;
    else if (Re>=490000 && Re<980000)
        alpha_m=-372*(k*A_amp*M_PI/180)+20;
    else if (Re>=980000 && Re<2500000)
        alpha_m=-372*(k*A_amp*M_PI/180)+22;
    else
        alpha_m=-372*(k*A_amp*M_PI/180)+25;
    /* calcolo del cl_minimo: */
    if(Re<=980000)
        cl_min=-15*(k*A_amp*M_PI/180)+0.75;
    else
        cl_min=-20*(k*A_amp*M_PI/180)+1.2;
    /* frequenza omega_6: */
    tau_omega_6=2*M_PI*(A_amp+A_med-alpha_m)/(4*A_amp*3*omega);
    /* frequenza omega_4: */
    omega_4=modW/3*tau;
	
    if(flag_debug!=0)
		fprintf(debug,"%d %f %f %f %f %f %f %f %f %f %f\n",N_TIME,alpha_ds,alpha_m,cl_min,omega_4,U_lev,tau_omega_6,(k*A_amp*M_PI/180),Re,alpha_max,alpha_min);
		
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Inizializzo  variabili allocate: */
    for(i=0;i<nt;i++)
    {
        al[i]=0;
        alpha_dot[i]=0;
        cls[i]=0;
        cl0s[i]=0;
        phase[i]=0;
        x_lev[i]=0;
        clsd[i]=0;
        cld[i]=0;
        der_cl[i]=0;
        cl_vort[i]=0;
        cl_shed[i]=0;
        cl_fin[i]=0;
        A[i]=0;
        i_vortex_start=0;
        vortex_start=0;
        reattached=0;
    }
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* inizio il calcolo per trovare la curva cl_dynamic-alpha */
    for(i=0;i<nt;i++)
    {
        /* moto sinusoidale */
        al[i]=A_med+A_amp*sin(omega*i*dt+0);
        /* derivata prima di alpha in rad/s */
        alpha_dot[i]=+A_amp*omega*cos(omega*i*dt+0)*M_PI/180;
        /* calcolo del coefficiente di portaza statico */
        segno=signum(al[i]);
        cls[i]=segno*cl_static(segno*al[i],Re);
        /* coefficiente di portanza lineare */
        cl0s[i]=cl_alpha*al[i];
        if(i==0)
        {
            cld[i]=cls[i];
            cl_vort[i]=0;
            cl_shed[i]=0;
        }
        else
        {
            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
            if(fabs(al[i])>fabs(al[i-1]))
            {
                reattached=0;
                if(vortex_start==0)
                {
                    x_lev[i]=0;
                    phase[i]=1;
                    cld[i]=cl0s[i];
                    der_cl[i]=cl0s[i]-cls[i];
                    cl_vort[i]=(cl_vort[i-1]+omega_3/tau*dt*der_cl[i])/(1+omega_3/tau*dt);
                    if(fabs(al[i])>alpha_ds)
                    {
                        vortex_start=1;
                        i_vortex_start=i;
                    }
                }
                else
                {
                    phase[i]=2;
                    cld[i]=(cld[i-1]+omega_4/tau*dt*cls[i])/(1+omega_4/tau*dt);
                    x_lev[i]=x_lev[i-1]+U_lev*dt;
                    if(x_lev[i-1]<1)
                    {
                        der_cl[i]=cld[i]-cls[i];
                        cl_vort[i]=(cl_vort[i-1]+omega_3/tau*dt*der_cl[i])/(1+omega_3/tau*dt);
                        cl_shed[i]=cl_vort[i]*sin(puls*(i-i_vortex_start)*dt);
                    }
                    else
                    {
                        der_cl[i]=0;
                        cl_vort[i]=(cl_vort[i-1]+omega_3/tau*dt*der_cl[i])/(1+omega_3/tau*dt);
                        cl_shed[i]=cl_vort[i]*sin(puls*(i-i_vortex_start)*dt);
                    }
                }/* su vortex start*/
            }
            else
            {
                if(fabs(al[i])>alpha_m || x_lev[i-1]<1 && phase[i-1]==2)
                {
                    phase[i]=3;
                    vortex_start=0;
                    reattached=0;
                    x_lev[i]=x_lev[i-1]+U_lev*dt;
                    tau_omega_6=2*M_PI*(A_amp+A_med-alpha_m)/(4*A_amp*3*omega);
                    omega_6=1/tau_omega_6;
                    segno=signum(al[i]);
                    cld[i]=(cld[i-1]+omega_6*dt*segno*cl_min)/(1+omega_6*dt);
                    der_cl[i]=0;
                    if((fabs(cl_shed[i-1])-fabs(cl_shed[i-2]))>0 && x_lev[i]>1)
                        cl_vort[i]=(cl_vort[i-1]+omega_3/tau*dt*der_cl[i])/(1+omega_3/tau*dt);
                    else
                        cl_vort[i]=(cl_vort[i-1]+8*omega_3/tau*dt*der_cl[i])/(1+8*omega_3/tau*dt);
					cl_shed[i]=cl_vort[i]*sin(puls*(i-i_vortex_start)*dt);
                }
                else
                {
                    if(reattached==0)
                        cld[i]=(cld[i-1]+omega_5/tau*dt*cls[i])/(1+omega_5/tau*dt);
                    else
                        cld[i]=cls[i];
                    x_lev[i]=x_lev[i-1]+U_lev*dt;
                    phase[i]=4;
                    der_cl[i]=0;
                    if(fabs(cld[i])>=fabs(cls[i]))
                    {
                        reattached=1;
                        cld[i]=cls[i];
                    }
                    cl_vort[i]=(cl_vort[i-1]+omega_3/tau*dt*der_cl[i])/(1+omega_3/tau*dt);
					cl_shed[i]=cl_vort[i]*sin(puls*(i-i_vortex_start)*dt);
                } /* su alpha_min */
            } /* su check up/down */
			
			if(flag_shed == 0)
				cl_shed[i] = 0;				
			
            cl_fin[i]=cld[i]+cl_shed[i];
            /*if(i>=(3/4*nt))
            *{
            *    fprintf(debug1,"%d %d %f %f %f %f %f %f\n",N_TIME,phase[i],modW,al[i],cl_fin[i],cld[i],cl_shed[i],cl_vort[i]);
            *    
            * } */
        }/* su esclusione primo giro */
    } /* su for i */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Smoothinge della curva cl_fin-alpha. faccio quest'operazione perche'
     * per ogni aggiornamento delle condizioni in ingresso della routine
     * dovrei ottimizzare la curva in uscita con un processo iterativo.
     * questo metodo non e' possibile durante il calcolo della soluzione
     * delle equazioni di NS. Per eliminare le discontinuita', dovute alla
     * non ottimizzazione dei parametri liberi del modello, applico uno
     * smooth della curva: */
    /* Smothing per 3 punti, lasciando gli estremi fuori: */
    /* lo faccio solo per l'ultimo giro perche' e' la curva che effetti-
     * vamente usero' per il calcolo del cl */
    nt_in=3*nt/4;
    for(i=(nt_in+1);i<(nt-1);i++)
        /*for(i=nt_in;i<nt;i++)*/
    {
        /* if(flag_debug!=0) */
        /* fprintf(debug2,"%d %d %f %f\n", N_TIME, i, al[i], cl_fin[i]); */
        
        delta_alpha=al[i+1]-al[i-1];
        if(fabs(delta_alpha)<0.001)
        {
            cl_fin[i]=0.5*(cl_fin[i+1]+cl_fin[i-1]);
        }
        else
        {
            rapp_incr=(cl_fin[i+1]-cl_fin[i-1])/(delta_alpha);
            cl_fin[i]=cl_fin[i-1]+rapp_incr*(al[i]-al[i-1]);
        }
        if(al[i]==(A_amp+A_med))
            clmax=cl_fin[i];
        else if(al[i]==(A_med-A_amp))
            clmin=cl_fin[i];
    }
    /* scrivo gli output: */
    segno=signum(alpha_in);
    foo=0.;
    if(alpha_in>=(A_med+A_amp))
        foo=clmax;
    else if (alpha_in<=(A_med-A_amp))
        foo=clmin;
    else if(alpha_in>=alpha_in_0) /* UPSTROKE: */
    {
        for(i=(nt_in);i<(nt+1);i++)
        {
            if(al[i-1]<=alpha_in && al[i]>alpha_in && al[i]>=al[i-1])
            {
                foo=((cl_fin[i]-cl_fin[i-1])/(al[i]-al[i-1])*(alpha_in-al[i-1])+cl_fin[i-1]);
                break;
            }
        }
    }
    else if(alpha_in<alpha_in_0) /* DOWNSTROKE: */
    {
        for(i=(nt_in);i<(nt+1);i++)
        {
            if(al[i]<=alpha_in && al[i-1]>alpha_in && al[i]<al[i-1])
            {
                foo=((cl_fin[i]-cl_fin[i-1])/(al[i]-al[i-1])*(alpha_in-al[i-1])+cl_fin[i-1]);
                break;
            }
        }
    }
    else
    {
        fprintf(debug2,"ERROR definition sign(a): N_TIME=%d alpha=%f alpha_in_0=%f out=%f cosTh=%f prova=%f\n",N_TIME, alpha_in, alpha_in_0, out, cosTh, alpha_max);
    }
    
    out=foo;
    
    if(out==0)
    {
        fprintf(debug2,"Error in the output evaluation, no interpolation for:\n");
        fprintf(debug2,"N_TIME=%d alpha=%f alpha_in_0=%f out=%f cosTh=%f prova=%f\n",N_TIME, alpha_in, alpha_in_0, out, cosTh, alpha_max);
        for(i=nt_in;i<nt;i++)
        {
            fprintf(debug2,"%f %f\n", al[i],cl_fin[i]);
        }
    }
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Libero le variabili: */
    free(al);
    free(alpha_dot);
    free(cls);
    free(cl0s);
    free(clns);
    free(phase);
    free(x_revers);
    free(x_lev);
    free(alpha_d);
    free(clsd);
    free(cl0sd);
    free(f);
    free(fd);
    free(cl0d);
    free(cld);
    free(der_cl);
    free(cl_vort);
    free(cl_shed);
    free(cl_fin);
    free(A);
    fclose(debug);
    fclose(debug1);
    fclose(debug2);
    
    /* Pointer da restituire per output */
    plhs[0] = mxCreateDoubleScalar(out);
    
    return;
}