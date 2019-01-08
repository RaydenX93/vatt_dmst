#include <math.h>
#include "mex.h"

/* Definizioni per controllo input */
#define IS_double_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_double_SCALAR(P) (IS_double_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)

/* ============== DICHIARAZIONE VARIABILI GLOBALI ============== */
double attack[59] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180};
double Re_levels[9] = {40000,80000,160000,360000,700000,1000000,2000000,5000000,10000000};
double cl_s[9][59]={{0,0.1054,0.2099,0.3078,0.4017,0.4871,0.5551,0.5730,0.4663,0.04330,-0.04130,-0.01440,0.02610,0.07410,0.1244,0.1756,0.2280,0.2815,0.3351,0.3889,0.4427,0.4966,0.5506,0.6045,0.6585,0.7125,0.7666,0.8222,0.8550,0.9800,1.035,1.050,1.020,0.9550,0.8750,0.7600,0.6300,0.5000,0.3650,0.2300,0.09000,-0.05000,-0.1850,-0.3200,-0.4500,-0.5750,-0.6700,-0.7600,-0.8500,-0.9300,-0.9800,-0.9000,-0.7700,-0.6700,-0.6350,-0.6800,-0.8500,-0.6600,0},{0,0.1100,0.2200,0.3300,0.4186,0.5180,0.6048,0.6760,0.7189,0.6969,0.5122,0.1642,0.07490,0.09670,0.1382,0.1861,0.2364,0.2873,0.3393,0.3927,0.4463,0.5001,0.5539,0.6078,0.6617,0.7156,0.7700,0.8277,0.8550,0.9800,1.035,1.050,1.020,0.9550,0.8750,0.7600,0.6300,0.5000,0.3650,0.2300,0.09000,-0.05000,-0.1850,-0.3200,-0.4500,-0.5750,-0.6700,-0.7600,-0.8500,-0.9300,-0.9800,-0.9000,-0.7700,-0.6700,-0.6350,-0.6800,-0.8500,-0.6600,0},{0,0.1100,0.2200,0.3300,0.4400,0.5500,0.6299,0.7150,0.7851,0.8311,0.8322,0.7623,0.5936,0.3549,0.2371,0.2376,0.2665,0.3098,0.3567,0.4066,0.4575,0.5087,0.5611,0.6148,0.6685,0.7224,0.7771,0.8382,0.8550,0.9800,1.035,1.050,1.020,0.9550,0.8750,0.7600,0.6300,0.5000,0.3650,0.2300,0.09000,-0.05000,-0.1850,-0.3200,-0.4500,-0.5750,-0.6700,-0.7600,-0.8500,-0.9300,-0.9800,-0.9000,-0.7700,-0.6700,-0.6350,-0.6800,-0.8500,-0.6600,0},{0,0.1100,0.2200,0.3300,0.4400,0.5500,0.6600,0.7390,0.8240,0.8946,0.9440,0.9572,0.9285,0.8562,0.7483,0.6350,0.5384,0.4851,0.4782,0.4908,0.5247,0.5616,0.6045,0.6528,0.7015,0.7511,0.8055,0.8788,0.8550,0.9800,1.035,1.050,1.020,0.9550,0.8750,0.7600,0.6300,0.5000,0.3650,0.2300,0.09000,-0.05000,-0.1850,-0.3200,-0.4500,-0.5750,-0.6799,-0.7600,-0.8500,-0.9300,-0.9800,-0.9000,-0.7700,-0.6700,-0.6350,-0.6800,-0.8500,-0.6600,0},{0,0.1100,0.2200,0.3300,0.4400,0.5500,0.6600,0.7483,0.8442,0.9260,0.9937,1.03630000000000,1.05080000000000,1.03020000000000,0.9801,0.9119,0.8401,0.7799,0.7305,0.7041,0.6990,0.7097,0.7298,0.7593,0.7961,0.8353,0.8838,0.9473,0.8550,0.9800,1.035,1.050,1.020,0.9550,0.8750,0.7600,0.6300,0.5000,0.3650,0.2300,0.09000,-0.05000,-0.1850,-0.3200,-0.4500,-0.5750,-0.6700,-0.7600,-0.8500,-0.9300,-0.9800,-0.9000,-0.7700,-0.6700,-0.6350,-0.6800,-0.8500,-0.6600,0},{0,0.1100,0.2200,0.3300,0.4400,0.5500,0.6600,0.7700,0.8504,0.9387,1.01410000000000,1.06860000000000,1.09710000000000,1.09570000000000,1.06560000000000,1.01450000000000,0.9567,0.8996,0.8566,0.8226,0.8089,0.8063,0.8189,0.8408,0.8668,0.9023,0.9406,0.9912,0.8550,0.9800,1.035,1.050,1.020,0.9550,0.8750,0.7600,0.6300,0.5000,0.3650,0.2300,0.09000,-0.05000,-0.1850,-0.3200,-0.4500,-0.5750,-0.6700,-0.7600,-0.8500,-0.9300,-0.9800,-0.9000,-0.7700,-0.6700,-0.6350,-0.6800,-0.8500,-0.6600,0},{0,0.11000,0.22000,0.33000,0.44000,0.55000,0.66000,0.77000,0.88000,0.95740,1.0433,1.1138,1.1667,1.1948,1.1962,1.1744,1.1356,1.0921,1.0510,1.0173,0.99540,0.98370,0.98270,0.99100,1.0078,1.0317,1.0591,1.0810,0.85500,0.98000,1.0350,1.0500,1.0200,0.95500,0.87500,0.76000,0.63000,0.50000,0.36500,0.23000,0.090000,-0.050000,-0.18500,-0.32000,-0.45000,-0.57500,-0.67000,-0.76000,-0.85000,-0.93000,-0.98000,-0.90000,-0.77000,-0.67000,-0.63500,-0.68000,-0.85000,-0.66000,0},{0,0.11000,0.22000,0.33000,0.44000,0.55000,0.66000,0.77000,0.88000,0.99000,1.0685,1.1553,1.2290,1.2847,1.3187,1.3298,1.3186,1.2917,1.2576,1.2242,1.1965,1.1771,1.1647,1.1611,1.1563,1.1322,1.1268,1.1397,0.85500,0.98000,1.0350,1.0500,1.0200,0.95500,0.87500,0.76000,0.63000,0.50000,0.36500,0.23000,0.090000,-0.050000,-0.18500,-0.32000,-0.45000,-0.57500,-0.67000,-0.76000,-0.85000,-0.93000,-0.98000,-0.90000,-0.77000,-0.67000,-0.63500,-0.68000,-0.85000,-0.66000,0},{0,0.11000,0.22000,0.33000,0.44000,0.55000,0.66000,0.77000,0.88000,0.99000,1.1000,1.1749,1.2591,1.3300,1.3825,1.4136,1.4233,1.4136,1.3897,1.3608,1.3325,1.3077,1.2767,1.1981,1.1538,1.1380,1.1374,1.1519,0.85500,0.98000,1.0350,1.0500,1.0200,0.95500,0.87500,0.76000,0.63000,0.50000,0.36500,0.23000,0.090000,-0.050000,-0.18500,-0.32000,-0.45000,-0.57500,-0.67000,-0.76000,-0.85000,-0.93000,-0.98000,-0.90000,-0.77000,-0.67000,-0.63500,-0.68000,-0.85000,-0.66000,0}};

int flag_debug = 0;
int myid = 0;
int CURRENT_TIME = 0;

double omega_1=0.1, omega_2=0.1, omega_3=0.1, omega_4=0.1, omega_f4=0.1, omega_f3=0.1, alpha_star=14, alpha_lev=18;

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
    if (Reynolds>Re_levels[8])
    {
        Reynolds=Re_levels[8];
    }
    
    /* Controllo che l'angolo sia compreso nel range. Se non lo e' viene limitato a +o-35 e viene prodotto un file di errore. Inoltre memorizzo il segno dell'angolo, di cui faccio il valore assoluto */
    segno_angolo=signum(angolo);
    angolo=fabs(angolo);
    
    if (angolo>180)
    {
        angolo=180;
    }
    
    /* Trovo f_Re, la posizione adimensionale rispetto all'intervallo di Reynolds in cui si trova il valore di input della funzione. Il primo intervallo e' 1 l'ultimo e' il 4*/
    for (j=1;j<9;j++)
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Leggi da mxArray in input */
    double modW = mxGetScalar(prhs[0]);
    
    double alpha_d_0 = mxGetScalar(prhs[1]);
    double f_d_0 = mxGetScalar(prhs[2]);
    double cl_forc_0 = mxGetScalar(prhs[3]);
    double cl_vortex_0 = mxGetScalar(prhs[4]);
    int flag_lev_0 = mxGetScalar(prhs[5]);
    double x_lev_0 = mxGetScalar(prhs[6]);
    double cl_circ_0 = mxGetScalar(prhs[7]);
    double alpha_0 = mxGetScalar(prhs[8]);
    
    double alpha = mxGetScalar(prhs[9]);
    double Re = mxGetScalar(prhs[10]);
    double n_ring = mxGetScalar(prhs[11]);
    double omega = mxGetScalar(prhs[12]);
    double chord = mxGetScalar(prhs[13]);
    
    /* Dichiarazione variabili */
    double tau = chord/(2*modW);
    /* double Re = rho*modW*chord/mu; */
    double dt=(2*M_PI)/(n_ring*fabs(omega));
    double cl_circ; /*cl ritardo circuitazione*/
    double cl_0d, cl_fs, cl_cd, f, f_d, k_f, k_a, cl_s1, cl_alpha;
    double alpha_d, cl_sd, cl_0sd;
    double x_lev;
    double cl_d, cl_forc, der_cl, cl_vortex, cl_fin;
    int phase, flag_lev;
    int i,j;
    int stronzo2;
    double f_attack, f_Re, segno_angolo, angolo;
    
    /* Calcoli preliminari */
    /* Per cl_statico */
    /* double f_attack, f_Re, segno_angolo, angolo; */
    /* Derivata pitch motion */
    double alpha_dot = (alpha-alpha_0)/dt;
    
    /* Dichiarazione file */
    char nname[]="debug_dyn_out.txt";
    FILE *deb=fopen(nname,"a");
    
    char filename[]="debug_dyn_warning.txt";
    FILE *debug=fopen(filename,"a");
    
    char filename1[]="debug_dyn_in.txt";
    FILE *debug1=fopen(filename1,"a");
    
    /* ============== INPUT AND OUTPUT ARGS VALIDATION ============== */
    /* Check for proper number of input and output arguments */
    if (nrhs != 14)
    {
        mexErrMsgIdAndTxt( "cldynmex:invalidNumInputs",
                "14 input arguments are required.");
    }
    if (nlhs != 9)
    {
        mexErrMsgIdAndTxt( "cldynmex:invalidNumOutputs",
                "9 output arguments are required.");
    }
    
    /* Check for proper type of input arguments */
    for (i = 0; i < 14; i++)
    {
        if(!IS_double_SCALAR(prhs[i]))
        {
            mexErrMsgIdAndTxt( "cldynmex:invalidTypeInputs",
                    "One or more inputs arguments are not scalar.");
        }
    }
    
    /* ============== END INPUT AND OUTPUT ARGS VALIDATION ============== */
	
    if (flag_debug!=0)
	printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", modW, alpha_d_0, f_d_0, cl_forc_0, cl_vortex_0, flag_lev_0, x_lev_0, cl_circ_0, alpha_0, alpha, Re, n_ring, omega, chord);
	
    /* if (flag_debug!=0)
     * {
     *    fprintf(debug1, "clfin0=%f, ad0=%f, fd0=%f, clforc0=%f, clvort0=%f, flaglev0=%d, xlev0=%f, clcirc0=%f, a0=%f\n", input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7], input[8]);
     * }
     * fclose(debug1); */
    
    /* Check valori in ingresso (probabilmente inutile) */
    stronzo2 = 0;
    if (alpha_d_0 != alpha_d_0 || alpha_d_0 > 10000 || alpha_d_0 < -10000)
    {
        alpha_d_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (f_d_0!= f_d_0 || f_d_0 > 10000 || f_d_0 < -10000)
    {
        f_d_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (cl_forc_0 != cl_forc_0 || cl_forc_0 > 10000 || cl_forc_0 < -10000)
    {
        cl_forc_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (cl_vortex_0 != cl_vortex_0 || cl_vortex_0 > 10000 || cl_vortex_0 < -10000)
    {
        cl_vortex_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (flag_lev_0 != flag_lev_0 || flag_lev_0 > 10000 || flag_lev_0 < -10000)
    {
        flag_lev_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (cl_vortex_0 != cl_vortex_0 || cl_vortex_0 > 10000 || cl_vortex_0 < -10000)
    {
        cl_vortex_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (x_lev_0 != x_lev_0 || x_lev_0 > 10000 || x_lev_0 < -10000)
    {
        x_lev_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (cl_circ_0 != cl_circ_0 || cl_circ_0 > 10000 || cl_circ_0 < -10000)
    {
        cl_circ_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (alpha_0 != alpha_0 || alpha_0 > 10000 || alpha_0 < -10000)
    {
        alpha_0 = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (alpha_dot < -10000 || alpha_dot > 10000)
    {
        alpha_dot = 0;
        stronzo2 = stronzo2 + 1;
    }
    
    if (stronzo2 != 0)
    {
        fprintf(debug,"time=%f, stronzo2=%i\n",CURRENT_TIME,stronzo2);
    }
    fclose(debug);
    
    /* Calcolo il Cl statico e il Cl_alpha */
	segno_angolo=signum(alpha);
    cl_s1 = segno_angolo*cl_static(alpha,Re);
    cl_alpha = 2*M_PI*M_PI/180;
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * CIRCULATION DELAY
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Ritardo nella formazione della circuitazione, equazione differenziale ricavata dall'integrale di convoluzione, proposto da Larsen */
    cl_circ = (cl_circ_0+(cl_alpha*alpha_dot*dt))/(1+(omega_1*dt/tau));
    cl_0d = cl_alpha*alpha-0.5*cl_circ;
    /* contributo al cl dello strato limite completamente separato */
    cl_fs = (M_PI*alpha_dot*chord)/(2*modW);
    /* CL totale per il ritardo nella circuitazione */
    cl_cd=cl_0d;
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * FASI DEL LEADING EDGE VORTEX:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    if (alpha_dot*signum(alpha)>=0 && alpha*signum(alpha)<alpha_lev)
    {
        phase=1;
        x_lev=0;
        flag_lev=0;
    }
    else if (alpha*signum(alpha)>alpha_lev && flag_lev_0==0)
    {
        phase=2;
        x_lev=x_lev_0+(modW*dt)/(3*chord);
        x_lev=minimum(1,x_lev);
        if (f_d_0<0.05 || x_lev==1) flag_lev=1;
        else flag_lev=0;
    }
    else if (flag_lev_0==1)
    {
        if (alpha*signum(alpha)>alpha_star)
        {
            phase=3;
            x_lev=1;
            flag_lev=1;
        }
        else if (alpha*signum(alpha)<alpha_star && alpha_dot*signum(alpha_dot)<0)
        {
            phase=4;
            x_lev=1;
            flag_lev=1;
        }
    }
    else
    {
        phase=4;
        x_lev=0;
        flag_lev=0;
    }
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * GRADO DI ATTACCO DELLO SL:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    if (phase==1 || phase==4)
    {
        k_a=omega_2;
        /* calcolo del ritardo su alpha */
        alpha_d=(alpha_d_0+(k_a*alpha*dt/tau))/(1+(k_a*dt/tau));
        /* cl statico calcolato con alpha ritardato */
		segno_angolo=signum(alpha_d);
        cl_sd = segno_angolo*cl_static(alpha_d,Re);
        cl_0sd=cl_alpha*alpha_d;
        /* grado di attacco dello strato limite statico, ma calcolato con alpha ritardato */
        f=(2*sqrt(fabs(cl_sd/cl_0sd))-1)*(2*sqrt(fabs(cl_sd/cl_0sd))-1);
        /* controllo dei limiti di f */
        f=minimum(f,1);
        f=maximum(f,0);
        /* calcolo del grado di attacco dinamico */
        if (phase==1) f_d=f;
        else if (phase==4)
        {
            k_f=omega_f4;
            f_d=(f_d_0+(k_f*dt/tau*f))/(1+k_f*dt/tau);
        }
        /* controllo dei limiti di f dinamico */
        f_d=minimum(f_d,1);
        f_d=maximum(f_d,0);
    }
    else if (phase==2)
        /* la frequenza in questa fase e' dettata dalla velocita' di convezione del LEV */
    {
        k_a=modW/(3*chord);
        alpha_d=(alpha_d_0+(k_a*alpha*dt))/(1+(k_a*dt));
        /* cl statico calcolato con alpha ritardato */
		segno_angolo=signum(alpha_d);
        cl_sd = segno_angolo*cl_static(alpha_d,Re);
        cl_0sd=cl_alpha*alpha_d;
        f=(2*sqrt(fabs(cl_sd/cl_0sd))-1)*(2*sqrt(fabs(cl_sd/cl_0sd))-1);
        f=minimum(f,1);
        f=maximum(f,0);
        /* HP) la separazione dinamica dello strato limite e' imposta dal moto del LEV. Questo ha velocita' costante, pertanto viene assunto lo strato limite si separi linearmente */
        f_d=f_d_0-k_a*dt;
        f_d=minimum(f_d,1);
        f_d=maximum(f_d,0);
    }
    else if (phase==3)
    {
        k_a=omega_3;
        alpha_d=(alpha_d_0+(k_a*alpha*dt/tau))/(1+(k_a*dt/tau));
        /* cl statico calcolato con alpha ritardato */
		segno_angolo=signum(alpha_d);
        cl_sd = segno_angolo*cl_static(alpha_d,Re);
        cl_0sd=cl_alpha*alpha_d;
        f=(2*sqrt(fabs(cl_sd/cl_0sd))-1)*(2*sqrt(fabs(cl_sd/cl_0sd))-1);
        f=minimum(f,1);
        f=maximum(f,0);
        k_f=omega_f3;
        f_d=(f_d_0+(k_f*dt/tau*f))/(1+k_f*dt/tau);
        f_d=minimum(f_d,1);
        f_d=maximum(f_d,0);
    }
    /* Calcolo del CL dinamico senza il contributo del LEV */
    cl_d=cl_cd*((1+sqrt(f_d))/2)*((1+sqrt(f_d))/2);
    
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * LEADING EDGE VORTEX
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Forzante che energizza il LEV */
    cl_forc=cl_0d-cl_d;
    
    if (phase==1 || phase==2) der_cl=cl_forc-cl_forc_0;
    else der_cl=0;
    
    /* calcolo della portanza data dal LEV */
    cl_vortex=(cl_vortex_0+der_cl)/(1+omega_4*dt/tau);
    /* Calcolo del Cl finale */
    cl_fin=cl_d+cl_vortex;
	
	if (flag_debug!=0)
	printf("%f %f %f %f %f %d %f %f %f\n", cl_fin, alpha_d, f_d, cl_forc, cl_vortex, flag_lev, x_lev, cl_circ, alpha);
    
    /* Scrittura output */
    plhs[0]=mxCreateDoubleScalar(cl_fin);
    plhs[1]=mxCreateDoubleScalar(alpha_d);
    plhs[2]=mxCreateDoubleScalar(f_d);
    plhs[3]=mxCreateDoubleScalar(cl_forc);
    plhs[4]=mxCreateDoubleScalar(cl_vortex);
    plhs[5]=mxCreateDoubleScalar(flag_lev);
    plhs[6]=mxCreateDoubleScalar(x_lev);
    plhs[7]=mxCreateDoubleScalar(cl_circ);
    plhs[8]=mxCreateDoubleScalar(alpha);
    
    /* if (flag_debug!=0)
     * {
     *     fprintf(deb, "clf=%f, ad=%f, fd=%f, clfor=%f, clvortx=%f, fl=%f, xl=%f, clcirc=%f, a=%f\n", foo[0], foo[1], foo[2], foo[3], foo[4], foo[5], foo[6], foo[7], foo[8]);
     * } */
    
    fclose(deb);
	fclose(debug);
	fclose(debug1);
    
    return;
}