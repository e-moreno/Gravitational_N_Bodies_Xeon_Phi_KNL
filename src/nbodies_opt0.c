#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <sys/time.h>

#ifndef DATATYPE
  #define DATATYPE float
#endif
#ifndef POW
    #define POW pow
#endif
#ifndef SQRT
    #define SQRT sqrt
#endif

#define VERSION "SEC-V1EM"
#ifdef PRINTBODIES
  #define PRINT1 printf("%ld, %ld\n",N,D);
#else
  #define PRINT1 ;
#endif
#ifdef PRINTBODIES
  #define PRINT2 printf("%.2lf, %.2lf, %.2lf \n", xpos[i],ypos[i],zpos[i]);
#else
  #define PRINT2 ;
#endif

#define EPS 20E3
//Constante universal de gravitación
const DATATYPE G = 6.674e-11;
//Diferencial de tiempo 1 segundo
const DATATYPE DT = 1;
const DATATYPE SOFT = 1e-20;

//Para calcular tiempo
double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

DATATYPE *masses;
DATATYPE *xpos;
DATATYPE *ypos;
DATATYPE *zpos;
DATATYPE *xvi;
DATATYPE *yvi;
DATATYPE *zvi;
DATATYPE *forcesx;
DATATYPE *forcesy;
DATATYPE *forcesz;

int main(int argc, char *argv[]){
    long N;
    long D;
    int num_threads=1;
    long t, i, j;
    if (argc < 3){
        printf("invalid arguments, use %s N Duration\n", argv[0]);
        return 0;
    }else{
        N = atoi(argv[1]);
        D = atoi(argv[2]);
    }
    masses = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    xpos = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    ypos = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    zpos = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    xvi = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    yvi = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    zvi = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    forcesx = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    forcesy = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    forcesz = (DATATYPE *)malloc(N * sizeof(DATATYPE));

    int sqrt_n = sqrt(N);
    if(!((sqrt_n * sqrt_n) == N)){
      sqrt_n++;
    }
    int dist = 100000;

    for (i = 0; i < N; i++){
        xvi[i] = 0;
        yvi[i] = 0;
        zvi[i] = 0;
        forcesx[i] = 0;
        forcesy[i] = 0;
        forcesz[i] = 0;
        xpos[i] = (i%sqrt_n)*dist;
        ypos[i] = dist * i;
        zpos[i] = 5000;
        if(i == 6){
          masses[i] = 5.97e20;
        }else{
          masses[i] = 5.97e20;
        }
    }
    PRINT1
    double timetick = dwalltime();
    for (t = 1; t <= D; t += DT){
        for (i = 0; i < N; i++){
            PRINT2
            forcesx[i] = 0.0;
            forcesy[i] = 0.0;
            forcesz[i] = 0.0;
            for (j = 0; j < N; j++){
                const DATATYPE dx = xpos[j] - xpos[i];
                const DATATYPE dy = ypos[j] - ypos[i];
                const DATATYPE dz = zpos[j] - zpos[i];
                // d²
                const DATATYPE dsquared = (dx*dx) + (dy*dy) + (dz*dz) + SOFT;
                // Newton 's gravitation.
                // F = (G*m1*m2) * 1/d²
                // Calculating (G*m1*m2) separatedly
                const DATATYPE F = G * masses[i] * masses[j];
                // Calculating 1/d³ separatedly
                const DATATYPE d32 = 1/POW(dsquared,3.0/2.0);

                // Unit vector x,y,z coordinates
                // dx/d, dy/d, dz/d
                // Force vector x,y,z coordinates
                // F * 1/d² * [dx | dy | dz ] * 1/d
                // F * dx * 1/d² * 1/d ==> F * dx * 1/d³
                forcesx[i] += F*dx*d32;
                // F * dy * 1/d² * 1/d ==> F * dy * 1/d³
                forcesy[i] += F*dy*d32;
                // F * dz * 1/d² * 1/d ==> F * dz * 1/d³
                forcesz[i] += F*dz*d32;
            }
            // Vector aceleración
            // |M = F * A| ---> |A = F / M|
            const DATATYPE ax = forcesx[i] / masses[i];
            const DATATYPE ay = forcesy[i] / masses[i];
            const DATATYPE az = forcesz[i] / masses[i];
            // Velocity verlet integrator
            const DATATYPE dvx = (xvi[i] + (ax*DT/2));
            const DATATYPE dvy = (yvi[i] + (ay*DT/2));
            const DATATYPE dvz = (zvi[i] + (az*DT/2));
            // Position difference delta p
            const DATATYPE dpx =  dvx * DT;
            const DATATYPE dpy =  dvy * DT;
            const DATATYPE dpz =  dvz * DT;
            // Update velocity
            xvi[i] = dvx;
            yvi[i] = dvy;
            zvi[i] = dvz;
            // Update position
            xpos[i] += dpx;
            ypos[i] += dpy;
            zpos[i] += dpz;
        }
    }
	const double tw = dwalltime() - timetick;
    printf("%s;%d;%d;%lf;%lf;%lf\n",argv[0],N,D, tw, ((double)D*(16*N*N+18*N))/(tw*1000000000), ((double)D*(20*N*N))/(tw*1000000000));
    free(masses);
    free(xpos);
    free(ypos);
    free(zpos);
    free(xvi);
    free(yvi);
    free(zvi);
    free(forcesx);
    free(forcesy);
    free(forcesz);
}
