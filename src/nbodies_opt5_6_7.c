#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>

#ifndef DATATYPE
  #define DATATYPE float
#endif
#ifndef POW
    #define POW powf
#endif
#ifndef SQRT
    #define SQRT sqrtf
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
DATATYPE *dpx;
DATATYPE *dpy;
DATATYPE *dpz;

int main(int argc, char *argv[]){
    long N;
    long D;
    int num_threads;
    long t, i, j;
    if (argc < 4){
        printf("invalid arguments, use %s N Iterations Threads\n", argv[0]);
        return 0;
    }else{
        N = atoi(argv[1]);
        D = atoi(argv[2]);
        num_threads = atoi(argv[3]);
    }
    omp_set_num_threads(num_threads);
  	omp_set_nested(1);

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
    dpx = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    dpy = (DATATYPE *)malloc(N * sizeof(DATATYPE));
    dpz = (DATATYPE *)malloc(N * sizeof(DATATYPE));

    int sqrt_n = sqrt(N);
    if(!((sqrt_n * sqrt_n) == N)){
      sqrt_n++;
    }
    int dist = 100000;

    for (i = 0; i < N; i++){
        xvi[i] = 0.0f;
        yvi[i] = 0.0f;
        zvi[i] = 0.0f;
        forcesx[i] = 0.0f;
        forcesy[i] = 0.0f;
        forcesz[i] = 0.0f;
        xpos[i] = (i%sqrt_n)*dist;
        ypos[i] = dist * i;
        zpos[i] = 5000.0f;
        if(i == 6){
          masses[i] = 5.97e20f;
        }else{
          masses[i] = 5.97e20f;
        }
    }
    PRINT1
    double timetick = dwalltime();
    for (t = 1; t <= D; t += DT){
        #pragma omp parallel for schedule(static) private(j)
        for(i = 0; i < N; i++){
            PRINT2
            forcesx[i] = 0.0f;
            forcesy[i] = 0.0f;
            forcesz[i] = 0.0f;
            #pragma omp simd
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

                const DATATYPE d12 = 1.0f/SQRT(dsquared);
                // Calculating 1/d³ separatedly
                const DATATYPE d32 = d12 * d12 * d12;

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
            const DATATYPE dvx = (xvi[i] + (ax*DT*0.5f));
            const DATATYPE dvy = (yvi[i] + (ay*DT*0.5f));
            const DATATYPE dvz = (zvi[i] + (az*DT*0.5f));
            // Position difference delta p
            dpx[i] =  dvx * DT;
            dpy[i] =  dvy * DT;
            dpz[i] =  dvz * DT;
            // Update velocity
            xvi[i] = dvx;
            yvi[i] = dvy;
            zvi[i] = dvz;
        }
        #pragma omp parallel for simd schedule(static)
        for(i = 0; i < N; i++){
          // Update position
          xpos[i] += dpx[i];
          ypos[i] += dpy[i];
          zpos[i] += dpz[i];
        }
    }
	const double tw = dwalltime() - timetick;
    printf("%s;%d;%d;%d;%lf;%lf;%lf\n",argv[0],N,D,num_threads, tw, ((double)D*(16*N*N+18*N))/(tw*1000000000), ((double)D*(20*N*N))/(tw*1000000000));

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
    free(dpx);
    free(dpy);
    free(dpz);
}
