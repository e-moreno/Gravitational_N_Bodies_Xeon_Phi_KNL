#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <malloc.h>

#ifndef ALIGNMENT
  #define ALIGNMENT 64
#endif

#ifndef BLOCKSIZE
  #define BLOCKSIZE 16
#endif

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
    long t, i, j, ii;
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

    masses = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    xpos = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    ypos = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    zpos = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    xvi = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    yvi = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    zvi = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    forcesx = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    forcesy = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    forcesz = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    dpx = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    dpy = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);
    dpz = (DATATYPE *)_mm_malloc(N * sizeof(DATATYPE), ALIGNMENT);

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

        #pragma omp parallel for schedule(static) private(i,j)
        for(ii = 0; ii < N; ii+=BLOCKSIZE){
			
			__declspec(align(ALIGNMENT)) DATATYPE forcesx[BLOCKSIZE] = {0.0};
			__declspec(align(ALIGNMENT)) DATATYPE forcesy[BLOCKSIZE] = {0.0};
			__declspec(align(ALIGNMENT)) DATATYPE forcesz[BLOCKSIZE] = {0.0};

			#pragma unroll(BLOCKSIZE)
			for(j = 0; j < N; j++){
                PRINT2
                #pragma omp simd aligned(xpos,ypos,zpos,masses,forcesx,forcesy,forcesz)
                for (i = ii; i < ii+BLOCKSIZE; i++){
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
                    const DATATYPE d32 = 1/POW(dsquared,3.0f/2.0f);

                    // Unit vector x,y,z coordinates
                    // dx/d, dy/d, dz/d
                    // Force vector x,y,z coordinates
                    // F * 1/d² * [dx | dy | dz ] * 1/d
                    // F * dx * 1/d² * 1/d ==> F * dx * 1/d³
                    forcesx[i-ii] += F*dx*d32;
                    // F * dy * 1/d² * 1/d ==> F * dy * 1/d³
                    forcesy[i-ii] += F*dy*d32;
                    // F * dz * 1/d² * 1/d ==> F * dz * 1/d³
                    forcesz[i-ii] += F*dz*d32;
                }
            }
            #pragma omp simd aligned(dpx,dpy,dpz,xvi,yvi,zvi,masses,forcesx,forcesy,forcesz)
            for (i = ii; i < ii+BLOCKSIZE; i++){
              // Vector aceleración
              // |M = F * A| ---> |A = F / M|
              const DATATYPE ax = forcesx[i-ii] / masses[i];
              const DATATYPE ay = forcesy[i-ii] / masses[i];
              const DATATYPE az = forcesz[i-ii] / masses[i];
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
        }

  //      #pragma omp parallel for //simd aligned(xpos: 64, ypos:64, zpos:64)
		#pragma omp parallel
		{
			#pragma unroll(BLOCKSIZE)
			#pragma omp for simd aligned(xpos,ypos,zpos,dpx,dpy,dpz)
			for(int i = 0; i < N; i++){
			  // Update position
			  xpos[i] += dpx[i];
			  ypos[i] += dpy[i];
			  zpos[i] += dpz[i];
			}
		}
    }
	const double tw = dwalltime() - timetick;
    printf("%s;%d;%d;%d;%d;%lf;%lf;%lf\n",argv[0],N,D,num_threads,BLOCKSIZE, tw, ((double)D*(16*N*N+18*N))/(tw*1000000000), ((double)D*(20*N*N))/(tw*1000000000));

    _mm_free(masses);
    _mm_free(xpos);
    _mm_free(ypos);
    _mm_free(zpos);
    _mm_free(xvi);
    _mm_free(yvi);
    _mm_free(zvi);
    _mm_free(forcesx);
    _mm_free(forcesy);
    _mm_free(forcesz);
    _mm_free(dpx);
    _mm_free(dpy);
    _mm_free(dpz);
}
