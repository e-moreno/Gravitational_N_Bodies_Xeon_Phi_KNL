#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>


//Constante universal de gravitación
const double G = 6.674e-11;
//Diferencial de tiempo 1 segundo
const double DT = 1;

int N;
int iterations = 100;
int *data;
double *xpos, *ypos, *zpos;
double *xv, *yv, *zv;
double *mass;
double timetick;
double *forcesx;
double *forcesy;
double *forcesz;
int *ids;
double EPS = 3E4;
//Para calcular tiempo
double dwalltime(){
        double sec;
        struct timeval tv;

        gettimeofday(&tv,NULL);
        sec = tv.tv_sec + tv.tv_usec/1000000.0;
        return sec;
}


int main(int argc, char *argv[]){
	if(argc == 3){
		N = atoi(argv[1]);
		iterations = atoi(argv[2]);
	}else{
		printf("Los parámetros del programa deben ser %s N Iterations\n",argv[0]);
		return 0;
	}
	printf("%d, %d\n", N, iterations);
	xpos = (double *)malloc(N*sizeof(double));
	ypos = (double *)malloc(N*sizeof(double));
  zpos = (double *)malloc(N*sizeof(double));
	xv = (double *)malloc(N*sizeof(double));
	yv = (double *)malloc(N*sizeof(double));
  zv = (double *)malloc(N*sizeof(double));
	mass = (double *)malloc(N*sizeof(double));
	forcesx = (double *)malloc(N*N*sizeof(double));
	forcesy = (double *)malloc(N*N*sizeof(double));
  forcesz = (double *)malloc(N*N*sizeof(double));
	int sqrt_n = sqrt(N);
	if(!((sqrt_n * sqrt_n) == N)){
		sqrt_n++;
	}
	//Distancia en cada eje entre los cuerpos
	int dist = 1000000;
	//Inicializo datos de los cuerpos
	for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      forcesx[i*N+j] = 0.0;
      forcesy[i*N+j] = 0.0;
      forcesz[i*N+j] = 0.0;
    }
		int posx = (i/sqrt_n),posy = (i%sqrt_n);
    xpos[i] = (i%sqrt_n)*dist;
		ypos[i] = dist * i;
    zpos[i] = 5000;
		xv[i] = 0.0;
		yv[i] = 0.0;
    zv[i] = 0.0;
    mass[i] = 5.97e25;
    if(i == 6){
      mass[i] = 5.97e25;
    }
	}
	timetick = dwalltime();

  double distance,fmagnitude;
	double xdistance,ydistance,zdistance;
	double xvers, yvers, zvers;
	double ax,ay,az,dvx,dvy,dvz,dpx,dpy,dpz,dax,day,daz;
	double fx,fy,fz;
  for(int it = 0; it < iterations; it++){
  	for(int i = 0; i < N; i ++){
		printf("%.0lf, %.0lf, %.0lf\n", xpos[i], ypos[i], zpos[i]);
  		for(int j = i+1; j < N; j++){
        //Vector de distancia con dirección Cuerpo_i ---> Cuerpo_j
        //Para obtener la dirección de la fuerza ejercida sobre Cuerpo_i por Cuerpo_j
  			xdistance = xpos[j] - xpos[i];
  			ydistance = ypos[j] - ypos[i];
        zdistance = zpos[j] - zpos[i];
  			//Magnitud de la distancia entre ambos cuerpos
  			distance = sqrt((xdistance*xdistance) + (ydistance*ydistance) + (zdistance*zdistance));
  			//Magnitud de la fuerza de atracción gravitacional
  			//Usando la ley de gravitación de Newton
  			// F = G * ((m1 * m2) / D²)
    			fmagnitude = G * ((mass[i] * mass[j]) / (distance * distance));
    			//Componentes del vector unidad, para obtener la dirección de la fuerza
    			xvers = xdistance/distance;
    			yvers = ydistance/distance;
          zvers = zdistance/distance;
          //Ignorar colisiones
          if(distance < EPS*2){
            //Componentes del vector fuerza ejercida por j sobre i.
      			forcesx[i*N+j] = 0.00;
      			forcesy[i*N+j] = 0.00;
            forcesz[i*N+j] = 0.00;
      			//Fuerzas opuestas (simétrica), ejercida por i sobre j, multiplico las componentes por -1.
      			forcesx[j*N+i] = 0.00;
      			forcesy[j*N+i] = 0.00;
            forcesz[j*N+i] = 0.00;
          }else{
      			//Componentes del vector fuerza ejercida por j sobre i.
      			forcesx[i*N+j] = (xvers)*fmagnitude;
      			forcesy[i*N+j] = (yvers)*fmagnitude;
            forcesz[i*N+j] = (zvers)*fmagnitude;
      			//Fuerzas opuestas (simétrica), ejercida por i sobre j, multiplico las componentes por -1.
      			forcesx[j*N+i] = ((-1)*xvers)*fmagnitude;
      			forcesy[j*N+i] = ((-1)*yvers)*fmagnitude;
            forcesz[j*N+i] = ((-1)*zvers)*fmagnitude;
          }
  		}
  	}
  	for(int i = 0; i < N; i++){
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
  		//Vector aceleración
  		//Cálculo usando |M = F * A| ---> |A = F / M|
  		for(int b = 0; b < N; b++){
  			if(b != i){
  				fx += forcesx[i*N+b];
  				fy += forcesy[i*N+b];
          fz += forcesz[i*N+b];
  			}
  		}
      //Diferencia de aceleración
  		ax = fx / mass[i] * DT;
  		ay = fy / mass[i] * DT;
      az = fz / mass[i] * DT;

  		//Diferencia de velocidad delta v
  		dvx = ax * DT;
  		dvy = ay * DT;
      dvz = az * DT;
  		//Diferencia de posición
  		//Leapfrog scheme (velocidad inicial + (1/2 * nueva velocidad)) * DT
  		dpx = ( xv[i] + (dvx / 2) ) * DT;
  		dpy = ( yv[i] + (dvy / 2) ) * DT;
      dpz = ( zv[i] + (dvz / 2) ) * DT;

  		//Actualizo posición
  		xpos[i] += dpx;
  		ypos[i] += dpy;
      zpos[i] += dpz;
  		//Actualizo velocidad
  		xv[i] = (xv[i] + dvx)/2;
  		yv[i] = (yv[i] + dvy)/2;
      zv[i] = (zv[i] + dvz)/2;

  	}
  }

	printf("Tiempo en segundos %f\n", dwalltime() - timetick);
  free(xpos);
	free(ypos);
  free(zpos);
	free(xv);
	free(yv);
  free(zv);
	free(mass);
	free(forcesx);
	free(forcesy);
  free(forcesz);
	return 0;
}
