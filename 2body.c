//Aurora Abbondanza
//Alessandro Agapito
//Andrea Belli Contarini
//Riccardo Caleno

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 2     // Number of particles
#define tEnd 100.0  // time at which simulation ends
#define dt 0.1  // timestep
#define softening 0.2  // softening length
#define G 1.0   // Newton's Gravitational Constant

void terminal_initial_prints(void){
    fprintf(stdout, "-------------------------------------------------\nSolving N-Body problem with %i bodies.\nDuration time of %.2lf and time step dt=%.2lf\nSoftening parameter lenght of %.2lf\n-------------------------------------------------\n", N, tEnd, dt, softening);
}

void getAcc(double pos[N][3], double mass[], double acc[N][3]){
  double r, dx, dy, dz;
  int i,j;
  for( i = 0; i < N; i++){
    acc[i][0] = 0.0;
    acc[i][1] = 0.0;
    acc[i][2] = 0.0;
    for(j = 0; j < N; j++){
      if(i != j) {
	dx = pos[j][0] - pos[i][0];
	dy = pos[j][1] - pos[i][1];
	dz = pos[j][2] - pos[i][2];
	r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
	acc[i][0] += G * mass[j] * dx / (r*r*r);
	acc[i][1] += G * mass[j] * dy / (r*r*r);
	acc[i][2] += G * mass[j] * dz / (r*r*r);
      } 
    }
  }
}

void getEnergy(double pos[][3], double vel[][3], double mass[], double KE, double PE) {
  double dx, dy, dz, inv_r;
  int i,j;
  KE = 0.0;
  PE = 0.0;
  for(i = 0; i < N; i++) {
    KE += 0.5 * mass[i] * (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2]);
    for( j = i+1; j < N; j++) {
      dx = pos[j][0] - pos[i][0];
      dy = pos[j][1] - pos[i][1];
      dz = pos[j][2] - pos[i][2];
      inv_r = sqrt(dx*dx + dy*dy + dz*dz);   //Includere il softening nel potenziale????
      if (inv_r > 0.0) {
	PE += -G * mass[i] * mass[j] / inv_r;
      }
    }
  }
}


void IntegrationFunction(double pos[][3],double vel[][3], double acc[][3]){
   // leap frog - Simulation Main Loop
  int j,i;
  for ( i = 0; i < N; i++) {
    for ( j = 0; j < 3; j++) {
      vel[i][j] += acc[i][j] * dt/2.0;
    }
  }
        
  // drift
  for (i = 0; i<N; i++) {
    for (j = 0; j < 3; j++) {
    pos[i][j] += vel[i][j] * dt;
    }
  }
}

       
     



/*void initialConditions(double pos[][], double vel[][])    *******INIZIALIZZAZIONE PROBLEMA********


  
 */

int main() {
  //INITIAL CONDITIONS
  double mass[N] = {1.0, 1.0};
  double pos[N][3] = {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}};
  double vel[N][3] = {{0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double acc[N][3] = {{0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  int i, j, k, Nt;
  double t,KE,PE;
   
  terminal_initial_prints();

  // set the random number generator seed
  // srand(17);
  
  /* // Convert to Center-of-Mass frame
  double vcm[3] = {0.0, 0.0, 0.0};
  for (i=0; i<N; i++) {
    for (j=0; j<3; j++) {
      vcm[j] += mass[i] * vel[i][j];
    }
  }
  for (j=0; j<3; j++) {
           vcm[j] /= (mass[0] + mass[1] + mass[2]);
       }
  for (i=0; i<N; i++) {
       for (j=0; j<3; j++) {
               vel[i][j] -= vcm[j];
       }
  }       TRASFERISCE IL PROBLEMA NEL SDR DEL CENTRO DI MASSA,INUTILE SE INIZIALIZZIAMO BENE LE VELOCITà INIZIALI*/     

       // calculate initial gravitational accelerations
   // getAcc(pos, mass, acc);

       // calculate initial energy of system
   // getEnergy(pos, vel, mass, KE, PE);

       // number of timesteps
       Nt = (int)ceil(tEnd/dt);
       t=0;
    
        // apertura file per OUTPUT 3D
      //FILE *file=fopen("output3D.txt","w");
      
       // apertura file per OUTPUT 2D
       FILE *file=fopen("output3D.txt","w");

       //Verifica apertura file
       if (file == NULL) {
        printf("Errore nell'apertura del file\n");
        return 1;
       }
       
       for(i=0; i<Nt; i++) {
         getAcc(pos, mass, acc);             
	 getEnergy(pos, vel, mass, KE, PE);   //Calcolo energia ad ogni step utile se vogliamo plottare 
         IntegrationFunction(pos, vel,  acc);

	 for(j=0; j<N; j++) {               //Ciclo per printare su file
	   //3D
        fprintf(file, "%g %g %g ", pos[j][0],pos[j][1],pos[j][2]);
        //2D
        //fprintf(file, "%g %g %d ", pos[j][0],pos[j][1], i);
	 }
	 fprintf(file,"\n");             
	 t+=dt;
       }
       fclose(file);
       return 0;
}

