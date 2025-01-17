
//Aurora Abbondanza
//Alessandro Agapito
//Andrea Belli Contarini
//Riccardo Caleno

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//STRUCT
typedef struct {
  double PE, KE, Q;
  double TotalM, Nout; //total mass and number of escapers
  //double Aipe; //Average initial particle energy for the normalization of the total energy
} param;

#define N 100 // Number of particles
#define tEnd 100.0 // Time at which simulation ends
#define dt 0.1 // Time-step
#define G 1.0 // Newton's Gravitational Constant (normalized)
#define RADIUS 10.0 // Sphere's radius
#define VOLUME ((4.0/3.0)*M_PI*RADIUS*RADIUS*RADIUS) //Sphere's volume
#define SOFTENING (RADIUS/(2.0*cbrt(N)))  // SOFTENING parameter length
#define Qini 0.001 //Initial value of Q


//EXTERNAL FUNCTIONS DEFINITION
void InitialPrints(void);
void InitialConditions(double pos[N][3],double vel[N][3], double mass[N], param* pointer);
void getAcc(double pos[N][3], double mass[], double acc[N][3]);
void getEnergy(double pos[N][3], double vel[][3], double mass[], param* pointer);
void getVirialValue(double pos[N][3], double vel[][3], double mass[], param* pointer);
void IntegrationFunction(double pos[N][3],double vel[][3], double acc[][3], double mass[N]);
void CountOutsider(double pos[N][3],double mass[N], param* pointer);


int main(){ /***********************************START MAIN********************************/
  //INITIAL CONDITIONS
  double mass[N] ;
  double pos[N][3];    
  double vel[N][3];
  double acc[N][3];
  int i, j, k, Nt;
  double t, E, E_ini_TOT, E_TOT_normalized, K_ini, PE_ini,KE_norm,PE_norm;
  
  param* pointer = (param*) malloc(sizeof(param)); //dynamic memory allocation
  pointer->KE=0.0;
  pointer->PE=0.0;
  pointer->TotalM=0.0;
  
  InitialPrints();
  
  //Masses' inizialization
  for (i=0; i<N; i++) {
    mass[i]=1.0;
    pointer->TotalM += mass[i];
  }
  
  //Initial conditions
  InitialConditions(pos, vel, mass, pointer);
    
  //Initial energy
  //getEnergy(pos, vel, mass, pointer);
  K_ini=pointer->KE;
  PE_ini=pointer->PE;
  E_ini_TOT=K_ini+PE_ini;

   //Opening files where we'll print the results
  FILE *file=fopen("positions.txt", "w");
  FILE *file2=fopen("energy.txt", "w");
  
  
  //Check correct files'opening
  if (file == NULL || file2 == NULL) {
    printf("Errore nell'apertura del file\n");
    return 1;
  }

   //Number of time-steps
  Nt = (int) ceil(tEnd/dt);
  t=0.0;

  fprintf(file2, "%g %g %g %g %g\n", t,(pointer->KE+pointer->PE)/(K_ini+PE_ini) , pointer->KE/K_ini, fabs((E_ini_TOT)/E_ini_TOT), (-2*pointer->KE/pointer->PE));
  

  //Initial average particle energy 
  //pointer->Aipe=E_ini_TOT/N;
    
 
  
  
 

 
  for(i=0; i<Nt; i++) {

     t+=dt;
     
    //Computing acceleration
    getAcc(pos, mass, acc);
    
    //Computing energy at every step (useful for further plotting)
    getEnergy(pos, vel, mass, pointer);
    E = pointer->KE+pointer->PE;

    //Computing Q taking account of the escapers
    getVirialValue(pos, vel, mass, pointer);

    //Integration funciton (LEAPFROG algorithm)
    IntegrationFunction(pos, vel, acc, mass);
    
    //Printing coordinates on file
    for(j=0; j<N; j++) {
      fprintf(file, " %g %g %g\n", pos[j][0],pos[j][1],pos[j][2]);
    }
    //fprintf(file,"\n");
    
    //Printing energies and virial ratio on file
    /*E_TOT_normalized = E/(-(double)(N-pointer->Nescapers)*pointer->Aipe);
    KE_norm = pointer->KE/((double)(N-pointer->Nescapers)*pointer->Aipe);
    PE_norm = pointer->PE/((double)(N-pointer->Nescapers)*pointer->Aipe);*/

    E_TOT_normalized = E/(-E_ini_TOT);
    KE_norm = pointer->KE/(-E_ini_TOT);
    PE_norm = pointer->PE/(-E_ini_TOT);
    
    fprintf(file2, "%g %g %g %g %g\n", t, E_TOT_normalized, KE_norm, fabs((E-E_ini_TOT)/E_ini_TOT), pointer->Q);

    //fprintf(file2, "%g\t %.4lf\t %.4lf\t\n", t, pointer->KE, pointer->PE);
    
   
  }
  CountOutsider(pos,mass,pointer);
  fclose(file2);
  fclose(file);
  
  free(pointer); //memory deallocation
} /***********************************END MAIN********************************/

//EXTERNAL FUNCTIONS
void InitialPrints(void){
  fprintf(stdout, "------------------------------------------------------\nThe following program will solve the N-Body problem\n------------------------------------------------------\nNumber of bodies N = %d \nSphere's RADIUS = %g  \nSphere's VOLUME = %g \nSOFTENING parameter's lenght = %g \n------------------------------------------------------\n", N, RADIUS, VOLUME, SOFTENING);
}

void InitialConditions(double pos[N][3],double vel[N][3], double mass[N], param* pointer){
  srand(time(NULL)); // Inizialization random numbers' generator
  
  double x, y, z, Vx, Vy, Vz, dx, dy, dz, r, V, sin_theta, phi,alpha, PEpart;
  int i, j;
  
  //******** POSITIONS'S INIZIALIZATION *************
  for (i=0; i <N; i++) {
    // Generiamo una distanza radiale casuale uniforme nell'intervallo [0,1] usando il metodo della trasformazione inversa
    r = cbrt((double) rand() /(double) RAND_MAX) * RADIUS;
    
    // Generiamo un angolo polare e uno azimutale casuale uniforme nell'intervallo [0, pi] e [0, 2*pi], rispettivamente
    sin_theta = ((double) rand() /(double) RAND_MAX);
    phi = ((double) rand() /(double) RAND_MAX) * 2 * M_PI;
    
    // Conversion from spherical to cartesian coordinates
    x = r * sin_theta * cos(phi);
    y = r * sin_theta * sin(phi);
    z = r * sqrt(1-sin_theta*sin_theta);
    
    //Salvo le nuove coordinate
    pos[i][0]=x;
    pos[i][1]=y;
    pos[i][2]=z;
  }

  pointer->KE = 0.0;
  pointer->PE = 0.0;
  
  //********  VELOCITY'S INIZIALIZATION *************
  for (i=0; i<N; i++) {
    //Zeroing potenzial energy
    PEpart=0.0;
    
      //Computing potential energy for every body
    for (j=0; j<N; j++) {
      
      //Computing potential energy
      dx = pos[j][0] - pos[i][0];
      dy = pos[j][1] - pos[i][1];
      dz = pos[j][2] - pos[i][2];
      r = sqrt(dx*dx + dy*dy + dz*dz+SOFTENING*SOFTENING);
      if (r > 0.0){

	//Compute PE of the single particle
        PEpart += -G * mass[i] * mass[j] / r;

	//Compute total PE for the Q
	pointer->PE += -G * mass[i] * mass[j] / r;
      }
      }
    
    //Virial Theorem
    //V = rand();
    V = sqrt(-2*(PEpart/mass[i]))*((double)rand()/RAND_MAX);

    //Compute total KE for the Q
    pointer->KE += 0.5 * mass[i] * V*V;
    
    
    //Generating new coordinates
    sin_theta = ((double) rand() /(double) RAND_MAX);
    phi = ((double) rand() /(double) RAND_MAX) * 2 * M_PI;
    
    //Converting velocities into components
    Vx = V * sin_theta * cos(phi);
    Vy = V * sin_theta * sin(phi);
    Vz = V * sqrt(1-sin_theta*sin_theta);
    
    //Saving velocities
    vel[i][0]=Vx;
    vel[i][1]=Vy;
    vel[i][2]=Vz;
  }

  //Compute Q for normalization
  pointer->Q=(2.0*pointer->KE/fabs(pointer->PE));
  alpha=sqrt(Qini/pointer->Q);

  pointer->KE=0.0;

  //Renormalization of velocities
  for (i=0; i<N; i++) {
    vel[i][0]*=alpha;
    vel[i][1]*=alpha;
    vel[i][2]*=alpha;
    V=sqrt( vel[i][0]* vel[i][0]+ vel[i][1]* vel[i][1]+ vel[i][2]* vel[i][2]);
    pointer->KE+= 0.5 * mass[i] * V*V;
  }

 
  pointer->Q=(2.0*pointer->KE/fabs(pointer->PE));

  printf("Il valore di Q iniziale è %g\n\n\n", pointer->Q);
  
}

void getAcc(double pos[N][3], double mass[], double acc[N][3]){
  double r, dx, dy, dz;
  int i, j;
  for(i=0; i<N; i++){
    acc[i][0] = 0.0;
    acc[i][1] = 0.0;
    acc[i][2] = 0.0;
    for(j = 0; j < N; j++){
      if(i != j) {
	dx = pos[j][0] - pos[i][0];
	dy = pos[j][1] - pos[i][1];
	dz = pos[j][2] - pos[i][2];
	r = sqrt(dx*dx + dy*dy + dz*dz+SOFTENING*SOFTENING);
	acc[i][0] += G * mass[j] * dx / (r*r*r);
	acc[i][1] += G * mass[j] * dy / (r*r*r);
	acc[i][2] += G * mass[j] * dz / (r*r*r);
      } 
    }
  }
}

void getEnergy(double pos[][3], double vel[][3], double mass[], param* pointer){
  double dx, dy, dz, r;
  int i, j;
  double average_r=0.0,d;
  
  
  pointer->PE = 0.0;
  pointer->KE = 0.0;


     
  for(i=0; i<N; i++) {
    
    pointer->KE += 0.5 * mass[i] * (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2]);
    
    for(j=i+1; j < N; j++) {
      dx = pos[j][0] - pos[i][0];
      dy = pos[j][1] - pos[i][1];
      dz = pos[j][2] - pos[i][2];
      r = sqrt(dx*dx + dy*dy + dz*dz + SOFTENING*SOFTENING);
      
      pointer->PE += -G  * mass[j] * mass[i] / r;
    }
  }
}
  

void getVirialValue(double pos[][3], double vel[][3], double mass[], param* pointer){
  double dx, dy, dz, r,K=0.0,P=0.0,k,p,COM[3]={0.0,0.0,0.0};
  int i, j;
  double average_r=0.0,d;
  
  
  pointer->Q = 0.0;
  

  //Compute average_r to count outsiders
  for (i=0; i< N; i++) {
    COM[0] += mass[i]*pos[i][0]/pointer->TotalM;
    COM[1] += mass[i]*pos[i][1]/pointer->TotalM;
    COM[2] += mass[i]*pos[i][2]/pointer->TotalM;
   }

  //distances from the COM
  for (i=0; i< N; i++) {
    d=sqrt((pos[i][0]-COM[0])*(pos[i][0]-COM[0])+(pos[i][1]-COM[1])*(pos[i][1]-COM[1])+(pos[i][2]-COM[2])*(pos[i][2]-COM[2]));
    average_r+=(d/(double)(N-(pointer->Nout)));
   }
  
  pointer->Nout=0.0;
     
  for(i=0; i<N; i++) {
    d=sqrt((pos[i][0]-COM[0])*(pos[i][0]-COM[0])+(pos[i][1]-COM[1])*(pos[i][1]-COM[1])+(pos[i][2]-COM[2])*(pos[i][2]-COM[2]));

 
    //Compute Energy of the single particle
    k = 0.5 * mass[i] * (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2]);
    p=0;
    
    for(j=i+1; j < N; j++) {
      dx = pos[j][0] - pos[i][0];
      dy = pos[j][1] - pos[i][1];
      dz = pos[j][2] - pos[i][2];
      r = sqrt(dx*dx + dy*dy + dz*dz + SOFTENING*SOFTENING);

      p += -G  * mass[j] * mass[i] / r;
    }
    if((k+p)<0 || d<1.5*average_r){

      //If not Escapers Compute energies
      K += k;
      P += p;
    }else{
      pointer->TotalM -= mass[i];
      //mass[i]=0.0;                   //**********uncomment if u want to campute better Q, comment if u want better Energies!****************
      pointer->Nout++;
    }  
  }
  pointer->Q = (2*K)/fabs(P);
}




void IntegrationFunction(double pos[][3],double vel[][3], double acc[][3], double mass[]){
  //LEAPFROG - Simulation Main Loop
  int j, i;
  for (i=0; i< N; i++) {
    for (j=0; j<3; j++) {
      vel[i][j] += acc[i][j] * dt/2.0;
    }
  }
  
  // Drift
  for (i = 0; i<N; i++) {
    for (j = 0; j < 3; j++) {
      pos[i][j] += vel[i][j] * dt;
 
    }
  }
  //Recompute velocity at the right time step 
  getAcc(pos, mass, acc);
  
  for (i=0; i< N; i++) {
    for (j=0; j<3; j++) {
      vel[i][j] += acc[i][j] * dt/2.0;
    }
  } 
}

void CountOutsider(double pos[N][3],double mass[N], param* pointer){
  int Nout=0,i;
  double average_r=0.0,d,COM[3]={0.0,0.0,0.0};

  //Compute average_r to count outsiders
  for (i=0; i< N; i++) {
    COM[0] += mass[i]*pos[i][0]/pointer->TotalM;
    COM[1] += mass[i]*pos[i][1]/pointer->TotalM;
    COM[2] += mass[i]*pos[i][2]/pointer->TotalM;
   }

  //distances from the COM
  for (i=0; i< N; i++) {
    d=sqrt((pos[i][0]-COM[0])*(pos[i][0]-COM[0])+(pos[i][1]-COM[1])*(pos[i][1]-COM[1])+(pos[i][2]-COM[2])*(pos[i][2]-COM[2]));
    average_r+=(d/(double)(N-Nout));
   
   

    //Count outsider as particle with r>3*average_r
    if(d>1.5*average_r){
      Nout++;
    }
  }

  fprintf(stdout, "------------------------------------------------------\nThe total number of Escapers is:%d\n------------------------------------------------------\n",Nout);
}
  

