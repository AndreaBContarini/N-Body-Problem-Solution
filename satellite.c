//Aurora Abbondanza
//Alessandro Agapito
//Andrea Belli Contarini
//Riccardo Caleno

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//STRUCT
typedef struct {
  int i;
  double x; double y; double z; double vx; double vy; double vz; double m_dot_ov_m;
} param;

//Constants definition
#define m_sat 11110 //[Kg] satellite mass (Hubble)
#define r_sat 4.2 //[m] satellite radius (Hubble)
#define M_e 5.972e24 //[Kg] Earth mass
#define R_e 6.364e6 //[m] Earth radius
#define G 6.674e-11 //[N*m^2*Kg^-1] Gravitational constant
#define r_geo 4.216e7 //[m] geostationary orbit
#define m (m_sat/M_e) //new useful variable for mass, adimensional
#define r_STOP (R_e/r_geo) //adimensional Earth's radius

#define lambda2_ov_m 1e-5 //[s^-1] atmospheric drag per unit mass on the initial mass m(0)
#define m_dot_ov_m_initial -1e-4 // adimensiona relative mass loss (dm'/dt')/m' respect to the initial mass m(0)
#define beta1 ((lambda2_ov_m + m_dot_ov_m_initial) * sqrt((r_geo*r_geo*r_geo) / (G*M_e))) //adimensional drag_coefficient

#define dt 1e-3 //size of the time step for RK4
#define T_circular_teo (1./(2*(lambda2_ov_m + m_dot_ov_m_initial))) * log(r_geo/R_e) //[s]
# define T_radial_teo (2./(3*sqrt(G*M_e)))*(pow(r_geo,(3./2))-pow(R_e,(3./2)))//[s] In the case of v(0)=-1

//Functions' definition
int problem_choice(void);
void terminal_initial_prints_circular(param* pointer);
void terminal_initial_prints_radial(param* pointer);
void terminal_final_prints_circular(param* pointer);
void terminal_final_prints_radial(param* pointer);
void RK4_CIRCULAR(param* pointer);
void RK4_RADIAL(param* pointer);
double a1(double x, double vx, int i, double m_dot_ov_m);
double ax(double x, double y, double z, double vx, int i, double m_dot_ov_m);
double ay(double x, double y, double z, double vy, int i, double m_dot_ov_m);
double az(double x, double y, double z, double vz, int i, double m_dot_ov_m);
double r_module(param* pointer);
double r_t(param* pointer);
double v_module2(param* pointer);
double energy(param* pointer);
double energy_t(param* pointer);

int main(){ /*-----------------------MAIN STARTS------------------------*/
  int input=problem_choice(), superata=0;
  FILE *file1;
  FILE *file2;
  FILE *file3;
  FILE *file4;
  FILE *file5;
  param* pointer = (param*) malloc(sizeof(param)); //dynamic memory allocation
  
  if (input==1) {  /*******************************CIRCULAR MOTION 3D**************************************/
    
    //INITIAL CONDITIONS
    /* 3 dimensions */
    /*
      pointer->x=1.; pointer->vx=0.; //Note: x = (r_0/r_geo) == 1
      pointer->y=0.; pointer->vy=1./sqrt(2);
      pointer->z=0.; pointer->vz=1./sqrt(2);
      pointer->i=0; */
    
    /* 2 dimensions */
    /* */
    pointer->x=1.; pointer->vx=0.; //x=r_0/r_geo==1
    pointer->y=0.; pointer->vy=0.;
    pointer->z=0.; pointer->vz=1.;
    pointer->i=0;
    
    double ris = r_module(pointer);
    terminal_initial_prints_circular(pointer);
    
    //File opening, where we'll print the results
    file1 = fopen("Circular2d.txt", "w");
    file2 = fopen("Rvst.txt", "w");
    file3 = fopen("Evst.txt", "w");
    
    // Check if the files have been correctly open
    if ((file1 == NULL)||(file2 == NULL)||(file3 == NULL)) {
      fprintf(stderr, "Error in files opening!\n");
      return 1;
    }
    
    //Print on file at t=0 (initial conditions)
    fprintf(file1, "%g %g\n", pointer->x, pointer->z );
    //fprintf(file1, "%g %g %g\n", pointer->x ,pointer->y, pointer->z );    fprintf(file2, "%g %g\n", ris, pointer->i*dt);
    fprintf(file3, "%g %g %g\n", energy(pointer), energy_t(pointer), (pointer->i)*dt);
    
    do { pointer->i++;
      double tempo = (pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e))/3600;
      if (tempo > 10.) {
	pointer->m_dot_ov_m=m_dot_ov_m_initial;
      } else {
	pointer->m_dot_ov_m=0.;
      }
      
      //Computing: Runge-Kutta 4
      RK4_CIRCULAR(pointer);
      
      //Computing: radius
      ris = r_module(pointer);
      
      //Printing coordinates if the satellite has not landed yet
      if (r_module(pointer)>r_STOP) {
	
	fprintf(file1, "%g %g\n", pointer->x, pointer->z);
	if (v_module2(pointer)>((2*r_geo)/ris) && superata==0) {
	  superata++;
	  fprintf(stdout, "Superata velocità di fuga v_escape = %g; tale velocità è v = %g \nIn: r = %g (x = %g; z = %g)\nTEMPO: t=%g",sqrt((2*r_geo)/ris), (sqrt(pointer->vx*pointer->vx  + pointer->vz*pointer->vz)), (sqrt(pointer->x*pointer->x  + pointer->z*pointer->z)), pointer->x, pointer->z, (pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e))/3600);
	}
	//fprintf(file1, "%g %g %g\n", pointer->x ,pointer->y, pointer->z );
	//fprintf(file2, "%g %g\n",ris, pointer->i*dt);
	//fprintf(file3, "%g %g %g\n", energy(pointer), energy_t(pointer), pointer->i*dt);
      } 
    } while(r_module(pointer)>r_STOP && pointer->i<1000000);
    
    fclose(file1);
    fclose(file2);
    fclose(file3);
    terminal_final_prints_circular(pointer);
    
  } else {  /*******************************RADIAL MOTION 3D**************************************/
    
    //INITIAL CONDITIONS
    pointer->x=1.; pointer->vx=1.;
    
    terminal_initial_prints_radial(pointer);
    file4 = fopen("xVSt4.txt", "w");
    file5 = fopen("velocity.txt", "w");
    
    // Check if the files have been correctly open
    if ((file4==NULL)||(file2 == NULL)) {
      fprintf(stderr, "Error in file opening!\n");
      return 1;
    }
    
    //Print on file at t=0 (initial conditions)
    fprintf(file4, "%g %g\n", pointer->i*dt, pointer->x);
    
    do {
      pointer->i++;
      
      //Computing: Runge-Kutta 4
      RK4_RADIAL(pointer);
      
      //Printing coordinates if the satellite has not landed yet
      if (pointer->x > r_STOP) {
        fprintf(file4, "%g %g\n", pointer->x, pointer->i*dt);
      }
      
      //Printing velocities to see when the satellite changes direction of motion
      fprintf(file5, "%g\t %g\t\n", (pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e))/3600, pointer->vx);
    } while (pointer->x > r_STOP);
    
    fclose(file4);
    fclose(file5);
    terminal_final_prints_radial(pointer);
  }
  free(pointer); //memory deallocation
} /*-----------------------------MAIN ENDS-----------------------------*/

/*----------------------- Functions' structures ----------------------- */
int problem_choice(void){
  int input;
  printf("------------------------------------------------------\nThe following program will solve the satellite problem\n------------------------------------------------------\nThe parameters of the problem are:\n drag = %g Hz\n beta1 = %g\n dt = %g \n (dm/dt)*1/m = %g Hz\n------------------------------------------------------\nType 1 to solve, starting with a circular motion.\nType 2 to solve, starting with a radial motion.\n",lambda2_ov_m , beta1, dt , m_dot_ov_m_initial );
  scanf("%i", &input);
  while (input != 1 && input != 2){
    printf("ERROR! Must type 1 or 2. Type again: ");
    scanf("%i", &input);
  }
  return(input);
}

void terminal_initial_prints_circular(param* pointer){
  fprintf(stdout, "------------------------------------------------------------\nCIRCULAR MOTION\n-Printing coordinates on file 'Circular.txt' [1:x  2:y  3:z]\nWith adimensional INITIAL CONDITIONS:\nx=%g, vx=%g;  y=%g, vy=%g;  z=%g, vz=%g;\n------------------------------------------------------------\n-Printing Radius and time on file 'Rvst.txt' [1:r, 2:t]\n------------------------------------------------------------\n-Printing Energy and time on file 'Evst.txt' [1:E, 2:t]\n", pointer->x,pointer->vx,pointer->y,pointer->vy,pointer->z,pointer->vz);
}

void terminal_initial_prints_radial(param* pointer){
  fprintf(stdout, "------------------------------------------------------------\nRADIAL MOTION\n-Printing radial coordinate on file 'Circular.txt' [1:r]\nWith adimensional INITIAL CONDITIONS:\nr=%g, vr=%g;\n------------------------------------------------------------\n-Printing Radius and time on file 'Rvst.txt' [1:r, 2:t]\n------------------------------------------------------------\n-Printing Energy and time on file 'Evst.txt' [1:E, 2:t]\n", pointer->x,pointer->vx);
}

void terminal_final_prints_circular(param* pointer){
  if(m_dot_ov_m_initial == 0){
    fprintf(stdout, "------------------------------------------------------------\nThe satellite has crushed on Earth!\nThe numerical impact time is: t_f = %g h\nThe theoretical impact time is: t_teo = %g h\n", (pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e))/3600, T_circular_teo/3600);
    
  } else {
    fprintf(stdout, "------------------------------------------------------------\nThe satellite has crushed on Earth!\nThe numerical impact time is: t_f = %g h\nThe relative loss of mass is: dm/m = %g\n", (pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e))/3600, -m_dot_ov_m_initial*(pointer->i*dt));
  }
}

void terminal_final_prints_radial(param* pointer){
  if(m_dot_ov_m_initial == 0){
    fprintf(stdout, "------------------------------------------------------------\nThe satellite has crushed on Earth!\nThe numerical impact time is: t_f = %g h\nThe theoretical impact time is: t_teo = %g h\n", ((pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e)))/3600, T_radial_teo/3600);
  } else {
    fprintf(stdout, "------------------------------------------------------------\nThe satellite has crushed on Earth!\nThe numerical impact time is: t_f = %g h\n", ((pointer->i*dt)*sqrt((r_geo*r_geo*r_geo)/(G*M_e)))/3600);
  } 
}

void RK4_RADIAL(param* pointer){
  double dx1, dx2, dx3, dx4, dv1, dv2, dv3, dv4;
  
  dv1 = a1(pointer->x, pointer->vx, pointer->i, pointer->m_dot_ov_m) * dt;
  dx1 = pointer->vx * dt;
  
  dx2 = (pointer->vx + 0.5 * dv1) * dt;
  dv2 = a1(pointer->x+0.5*dx1, pointer->vx + 0.5 * dv1, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dx3 = (pointer->vx + 0.5 * dv2) * dt;
  dv3 = a1(pointer->x+0.5*dx2, pointer->vx+0.5*dv2, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dx4 = (pointer->vx + dv3) * dt;
  dv4 = a1(pointer->x+dx3, pointer->vx+dv3, pointer->i, pointer->m_dot_ov_m) * dt;
  
  pointer->vx += ((dv1 + 2*dv2 + 2*dv3 + dv4)/6);
  pointer->x += ((dx1 + 2*dx2 + 2*dx3 + dx4)/6);
}

void RK4_CIRCULAR(param* pointer){
  double dx1, dx2, dx3, dx4, dvx1, dvx2, dvx3, dvx4;
  double dy1, dy2, dy3, dy4, dvy1, dvy2, dvy3, dvy4;
  double dz1, dz2, dz3, dz4, dvz1, dvz2, dvz3, dvz4;
  
  /*----- 1st step -----*/
  dvx1 = ax(pointer->x, pointer->y, pointer->z, pointer->vx, pointer->i, pointer->m_dot_ov_m) * dt;
  dx1 = pointer->vx * dt;
  
  dvy1 = ay(pointer->x, pointer->y, pointer->z, pointer->vy, pointer->i, pointer->m_dot_ov_m) * dt;
  dy1 = pointer->vy * dt;
  
  dvz1 = az(pointer->x, pointer->y, pointer->z, pointer->vz, pointer->i, pointer->m_dot_ov_m) * dt;
  dz1 = pointer->vz * dt;
  
  /*----- 2nd step -----*/
  dx2 = (pointer->vx + 0.5 * dvx1) * dt;
  dvx2 = ax(pointer->x+0.5*dx1, pointer->y+0.5*dy1, pointer->z+0.5*dz1, pointer->vx + 0.5 * dvx1, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dy2 = (pointer->vy + 0.5 * dvy1) * dt;
  dvy2 = ay(pointer->x+0.5*dx1, pointer->y+0.5*dy1, pointer->z+0.5*dz1, pointer->vy + 0.5 * dvy1, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dz2 = (pointer->vz + 0.5 * dvz1) * dt;
  dvz2 = az(pointer->x+0.5*dx1, pointer->y+0.5*dy1, pointer->z+0.5*dz1, pointer->vz + 0.5 * dvz1, pointer->i, pointer->m_dot_ov_m) * dt;
  
  /*----- 3rd step -----*/
  dx3 = (pointer->vx + 0.5 * dvx2) * dt;
  dvx3 = ax(pointer->x+0.5*dx2, pointer->y+0.5*dy2, pointer->z+0.5*dz2, pointer->vx+0.5*dvx2, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dy3 = (pointer->vy + 0.5 * dvy2) * dt;
  dvy3 = ay(pointer->x+0.5*dx2, pointer->y+0.5*dy2, pointer->z+0.5*dz2, pointer->vy+0.5*dvy2, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dz3 = (pointer->vz + 0.5 * dvz2) * dt;
  dvz3 = az(pointer->x+0.5*dx2, pointer->y+0.5*dy2, pointer->z+0.5*dz2, pointer->vz+0.5*dvz2, pointer->i, pointer->m_dot_ov_m) * dt;
  
  /*----- 4th step -----*/
  dx4 = (pointer->vx + dvx3) * dt;
  dvx4 = ax(pointer->x+dx3, pointer->y+dy3, pointer->z+dz3, pointer->vx+dvx3, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dy4 = (pointer->vy + dvy3) * dt;
  dvy4 = ay(pointer->x+dx3, pointer->y+dy3, pointer->z+dz3, pointer->vy+dvy3, pointer->i, pointer->m_dot_ov_m) * dt;
  
  dz4 = (pointer->vz + dvz3) * dt;
  dvz4 = az(pointer->x+dx3, pointer->y+dy3, pointer->z+dz3, pointer->vz+dvz3, pointer->i, pointer->m_dot_ov_m) * dt;
  
  /*----- RESULTS -----*/
  pointer->vx += ((dvx1 + 2*dvx2 + 2*dvx3 + dvx4)/6);
  pointer->x += ((dx1 + 2*dx2 + 2*dx3 + dx4)/6);
  
  pointer->vy += ((dvy1 + 2*dvy2 + 2*dvy3 + dvy4)/6);
  pointer->y += ((dy1 + 2*dy2 + 2*dy3 + dy4)/6);
  
  pointer->vz += ((dvz1 + 2*dvz2 + 2*dvz3 + dvz4)/6);
  pointer->z += ((dz1 + 2*dz2 + 2*dz3 + dz4)/6);
}

double a1(double x, double vx, int i, double m_dot_ov_m){
  double acceleration;
  double beta=((lambda2_ov_m + m_dot_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)));
  if(i*dt < (-0.1/m_dot_ov_m)){
    acceleration= (-1./(x*x)) - (beta/(1+m_dot_ov_m*(i*dt)))*vx;
    return (acceleration);
  } else {
    acceleration = (-1./(x*x)) - (lambda2_ov_m)*vx;
    return (acceleration);
  }
}

double ax(double x, double y, double z, double vx, int i, double m_dot_ov_m){
  double a;
  double beta=((lambda2_ov_m + m_dot_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)));
  if(i*dt < (-0.1/m_dot_ov_m)){
    //choice of mass'loss starting from the initial instant
    a = -(pow((x*x + y*y + z*z),-3./2))*x - (beta/(1+m_dot_ov_m*(i*dt)))*vx;
    return (a);
  } else {
    a = -(pow((x*x + y*y + z*z),-3./2))*x - ((lambda2_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)))*vx ;  //maximum loss of 10%
    return (a);
  }
}

double ay(double x, double y, double z, double vy, int i, double m_dot_ov_m){
  double a;
  double beta=((lambda2_ov_m + m_dot_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)));
  if(i*dt < (-0.1/m_dot_ov_m)){
    a = -(pow((x*x + y*y + z*z),-3./2))*y - (beta/(1+m_dot_ov_m*(i*dt)))*vy;
    return (a);
  } else {
    a = -(pow((x*x + y*y + z*z),-3./2))*y - ((lambda2_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)))*vy ;
    return (a);
  }
}

double az(double x, double y, double z, double vz, int i, double m_dot_ov_m){
  double a;
  double beta=((lambda2_ov_m + m_dot_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)));
  if(i*dt < (-0.1/m_dot_ov_m)){
    a = -(pow((x*x + y*y + z*z),-3./2))*z - (beta/(1+m_dot_ov_m*(i*dt)))*vz;
    return (a);
  } else {
    a = -(pow((x*x + y*y + z*z),-3./2))*z - ((lambda2_ov_m) * sqrt((r_geo*r_geo*r_geo) / (G*M_e)))*vz ;
    return (a);
  }
}

double r_module(param* pointer){
  double r = sqrt((pointer->x*pointer->x)+(pointer->y*pointer->y)+(pointer->z*pointer->z));
  return(r);
}

double v_module2(param* pointer){
  double v = ((pointer->vx*pointer->vx)+(pointer->vy*pointer->vy)+(pointer->vz*pointer->vz));
  return(v);
}

double r_t(param* pointer){
  if(pointer->i*dt < (-0.1/m_dot_ov_m_initial)){
    double r_t = exp(-2*(beta1/(1+m_dot_ov_m_initial*(pointer->i*dt)))*(pointer->i*dt));
    return(r_t);
  } else {
    double r_t = exp(-2*((beta1/(1+m_dot_ov_m_initial))*pointer->i*dt));
    return(r_t);
  }	   
}

double energy(param* pointer){ 
  if(pointer->i*dt < (-0.1/m_dot_ov_m_initial)){
    double E = (0.5 * v_module2(pointer) - (1./r_module(pointer)))/(1+m_dot_ov_m_initial*(pointer->i*dt));
    return(E);
  } else {
    double E = (0.5 * v_module2(pointer) - (1./r_module(pointer)))/(1+m_dot_ov_m_initial);
    return(E);
  }
}

double energy_t(param* pointer){
  if(pointer->i*dt < (-0.1/m_dot_ov_m_initial)){
    double E = (- 0.5 * (1./r_module(pointer)))/(1+m_dot_ov_m_initial*(pointer->i*dt));
    return(E);
  } else {
    double E = (- 0.5 * (1./r_module(pointer)))/(1+m_dot_ov_m_initial);
    return(E);
  }
}
