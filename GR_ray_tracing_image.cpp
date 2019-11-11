#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>

using namespace std;
using std::ifstream;
#define PI 3.14159265358979323846
#define N 6

#define t0 0   //Always start with at time t0=0

//save constants of motion globally 
double L, energy1;
double kappa;

//Prepare to read initail conditions in from file IP.txt
ifstream inFile;
#include <iostream>
#include <fstream>

std::ifstream input("IP.txt"); //Open IP.txt containing initial conditions
double data[7];  //define 7-vector to save the initial conditions 




void initial(double *y0, double *ydot0, double x, double y)
{
      for (int i = 0; i < 6; i++) {
           input >> data[i];      //save initial conditions in data 

           }
      //redefine the inital conditions
      double r0=data[0];     //radius
      double theta0=data[1]; //polar angle
      double phi0=data[2];  // azimuthal angle
      double a=data[3];   //spin
      double a2=a*a;
      double RHor=data[4];  //ISCO, which depend on spin.


      //Save the initial conditions in a vector y0
      y0[0]=r0;                   
      y0[1]=theta0; 
      y0[2]=phi0; 

      // Define the initial velocity (see thesis).
      double rdot0=-cos(y)*cos(x); //radial 
      double thetadot0=sin(y)/r0;  // polar
      double phidot0= cos(y)*sin(x)/(r0*sin(theta0)); //azimuthal

      //Define some parameters which will be used later on:
      double r2=r0*r0;
      double costheta2=cos(theta0)*cos(theta0);
      double sintheta2=sin(theta0)*sin(theta0);
      double sum=r2+a2*costheta2;
      double delta=r2-2.0*r0+a2;

      // Use the initial velocity to define the momentum vector
      y0[4]= rdot0*sum/delta;  //define pr
      y0[5]= thetadot0*sum;    //define ptheta

      // save the four velocity in vector ydot0
      ydot0[0] = rdot0;        //rdot
      ydot0[1] = thetadot0;    //thetadot
      ydot0[2] = cos(y)*sin(x)/(r0*sin(theta0)); //phidot
      
      //Calculate the two constants of motion:
      double energy2=(sum-2.0*r0)*(rdot0*rdot0/delta+thetadot0*thetadot0)+delta*sintheta2*phidot0*phidot0; // Calculate the energy squared (see thesis)
      double energy = sqrt(energy2);
      energy1=energy;
      L=(sum*delta*phidot0-2.0*a*r0*energy)*sintheta2/(sum-2.0*r0); //angluar momentum in phi


      //Normalize four momentum
      y0[4]=y0[4]/energy;
      y0[5]=y0[5]/energy;  
      L=L/energy; 

      //Define kappa, which is used to evolve trajectory (see thesis)
      kappa=y0[5]*y0[5]+a2*sintheta2+L*L/sintheta2;
} 


// setup the 6 parameters for the Fouth order Runge-Kutta integrator
void geodesic(double *y, double *dydlamda) 
{

  // Read data from the inital conditions vector (data)
  for (int i = 0; i < 6; i++) {
       input >> data[i];               
       }
  //Define parameters, so it is easier to remember, which is which.
  double a=data[3];
  double a2=a*a;
  double r, theta, phi, t, pr, ptheta;
  r=y[0];
  theta= y[1];
  phi = y[2];
  t = y[3];
  pr = y[4];
  ptheta = y[5];
  double r2=r*r;
  double twor=2.0*r;
  double cos2=cos(theta)*cos(theta);
  double sin2=sin(theta)*sin(theta);
  double sum=r2+a2*cos2;
  double delta=r2-twor+a2;
  double sd=sum*delta;


  // Calculate the four velocity (see thesis)
  dydlamda[0] = pr*delta/sum;  //dot r
  dydlamda[1] = ptheta/sum;    //dot theta
  dydlamda[2] = (twor*a+(sum-twor)*L/sin2)/sd;   // dot phi
  dydlamda[3] = 1.0+(twor*(r2+a2)-twor*a*L)/sd; // dot t and use normalization so E=1

  // Calculate the derviative of the four momentum (see thesis)
  dydlamda[4] = ((r-1.0)*(-kappa)+twor*(r2+a2)-2.0*a*L)/sd-2.0*pr*pr*(r-1.0)/sum; // dot p_r and use H=0 for photons 
  dydlamda[5] = sin(theta)*cos(theta)*(L*L/(sin2*sin2)-a2)/sum; // dot p_theta use H=0, 

}


// Evolve the 6 parameters with a standard fourth order Runge-kutta integrator. Here the paramers are found in the book Numerical Recipes in C by William Press.

void rkck (double *y, double *dydx, double h, double *yout, double *yerr)
{

  // Constants from  Cash and Karp in the book Numerical Recipes in C by William Press.
  static const double b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0, b41 = 0.3, b42 = -0.9,
    b43 = 1.2, b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0,
    b54 = 35.0/27.0, b61 = 1631.0/55296.0, b62 = 175.0/512.0,
    b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
    c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0,
    dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25; 

  
  int i;
  int n=N;

  double ak2[6], ak3[6], ak4[6], ak5[6], ak6[6], ytemp[6];


  for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i];  //first step
  geodesic(ytemp, ak2); // update integration and get derivatives in ak2.
  
  
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]); //second step: use result from derivative ak2 to get better estimate ak3
  geodesic(ytemp, ak3);
 
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]); // third step
  geodesic(ytemp, ak4);

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);  //fouth step
  geodesic(ytemp, ak5);  

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]); // fith step can be used to adjust the step-size so we can setup adaptive step-size
  geodesic(ytemp, ak6); 

  for (i=0;i<n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); // accumalete incremenets with proper weights. (using fourth order)
 
  for (i=0;i<n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]); // calculate/estimate the error between 5th and 4th order to see if it is below a treshhold and smaller stepsize are needed

}


// Define functions max and min to find which of two parameters are larger/smaller. This will be used later to define adaptive step-sizes.
double MAX(double x, double y)
{
    if (x >= y)
      return x;
    else 
      return y;
}

double MIN(double x, double y)
{
    if (x <= y)
      return x;
    else
      return y;
}


// Give the following error message if the system cannot update the null geodesics
void nrerror(const string error_text)
{
  cerr << "Numerical Recipes run-time error..." << endl;
  cerr << error_text << endl;
  cerr << "...now exiting to system..." << endl;
  exit(1);
}

#include <cmath>


// This give us adaptive size for the integrator. 
void rkqs (double *y, double *dydx, double &x, double htry, double eps, double *yscal, double &hdid, double &hnext)
{
  const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
  
  int i;

  double errmax, h, htemp, xnew;

  int n=N;
  h=htry;
  double yerr[6], ytemp[6];


  for (;;) 
    {   
      rkck(y, dydx, h, ytemp, yerr);

      errmax=0.0;
      for(i=0;i<n;i++) errmax = MAX(errmax, fabs(yerr[i]/yscal[i])); //errmax is the error tolotance and fabs(yerr[i]/yscal[i]) is the relative error, which ever is bigger will be errmax

      errmax/=eps;  
      if(errmax<=1.0) break; //if errmax is less than our treshhold, then keep

      // Shrink step-size 
      htemp=SAFETY*h*pow(errmax, PSHRNK);  // Shrinks stepsize h with by a factor (1/err)^(1/5)
      h=(h>=0.0 ? MAX(htemp,0.1*h) : MIN(htemp, 0.1*h)); //cannot vary stepsize to rapidly
    
      xnew=x+h; // check if stepsize is updated
      if(xnew==x) nrerror("stepsize underflow in rkqs"); //error if stepsize gets too low
    }

  if (errmax>ERRCON) hnext=SAFETY*h*pow(errmax, PGROW); //if stepsize is too small then grow stepsize
  else hnext=5.0*h;

  x+=(hdid=h); // new stepsize increment x + the old stepsize 

  for (i=0;i<n;i++) y[i]=ytemp[i];

}



int main() //our main program
{
  //Read the initial conditions in
  for (int i = 0; i < 7; i++) {
       input >> data[i];                
       }
  double r0=data[0];
  double a=data[3];
  double RHor=data[4];
  double range00=data[5];
  double a2=a*a;
  int m=data[6];
  int M=m*m;
  int n = N;

  //Setup vectors to do the integration of the 6+6 equations of motion plus change in y (yscal=y+dydx*h). 
  double y[6], dydx[6], yscal[6];
  double ylaststep[6]; 
  double xlaststep=0;

  // Define some parameters I will use later. Zmodel is the height of the reflection surface and y01 is the coordinates before intersection with the reflection surface and their momenta.
  double zmodel, y01[7];

  //Open file to save data in
  std::ofstream datafile;  
  datafile.open("position.dat");
  

  //Loop over the number of pixels in image plane
  for(int j=0; j<M; j++)  
    { cout <<"j=  "<< j << endl;

      // Set up intial try for step-size
      double htry=0.5, eps=1e-11, hdid=0.0, hnext=0.0;
      int i=0, q=0;
      const double TINY=1.0e-3;
      double x=0.0;
      

      //Restrict the emitted angles and therefore how big an area around the BH, which needs to be propagated 
      double range =range00/double(m-1);  
      
      // Definde the Image plane (x,y) or (p,t)
      // The x and y we use in initial conditions (see thesis for initial conditions for backwards ray-tracing) 
      double p, t; 

      t = -(j/m-m/2+0.5)*range;
      if(j/m != m) p = (j%m-m/2+0.5)*range; 
      else p = (0.5-m/2)*range; 


      //run the initial function to calculate the four velocity/momentum 
      initial(y, dydx, p, t); //I get dy and dydx out and I put coordinate in image plane in (p,t)

 
      double r2=y[0]*y[0];
      // define reflection surface, which is z=a1[0] rho^2+ a1[1]*rho (second order polynomial)
      double a1[2];
      a1[0]=0.02991911;
      a1[1]=0.03041085; 
      double x1=sqrt(r2 + a2)*sin(y[1])*cos(y[2]) ; //Cartesian x-coordinate
      double y1=sqrt(r2 + a2)*sin(y[1])*sin(y[2]) ; //Cartesian y-coordinate
      double xy=sqrt(x1*x1 + y1*y1) ;               //  cylindrical coordinate rho
      double z1=y[0]*cos(y[1]) ;  // cylindrical height z


      //Setting up stopping condition: if the photon is within a distance R_max of the BH, then calculate where the reflection surface is. If the photon is further away only let the reflection surface be negative. Stopping condition is that the photon for a given cylindrical radius xy is always above the reflection surface. Which is set to start/end at Rmax=400 R_g. 
      double Rmax=3000; 
      if(y[0]< Rmax) zmodel=a1[0]*pow(xy,2) + a1[1]*pow(xy,1) ; //Height of reflection surface when photons are within the Funnel.
      else zmodel=-5.0; //Height of reflection surface for photons not close enough to th BH, so they can always be propagated


// Kepp integrating as long as the photons is moving towards the BH, and it is propagating through the optically thin funnel and not the optically thick winds.
      while (y[0]<r0+100 &&z1>zmodel)

      // For a thin disk keep integrating untill photon hits the equitorial plane
      //while (y[0]<r0+100.0&&y[1]<PI/2)     // go till the particle falls into the BH or fly faraway
	{      //cout <<"zmodel=  "<< zmodel <<" z=  "<< z1 << endl;


	  //prepare to calculate next step and save the current position and step as ylaststep and xlaststep
	  for (i=0;i<n; i++) ylaststep[i]=y[i]; 
	  xlaststep=x; 

          //Calculate the four momenta/velocity in this place.
	  geodesic(y, dydx);  // calculate dy and dydlambda and call them y and dydx

          // Save the current position and four momenta, so I can do linear interpolation later
          y01[0]=y[0];
          y01[1]=y[1];
          y01[2]=y[2];
          y01[3]=y[3];
          y01[4]=L;
          y01[5]=y[4];
          y01[6]=y[5];


	  for (i=0; i<n; i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*htry)+TINY; 

          //the next y-value with the Runge-Kutta integrator
	  rkqs(y, dydx, x, htry, eps, yscal, hdid, hnext); // calculate the new position and velocity and next stepsize

	  //q++; 

	 //If the photon is close the BH horizon abort, since numerical under-flow can happen, and we won't see any emission
	  if(y[0]<RHor+0.001) 
	    {
	      y[0]=1000001;  
	      break;
	    }


	  htry=hnext;  //use the next stepsize to when calculating the new integration (adaptive step-size)
	  

	  //Define parameters as we did before while loop
          double r2=y[0]*y[0];
          double a1[2];
          a1[0]=0.02991911;  //tau=10:0.0125194
          a1[1]=0.03041085; //tau=10:0.13042059
          x1=sqrt(r2 + a2)*sin(y[1])*cos(y[2]) ;
          y1=sqrt(r2 + a2)*sin(y[1])*sin(y[2]) ;
          xy=sqrt(x1*x1 + y1*y1) ;
          z1=y[0]*cos(y[1]) ;

	  //Update zmodel to test if we are still inside the funnel or if we have hit the reflection surface, and therefore needs to stop 
          if(y[0]< Rmax+100.0) zmodel=a1[0]*pow(xy,2) + a1[1]*pow(xy,1) ; 
          else zmodel=-5.0;

            


	}

      // Check if the integration stopped because the photon hit the reflection surface through the funnel. 
      if (y[0]<Rmax-50&&y[0]>RHor+0.01 && y01[0]<Rmax-50&&y01[0]>RHor+0.01) 
	{
	  //paramters after hit
          double r_after=y[0]; 
          double theta_after=y[1];
          double phi_after=y[2];
	  double t_after=y[3];
          double pr_after=y[4];
          double ptheta_after=y[5];
          double pphi_after=L; //Constant of motion L
          double pt_after=energy1; //Energy (pt) constant of motion and normalized to E=-1


  	  //parmaters before hit
          double r_before=y01[0]; 
          double theta_before=y01[1];
          double phi_before=y01[2];
	  double t_before=y01[3];
          double pphi_before=y01[4];  
          double pr_before=y01[5];
          double ptheta_before=y01[6];





          //Save data before and after it hits the reflection surface, so I can do linear interpolation later.
	  datafile<<r_after<<" "<<theta_after<<" "<<phi_after<<" "<<pr_after<<" "<<ptheta_after<<" "<<pphi_after<<" "<<r_before<<" "<<theta_before<<" "<<phi_before<<" "<<pphi_before<<" "<<pr_before<<" "<<ptheta_before<<" "<<t_after<<" "<<t_before<<endl;
	}

      }
  datafile.close();
  return 0;
}
