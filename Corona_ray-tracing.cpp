#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>

using namespace std;
using std::ifstream;
#define PI 3.14159265358979323846
#define PI 3.14159265358979323846

#define N 6



#define a 0.8
//double a = 0.0;F
double a2 = a*a;

#define t0 0

/*inclined e=0.7 a=0.5 case       
double r0 = 7.06972;                 
double theta0 = PI/180.0*0.1;
double phi0 = 0.0;
double tdot0 = 1.0;
double rdot0 = 0.0;
double thetadot0 = 0.05999267;
double phidot0 = 0.0;


//inclined e=0.9 a=0.5 case

double r0 = 7.06972;
double theta0 = PI/180.0*0.1;
double phi0 = 0.0;
double tdot0 = 1.0;
double rdot0 = 0.0;
double thetadot0 = 0.0;
double phidot0 = 0.0;

*/






#define H 0.0  //0 for fotons and -1 for particles
double RHor=1+sqrt(1-a2)+0.00001; 
double RHor2=RHor*RHor ;
     //horizon event radius
#define RMstable 1.24
#define RDisc 200

double r0 = 10.06972;
double theta0 = PI/2.0;
double phi0 = PI;
double tdot0 = 1.0;
double rdot0 = 0.0;
double thetadot0 = 0.0;
double phidot0 = 0.0;


ifstream inFile;
#include <iostream>
#include <fstream>



double E;
double L;
double kappa;

void initial(double *y0, double *ydot0)
{

  double data[7];

  std::ifstream input("IP.txt"); //read data

  for (int i = 0; i < 7; i++) {
      input >> data[i];                 //save data 
      std::cout<< data[i]<<std::endl;
      }

  r0=data[0];
  theta0=data[1];
  phi0=data[2];
  tdot0=data[3];
  rdot0=data[4];
  thetadot0=data[5];
  phidot0=data[6];
  y0[0]=r0; //r0                       //we define the parameters radial
  y0[1]=theta0; //theta 0
  y0[2]=phi0; //phi 0
  y0[3]=t0;
    
  double r2=r0*r0;                //we have r²
  double r32=r0*sqrt(r0);         //r^3/2 

//trigonomic identities
  double costheta2=cos(theta0)*cos(theta0);
  double sintheta2=sin(theta0)*sin(theta0);
  double sum=r2+a2*costheta2;
  double delta=r2-2.0*r0+a2;


  double uaua=-(1-2.0*r0/sum)*tdot0*tdot0-4.0*a*r0*sintheta2*tdot0*phidot0/sum+sum*rdot0*rdot0/delta+sum*thetadot0*thetadot0+(r2+a2+2.0*r0*a2*sintheta2/sum)*sintheta2*phidot0*phidot0; //L  

  cout<<"uaua is "<< uaua << endl;
  cout<<"r0 "<< data[0] << "theta 0 " << data[1]<< "theta phi " << data[2] << "tdot 0 " << data[3] << "rdot 0" << data[4] << "theta_dot 0 " << data[5] << "phi_dot 0 " << data[6] << endl;


  double norm=sqrt(fabs(uaua));   //ds
  rdot0=rdot0/norm;          //dr/ds
  thetadot0=thetadot0/norm;  // dth / ds
  phidot0=phidot0/norm;      //dphi / ds
  tdot0=tdot0/norm;          //dt / ds
                                                          
  y0[4]= rdot0*sum/delta;    //pr       impuls in r
  y0[5]= thetadot0*sum;      //p_theta  impuls in theta

//set up dotted in new array
  ydot0[0] = rdot0; // rdot 0         
  ydot0[1] = thetadot0; // theta_dot
  ydot0[2] = phidot0; // phi_dot
  ydot0[3] = tdot0; //tdot 0
  cout << "tdot is " << tdot0 << " , and phidot is " <<phidot0<< endl <<" , and ds is " <<norm<< endl; 

  E = 2.0*a*r0*sintheta2*phidot0/sum-(-1.0+2.0*r0/sum)*tdot0;  //energy constant
  L=(sum*delta*phidot0-2.0*a*r0*E)*sintheta2/(sum-2.0*r0);     //angular momentum constant


//Print E and L
  cout <<r2*y0[0]-3.0*r2+2.0*a*r32<<endl;
  cout << "E is " <<E <<"; L is "<< L <<endl;
    

  kappa=y0[5]*y0[5]+a2*sintheta2*(E*E+H)+L*L/sintheta2;    // make kappa used to pr and p-theta
 
} 



void geodesic(double *y, double *dydlamda)  //make a geodesic function
{
  double r, theta, phi, t, pr, ptheta;
  
  r=y[0];
  theta= y[1];
  phi = y[2];
  t = y[3];
  pr = y[4];
  ptheta = y[5];

  //if(r>1000.0) r=0.001;
  double r2=r*r;
  double twor=2.0*r;
  double cos2=cos(theta)*cos(theta);
  double sin2=sin(theta)*sin(theta);
  double sum=r2+a2*cos2;
  double delta=r2-twor+a2;
  double sd=sum*delta;

  dydlamda[0] = pr*delta/sum;     //rdot
  dydlamda[1] = ptheta/sum;       //theta dot
  dydlamda[2] = (twor*a*E+(sum-twor)*L/sin2)/sd;   //phi dot
  dydlamda[3] = E+(twor*(r2+a2)*E-twor*a*L)/sd;    // t dot
  dydlamda[4] = ((r-1.0)*((r2+a2)*H-kappa)+r*delta*H+twor*(r2+a2)*E*E-2.0*a*E*L)/sd-2.0*pr*pr*(r-1.0)/sum;  // pr dot dot
  dydlamda[5] = sin(theta)*cos(theta)*(L*L/(sin2*sin2)-a2*(E*E+H))/sum;  //p-theta dot dot
  // if (y[1]>PI/2-0.001 &&  y[1]<PI/2+0.001)cout<<"theta is "<<y[1]<<"; ptheta is " << y[5]<<"; pthetadot is " <<dydlamda[5]<<endl; 
}

void rkck (double *y, double *dydx, double h, double *yout, double *yerr) //make integration RK method adaptive stepsize integral ck=Cash-Karp
{
  static const double b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0, b41 = 0.3, b42 = -0.9,
    b43 = 1.2, b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0,
    b54 = 35.0/27.0, b61 = 1631.0/55296.0, b62 = 175.0/512.0,
    b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
    c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0,
    dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25; //matrix of coefficient 

  
  int i;
  int n=N;
  double ak2[6], ak3[6], ak4[6], ak5[6], ak6[6], ytemp[6]; // we give y and ak a 6 dims array as defult

  for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i];    //First step
  geodesic(ytemp, ak2);
  
  
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);   //second step
  geodesic(ytemp, ak3);
 
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);  //3
  geodesic(ytemp, ak4);

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);  // 4 
  geodesic(ytemp, ak5);  

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);  // 5 last step
  geodesic(ytemp, ak6); 

  for (i=0;i<n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]); // output geodeisic
 
  for (i=0;i<n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]); // error between 5 and 6 order

  //if(yout[1]>2*PI) yout[1]=yout[1]-2*PI;
  //if(yout[1]<0) yout[1]=-yout[1];
  //if(yout[1]>PI) yout[1]=2*PI-yout[1];

}



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

void nrerror(const string error_text)
{
  cerr << "Numerical Recipes run-time error..." << endl;
  cerr << error_text << endl;
  cerr << "...now exiting to system..." << endl;
  exit(1);
}

#include <cmath>

void rkqs (double *y, double *dydx, double &x, double htry, double eps, double *yscal, double &hdid, double &hnext)  // if error small make h (step) larger
{
  const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
  
  int i;

  double errmax, h, htemp, xnew;

  int n=N;
  h=htry;
  double yerr[6], ytemp[6];  


  for (;;)   //break if errmax < 1
    {   
      rkck(y, dydx, h, ytemp, yerr);

      errmax=0.0;
      for(i=0;i<n;i++) errmax = MAX(errmax, fabs(yerr[i]/yscal[i]));  //errmax=eps1 (take highest error of vector)

      errmax/=eps;  
      if(errmax<=1.0) break; //if eps1=eps then we have obtained right h value. step succeded compute next step

      htemp=SAFETY*h*pow(errmax, PSHRNK);  //htemp=0.9*h*errmax^(-0.25). 
      h=(h>=0.0 ? MAX(htemp,0.1*h) : MIN(htemp, 0.1*h)); //if h>0 then h=htemp
    
      xnew=x+h;
      if(xnew==x) nrerror("stepsize underflow in rkqs");  //if xnew = x then it use h=0 and we get error
    }

  if (errmax>ERRCON) hnext=SAFETY*h*pow(errmax, PGROW); //if errmax > 2e-4 then use higher h
  else hnext=5.0*h;  //if errmax < 2e-4 then use hnext=5*h

  x+=(hdid=h);    //x= h_old + h_new

  for (i=0;i<n;i++) y[i]=ytemp[i];  //y_temporary=y

}



int main() //main
{
  int n = N;
  double y[6], dydx[6], yscal[6]; //make a 7th element in array

  double SIDE;
  double ylaststep[6], ynextstep[6];
  double xlaststep=0, xnextstep=0;

  std::ofstream datafile;
  datafile.open("geodesic_e7.dat");  //make/open datafile
 
  std::ofstream datafile2;
  datafile2.open("hitpoint.dat");
   

  double htry=0.5, eps=1e-11, hdid=0.0, hnext=0.0; //try, did(actual), next stepsize
  int i=0, q=0;
  const double TINY=1.0e-3;
  double x=0.0;
 
  //cout <<"**1";
  
  initial(y, dydx);
  // cout <<"**2";
  //datafile2<<E<<endl<<L<<endl;

  for (i=0; i<N; i++) cout << "y["<< i <<"] is " <<y[i]<<endl;     //print y and dy
  for (i=0; i<N; i++) cout << "ydot["<< i <<"] is " <<dydx[i]<<endl;
  int indicator=0;    
  double r1=0;
  double turned = 0;

  while ((y[0])<500&&y[0]>RHor&&fabs(y[1])<PI/2) //the radius is bound (not at r=0)
    {	  
      for (i=0;i<n; i++) ylaststep[i]=y[i];
      xlaststep=x; 
       
      geodesic(y, dydx);
      
      for (i=0; i<n; i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*htry)+TINY;
      


     
      if (y[1]>PI/2) SIDE= 1.0; //we work on two sides theta>2pi/2 (side 0) and theta<2pi/2(side -1)
      else if (y[1]<PI/2) SIDE = -1.0;
      else SIDE = 0.0;
	

      rkqs(y, dydx, x, htry, eps, yscal, hdid, hnext);
      q++;

      if(y[0]<RHor+0.00001) //if r < horizon radius (particle are absorped at event horizon)=case 3
	{
	  y[0]=100001;  
	  break;
	}
     
      else if (q>4000) //case 4 = stable over 1000 integration and then stop
	{
	  //cout<<"CASE 4: stable"<<endl;  //stable
	  break;
	}
        
      htry=hnext;
      double costheta2=cos(y[1])*cos(y[1]);
      double sintheta2=sin(y[1])*sin(y[1]);
      double sum=y[0]*y[0]+a2*costheta2;
      double delta=y[0]*y[0]-2.0*y[0]+a2;


      
      double X=sqrt(y[0]*y[0]+a*a)*sin(y[1])*cos(y[2]);  //back to cartesian easy to plot
      double Y=sqrt(y[0]*y[0]+a*a)*sin(y[1])*sin(y[2]);
      double Z=y[0]*cos(y[1]);
      double r = sqrt(X*X+Y*Y+Z*Z);
      double E1 = 2.0*a*y[0]*sintheta2*dydx[2]/sum-(-1.0+2.0*y[0]/sum)*dydx[3];  //energy constant
      double L1=(sum*delta*dydx[2]-2.0*a*y[0]*E1)*sintheta2/(sum-2.0*y[0]);     //angular momentum constant


      datafile2<<X<<" "<<Y<<" "<<Z<<" "<<y[3]<<" "<<y[4]<<" "<<y[5]<<" "<<E<<" "<<L<<" "<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl; 
      datafile<<E1<<" "<<L1<<" "<<y[0]<<" "<<y[1]<<" "<<y[2]<<endl;

	//datafile << y[0] << " " << y[1] << " " << y[2]<<" " << dydx[0]<< " " <<dydx[1] << " " <<dydx[2]<< " " <<dydx[3]<<endl;

	//datafile<<y[5]<<" " <<y[1]<<" " <<Y<<" "<<X<<" "<<Z<<endl;
    }
  
  datafile<<endl;
  datafile.close();
  datafile2.close();
  return 0;
}




