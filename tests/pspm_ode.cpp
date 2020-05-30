#include "pspm_ode_solver2.h"
#include <fstream>

void simple_pendulum(double x, const vector<double>& y, vector<double>& dydx){ 
	// use time units omega*t->t
  // d^2 u/dt^2 = - omega^2 u    is written as
  // y[0,1] = [u(t), du/dt]
  // d y/dt = [du/dt, d^2u/dt] = [du/dt, -u] = [y[1],-y[0]]
  // Exact solution is y[0,1] = [xmax*sin(t), xmax*cos(t)]
  dydx[0] = y[1];
  dydx[1] = -y[0];
}


double Energy(const vector<double>& y){
	// Energy of the pendulum is E = 1/2*xmax ( u^2 +  (du/dt)^2 )
  return 0.5*(y[0]*y[0]+y[1]*y[1]);
}


int main(){
  double E = 1;                 // our choice for energy in dimensionless units
  double t_start=0, t_stop=200; // elapsed time in dimensionless units
  double dh=0.1;                // time step

  double xmax = sqrt(2*E);      // Energy uniquely determines maximum amplitude
  vector<double> y(2);
  y[0]=0;                       // Pendulum is initially at zero but has maximum momentum
  y[1]=xmax;
  // The exact solution of this problem is:
  //    y[0,1] = [xmax*sin(t), xmax*cos(t)]


  RKCK45<vector<double> > rk(t_start, t_stop, 1e-8); // RK class for adaptive step
  
  cout.precision(15);
  int M = static_cast<int>((t_stop-t_start)/dh+0.5); // Number of timesteps
  double ti = t_start;
  
  /*
  for (int i=0; i<M; i++, ti+=dh) {
    cout<<setw(25)<<ti<<" "<<setw(25)<<y[0]<<" "<<setw(25)<<y[1]<<" "<<setw(25)<<y[0]-xmax*sin(ti)<<" "<<setw(25)<<Energy(y)<<endl;
    //Euler(ti, dh, y, simple_pendulum);
    RK4(ti, dh, y, simple_pendulum);
  }
  */
  

//	ofstream fout("pendulum.txt");
  for (; ; ) {
//	fout<<ti<<" "<<y[0]<<" "<<y[1]<<" "<<y[0]-xmax*sin(ti)<<" "<<Energy(y)<<endl;
    if (!rk.Step(ti, y,  simple_pendulum)) break;
    
	if (fabs(y[0] - xmax*sin(ti)) > 1e-6) return 1;
    if (fabs(y[1] - xmax*cos(ti)) > 1e-6) return 1;

  }

  if (rk.size() != 2) return 1;
//  fout<<ti<<" "<<y[0]<<" "<<y[1]<<" "<<y[0]-xmax*sin(ti)<<" "<<Energy(y)<<endl;
  
  return 0;
}

