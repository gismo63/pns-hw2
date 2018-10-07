#include <iostream>
#include <vector>
#include <cmath>

typedef std::vector<double> Vec;
using std::cout;


//Here we overload the * operator so that a Vec type can be multiplied by a scalar on the right
Vec operator*(Vec x, double factor) {
  for (size_t i=0;i<x.size();i++) { // Here we use size_t because x.size() is an unsigned integer so we also make i and unsigned integer for comparison purposes
    x[i] = (x[i]*factor);
  }
  return x;
}

//same as above but by scalar on left
Vec operator*(double factor, Vec x) {
  for (size_t i=0;i<x.size();i++) {
    x[i] = (x[i]*factor);
  }
  return x;
}

//Here we overload the / operator so that a Vec type can be divided by a scalar
Vec operator/(Vec x, double denom) {
  for (size_t i=0;i<x.size();i++) {
    x[i] = (x[i]/denom);
  }
  return x;
}

//Here we overload the + operator so that two Vecs can be added
Vec operator+(Vec x, const Vec& y){
  for (size_t i=0;i<x.size();i++) {
    x[i] = (x[i] + y[i]);
  }
  return x;
}

//This function takes all x_i in vector form and returns f_i(x,t) in vector form where f_i = d(x_i)/dt
Vec f(const Vec& x, double t, double a, double b){
  Vec y(2);
  y[0] = x[1];
  y[1] = -(a*x[1]+pow(x[0],3)*exp(-b*t*x[0]*x[0]));
  return y;
}

//This function runs one iteration of the runge-kutta algorithm for any number of first order differential equations, takes x[t] and returns x[t+h]
Vec steps(Vec x, double t, double h, double a, double b) {
  Vec k1(2), k2(2), k3(2), k4(2);
  k1 = h*f(x,t,a,b);
  k2 = h*f((x+k1/2),(t+h/2),a,b);
  k3 = h*f((x+k2/2),(t+h/2),a,b);
  k4 = h*f((x+k3),(t+h),a,b);
  x = x + (k1 + 2*k2 +2*k3 + k4)/6;
  return x;
}

double zeros(Vec variables, double t, double h, int n, double a, double b) {
  for (int i=0; i<n; i++) { //run the runge kutta algoithm n times
    variables = steps(variables,t,h,a,b);
    t += h;
  }
  return variables[0]+1;

}

double solution(double accuracy,double low_lim,double up_lim, Vec variables, double t, double h, int n, double a, double b) {
  variables[1] = low_lim;
  double c_neg, c_pos;
  if (std::signbit(zeros(variables, t, h, n, a, b))) {
    c_neg = low_lim;
    c_pos = up_lim;
  }
  else {
    c_pos = low_lim;
    c_neg = up_lim;
  }
  while (fabs(c_pos-c_neg) > accuracy) {
    variables[1] =  (c_pos+c_neg)/2;
    double avg_out = zeros(variables, t, h, n, a, b);
    if (std::signbit(avg_out)) {
      c_neg = variables[1];
    }
    else {
      c_pos = variables[1];
    }
  }
  cout << abs(c_pos-c_neg) <<"\n";
  return (c_pos+c_neg)/2;
}


int main(int argc, char const *argv[]) {
  Vec variables(2);
  variables[0]=0;//initial conditions
  double a=0.09, b=0.1, t_f=20, t = 0,c_max = 5, c = -5, z;
  int n = 20000, n2 = 1000;//number of iterations
  double h = (t_f-t)/n, h2 = (c_max-c)/n2;//step size
  /**
  for (int i=0; i<n2; i++) { //run the runge kutta algoithm n times
    variables[1] = c;
    z = zeros(variables, t, h, n, a, b);
    cout << c << ' ' << z << '\n';
    c += h2;
  }
  **/
  double acc = 0.0000000001;
  double sol = solution(acc,1.94,1.95,variables, t, h, n, a, b);
  variables[1] = sol;
  cout.precision(10);
  cout << "v(0) of " << sol << " gives x(20) of " << zeros(variables, t, h, n, a, b) << "\n";

  return 0;
}
