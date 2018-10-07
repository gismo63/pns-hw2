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



class shoot {
  private:
    double a, b, t_f, c_max, h1, h2;
    int n1, n2;
    double ans, z;
  public:
    double t, c;
    Vec variables = Vec(2);
    //Vec zeros;
    shoot(double,double,double,double,double,double,int,int);
    Vec f(const Vec& x, double t);
    void steps();
    double guess();
    void graph();
};

shoot::shoot(double x_i_,double a_,double b_,double t_i_,double t_f_,double c_max_,int n1_,int n2_) {
  a=a_;
  b=b_;
  t=t_i_;
  t_f=t_f_;
  c_max=c_max_;
  n1=n1_;
  n2=n2_;
  c = -c_max_;
  h1 = (t_f_-t_i_)/n1_;
  h2 = (2*c_max_)/n2_;
  variables[0] = x_i_;
}

Vec shoot::f(const Vec& x,double t) {
  Vec y(2);
  y[0] = variables[1];
  y[1] = -(a*variables[1]+pow(variables[0],3)*exp(-b*t*variables[0]*variables[0]));
  return y;
}

void shoot::steps() {
  Vec k1(2), k2(2), k3(2), k4(2);
  k1 = h1*f(variables,t);
  k2 = h1*f((variables+k1/2),(t+h1/2));
  k3 = h1*f((variables+k2/2),(t+h1/2));
  k4 = h1*f((variables+k3),(t+h1));
  variables = variables + (k1 + 2*k2 + 2*k3 + k4)/6;
}

double shoot::guess() {
  for (int i=0; i<n1; i++) { //run the runge kutta algoithm n times
    steps();
    t += h1;
  }

  return variables[0]+1;
}

void shoot::graph() {
  for (int i=0; i<n2; i++) { //run the runge kutta algoithm n times
    variables[1] = c;
    z = guess();
    cout << c << ' ' << z << '\n';
    c += h2;
  }
}


int main(int argc, char const *argv[]) {
  Vec variables(2);
  variables[0]=0;//initial conditions
  double a=0.09, b=0.1, t_f=20, t_i = 0,c_max = 10;
  int n1 = 2000, n2 = 1000;//number of iterations

  shoot test(variables[0],a,b,t_i,t_f,c_max, n1, n2);
  test.graph();

  return 0;
}
