#include <iostream>
#include <vector>
#include <cmath>


typedef std::vector<double> Vec;
using std::cout;

#include "Planet.h"

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



std::vector<Planet> initial_cond(const std::vector<Vec>& x, const std::vector<Vec>& v, const Vec& m, int n) {
  std::vector<Planet> p;
  for(int i = 0; i<n; i++){
    p.push_back(Planet(m[i], x[i], v[i]));
  }
  return p;
}


Vec f(const std::vector<Planet>& p, int i, int d){
  Vec a(d);
  for(int k = 0; k<d; k++){
    a[k] = 0;
  }
  int n = p.size();
  for(int j = 0; j<n; j++){
    double dist_squared = 0;
    if(j != i) {
      for(int k = 0; k<d; k++){
        dist_squared += (p[i].x[k]-p[j].x[k])*(p[i].x[k]-p[j].x[k]);
      }
      for(int k = 0; k<d; k++){
        a[k] += p[j].mass()*(p[j].x[k]-p[i].x[k])/pow(dist_squared,1.5);
      }
    }
  }
  return a;
}

std::vector<Planet> v_step(double h, std::vector<Planet> p,int d){
  int n = p.size();
  for(int i = 0; i<n; i++){
    p[i].v = p[i].v + h*f(p, i, d);
  }
  return p;
}

std::vector<Planet> x_halfstep(double h, std::vector<Planet> p, int d){
  int n = p.size();
  for(int i = 0; i<n; i++){
    for(int j = 0; j<d; j++){
      p[i].x[j] += (h/2)*p[i].v[j];
    }
  }
  return p;
}


int main(int argc, char const *argv[]) {
  int n = 4; //number of planets
  int d = 2; //number of dimesions
  std::vector<Planet> p;
  std::vector<Vec> x;
  std::vector<Vec> v;
  x.resize(n);;
  v.resize(n);;
  for(int i = 0 ; i < n ; ++i){
    x[i].resize(d);
    v[i].resize(d);
  }
  Vec m(n);
  x[0][0] = -0.5, x[0][1] = 0.1, x[1][0] = -0.6, x[1][1] = -0.2, x[2][0] = 0.5, x[2][1] = 0.1, x[3][0] = 0.5, x[3][1] = 0.4;
  v[0][0] = -0.84, v[0][1] = 0.65, v[1][0] = 1.86, v[1][1] = 0.7, v[2][0] = -0.44, v[2][1] = -1.5, v[3][0] = 1.15, v[3][1] = -1.6;
  m[0] = 2.2, m[1] = 0.8, m[2] = 0.9, m[3] = 0.4;

  p = initial_cond(x,v,m,n);

  int steps = 40000;
  double t = 0, t_f = 4;

  double h = (t_f-t)/steps;

  for(int i = 0; i<steps; i++) {
    p = x_halfstep(h,p,d);
    p = v_step(h,p,d);
    p = x_halfstep(h,p,d);
    for(int i = 0; i<n; i++){
      for(int k = 0; k<d; k++){
        cout << p[i].x[k] << ' ';
      }
    }
    cout << '\n';
  }
  /**
  for(int i = 0; i<n; i++){
    for(int k = 0; k<d; k++){
      cout << p[i].x[k] << ' ';
    }
    cout << '\n';
  }
  **/
  return 0;
}
