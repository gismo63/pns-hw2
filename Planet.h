class Planet{
  private:
    double m;
  public:
    Vec x;
    Vec v;

    Planet(double m_, Vec x_, Vec v_);


    double mass() const { return m; }
};

Planet::Planet(double m_, Vec x_, Vec v_) {
  m = m_;
  int d = x_.size();
  x.resize(d);
  v.resize(d);
  for(int i = 0; i<d; i++){
    x[i] = x_[i];
    v[i] = v_[i];
  }
}
