#ifndef VEHICLE_H
#define VEHICLE_H
#include <iostream>
#include <random>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <iterator>

using namespace std;

class Vehicle {
public:
  
  /**
  * Constructor
  */
  Vehicle();

  /**
  * Copy Constructor
  */
  Vehicle(const Vehicle& orig);

  /**
  * Destructor
  */
  virtual ~Vehicle();

  /*
   * We will be using Frenet Coordinates for the car's position velocity and acceleration
   */
  void set_s(double offset, double vel, double acc);
  void set_d(double offset, double vel, double acc);
  vector<double> get_s() const;
  vector<double> get_d() const;
  vector<double> state_at(double t) const;
  double get_forecasted_param(int lag, int idx);
  void set_forecasted_param(int idx, double s_v1, double s_a, double d_v1, double d_a);
 private:
  double _s_s;
  double _s_v;
  double _s_a;
  double _d_d;
  double _d_v;
  double _d_a;
  vector<vector<double>> _forecasted_states;
 
};

#endif
