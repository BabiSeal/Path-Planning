#include <iostream>
#include <iostream>
#include <math.h>
#include <map>
#include <string>
#include <iterator>

#include "Vehicle.h"

Vehicle::Vehicle()
{
  _s_s = 0.0;
  _s_v = 0.0;
  _s_a = 0.0;
  _d_d = 0.0;
  _d_v = 0.0;
  _d_a = 0.0;
  for (int i = 0; i < 10; i++)
    {
      _forecasted_states.push_back({0.0, 0.0, 0.0, 0.0});
    }
}
Vehicle::Vehicle(const Vehicle& orig) {
}

Vehicle::~Vehicle() {
}

void Vehicle::set_s(double offset, double vel, double acc)
{
  _s_s = offset;
  _s_v = vel;
  _s_a = acc;
}

void Vehicle::set_d(double offset, double vel, double acc)
{
  _d_d = offset;
  _d_v = vel;
  _d_a = acc;

}

vector<double> Vehicle::get_s() const {
  return {_s_s, _s_v, _s_a};
}

vector<double> Vehicle::get_d() const {
  return {_d_d, _d_v, _d_a};
}

// returns frenet coordinates of predicted position at time t
// this is for the neighboring cars and not the ego. vehicle
vector<double> Vehicle::state_at(double t) const {
  double new_s = _s_s + t * _s_v;
  return {new_s, _d_d};
}
double Vehicle::get_forecasted_param(int lag, int idx) {
  return _forecasted_states[lag][idx];
}
void Vehicle::set_forecasted_param(int idx, double s_v1, double s_a, double d_v1, double d_a) {
  _forecasted_states[idx][0] = s_v1;
  _forecasted_states[idx][1] = s_a;
  _forecasted_states[idx][2] = d_v1;
  _forecasted_states[idx][3] = d_v1;
}
