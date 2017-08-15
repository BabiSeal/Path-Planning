#ifndef PTG_H
#define PTG_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <math.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/LU"
#include <random>
#include "Vehicle.h"

using namespace std;


class PTG {
 public:

  PTG();

  ~PTG();

  double poly_eval(vector<double> coeffs, double x);
  double poly_differentiate_and_eval(vector<double> coeffs, double x);
  double poly_differentiate_twice_and_eval(vector<double> coeffs, double x);
  double poly_differentiate_thrice_and_eval(vector<double> coeffs, double x);
  
  int d_to_lane(double d);
  
  void start_gracefully_from_complete_stop(double start_velocity);

  void realize_keep_lane_state_trajectories(int curr_lane,
                                            bool follow_leader,
                                            vector<double> const &start_s,
                                            vector<double> const &lead_s,
                                            vector<Vehicle> const &vehicles,
                                            vector<vector<double>>& all_goals);

  void realize_hard_stop_state_trajectories(int curr_lane,
                                            bool follow_leader,
                                            vector<double> const &start_s,
                                            vector<double> const &lead_s,
                                            vector<Vehicle> const &vehicles,
                                            vector<vector<double>>& all_goals);
  
  void realize_lane_change_left_trajectories(int curr_lane,
                                             bool follow_leader,
                                             vector<double> const &start_s,
                                             vector<double> const &lead_s,
                                             vector<Vehicle> const &vehicles,
                                             vector<vector<double>>& all_goals);

  void realize_lane_change_right_trajectories(int curr_lane,
                                              bool follow_leader,
                                              vector<double> const &start_s,
                                              vector<double> const &lead_s,
                                              vector<Vehicle> const &vehicles,
                                              vector<vector<double>>& all_goals);

  void generate_next_state_trajectories(vector<double> const &start,
                                        vector<Vehicle> const &vehicles,
                                        vector<vector<double>>& all_goals);
  
  vector<vector<double>> generate_sd_path(vector<double> const &start,  double max_speed,
                                          double horizon, vector<Vehicle> const &vehicles);

  void perturb_goal(vector<double> goal, vector<vector<double>> &all_goals); 

  vector<double> jmt(vector<double> const &start, vector<double> const &goal, int t);
  
  int find_closest_behind_vehicle_for_lane(int lane,
                                           vector<double> const &start,
                                           vector<Vehicle> const &vehicles);

  bool check_behind_vehicle_approach(vector<double> const &start,
                                     Vehicle const &from_back);


  int find_closest_ahead_vehicle_for_lane(int lane,
                                    vector<double> const &start,
                                    vector<Vehicle> const &vehicles);
  
  vector<int> find_closest_ahead_vehicle_in_each_lane(vector<double> const &start,
                                                vector<Vehicle> const &vehicles);
  
  double exceeds_speed_limit_cost(pair<vector<double>, vector<double>> const &traj,
                                  vector<double> const &goal,
                                  vector<Vehicle> const &vehicles);
  
  double exceeds_accel_cost(pair<vector<double>, vector<double>> const &traj,
                            vector<double> const &goal,
                            vector<Vehicle> const &vehicles);
  
  double exceeds_jerk_cost(pair<vector<double>, vector<double>> const &traj,
                           vector<double> const &goal,
                           vector<Vehicle> const &vehicles);
  
  double collision_cost(pair<vector<double>, vector<double>> const &traj,
                        vector<double> const &goal,
                        vector<Vehicle> const &vehicles);

  double lane_selection_cost(pair<vector<double>, vector<double>> const &traj,
                             vector<double> const &goal,
                             vector<Vehicle> const &vehicles);

  double buffer_cost(pair<vector<double>, vector<double>> const &traj,
                     vector<double> const &goal,
                     vector<Vehicle> const &vehicles);

  double efficiency_cost(pair<vector<double>, vector<double>> const &traj,
                         vector<double> const &goal,
                         vector<Vehicle> const &vehicles);

  double total_accel_s_cost(pair<vector<double>, vector<double>> const &traj,
                            vector<double> const &goal,
                            vector<Vehicle> const &vehicles);
  
  double total_accel_d_cost(pair<vector<double>, vector<double>> const &traj,
                            vector<double> const &goal,
                            vector<Vehicle> const &vehicles);

  double total_jerk_cost(pair<vector<double>, vector<double>> const &traj,
                         vector<double> const &goal,
                         vector<Vehicle> const &vehicles);

  double logistic(double x);
   
  double calculate_cost(pair<vector<double>, vector<double>> const &traj,
                        vector<double> const &goal,
                        vector<Vehicle> const &vehicles,
                        vector<vector<double>> &all_costs);

  string get_current_state_category();

  const string kKeepLane = "Keep Lane";
  const string kKeepLaneHardStopAhead = "Slow down immediately";
  const string kLaneChange = "Changing Lanes";
  const string kNudgeToMiddleLane = "Nudge to Middle Lane";

 private:
  string _current_state_category;

  enum State { KEEP_LANE = 0,
               KEEP_LANE_HARD_STOP_AHEAD,
               LANE_CHANGE_LEFT,
               LANE_CHANGE_RIGHT };

  vector <string> _state_names = {"Keep Lane",
                                  "Keep Lane, Hard Stop Ahead",
                                  "Lane Change Left",
                                  "Lane Change Right"};
    
  int _num_steps_ahead;
  double _max_vel_in_ts = 0.0;
  double _delta_s_with_max_vel  = 0.0;

  const double kMaxLookAhead = 200;

  const double kMaxVelPerUnitTimeStep = 0.00894 * 47.5; 
  const double kMaxAccPerUnitTimeStep = (10.0 / 50.0)*0.90; // 10 m/s
  const double kMaxJerkPerUnitTimeStep = 10.0 / 50.0; // 10 m/s

  const int kNumSamplesForPerturb = 10;

  const double kCarWidth = 2.0;
  const double kCarLength = 5.0;
  const double kCarCollisionWidth = 0.5 * kCarWidth;
  const double kCarCollisionLength = 0.5 * kCarLength;
  const double kCollisionBufferWidth = kCarWidth;
  const double kCollisionBufferLength = 6 * kCarLength;

 
  
  std::default_random_engine _generator;

  std::map<std::string, double> _cost_weights = {
     {"lane_selection_cost",  30.0}, 
     {"buffer_cost",         190.0},
     {"eff_cost",            110.0}, 
     {"tot_acc_s_cost",       30.0},
     {"tot_acc_d_cost",       20.0},
     {"tot_jerk_cost",        10.0},
   }; 
};

#endif /* PTG_H */
