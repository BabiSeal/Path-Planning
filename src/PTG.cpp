#include "PTG.h"

PTG::PTG() {
}

PTG::~PTG() {
}

double PTG::poly_eval(vector<double> coeffs, double x)
{
  double total = 0;
  for (int i = 0; i < coeffs.size(); i++)
    total = total + coeffs[i]*pow(x, i);
  return total;
}

double PTG::poly_differentiate_and_eval(vector<double> coeffs, double x)
{
  double total = 0;
  vector<double> d_coeffs;
  for (int i = 0; i < coeffs.size(); i++) {
    if (i > 0) {
      double d = i * coeffs[i];
      d_coeffs.push_back(d);
    }
  }
  for (int i = 0; i < d_coeffs.size(); i++)
    total = total + d_coeffs[i]*pow(x, i);
  return total;
}

double PTG::poly_differentiate_twice_and_eval(vector<double> coeffs, double x)
{
  double total = 0;
  vector<double> dd_coeffs;
  for (int i = 0; i < coeffs.size(); i++) {
    if (i > 0) {
      double d = i * coeffs[i];
      if (i > 1) {
        dd_coeffs.push_back((i - 1) * d);
      }
    }
  }
  for (int i = 0; i < dd_coeffs.size(); i++)
    total = total + dd_coeffs[i]*pow(x, i);
  return total;
}

double PTG::poly_differentiate_thrice_and_eval(vector<double> coeffs, double x)
{
  double total = 0;
  vector<double> ddd_coeffs;
  for (int i = 0; i < coeffs.size(); i++) {
    if (i > 0) {
      double d = i * coeffs[i];
      if (i > 2) {
        ddd_coeffs.push_back((i - 2) * d);
      }
    }
  }
  for (int i = 0; i < ddd_coeffs.size(); i++)
    total = total + ddd_coeffs[i]*pow(x, i);
  return total;
}


// creates randomly generated variations of goal point
void PTG::perturb_goal(vector<double> goal, vector<vector<double>> &perturbed_goals_vec) {

  vector<double> perturbed_goal(6);
 
  double sigma_sq = 0.1;
  std::normal_distribution<double> distrib(0.0, sigma_sq);

  for (int i = 0; i < kNumSamplesForPerturb; i++) {
    // Just do not want to exceed max parameters for s and v
    double factor = -1.0*distrib(_generator);
    perturbed_goal[0] = goal[0] + (_delta_s_with_max_vel * factor);
    perturbed_goal[1] = goal[1] + (_max_vel_in_ts * factor);
    perturbed_goal[2] = 0.0;

    // Perturbation in d. Do not mess with d_v or d_a
    factor = distrib(_generator);
    perturbed_goal[3] = goal[3] + factor;
    perturbed_goal[4] = 0.0;
    perturbed_goal[5] = 0.0;
    perturbed_goals_vec.push_back(perturbed_goal);
  }
}

vector<double> PTG::jmt(vector<double> const &start, vector<double> const &goal, int t) {
  
  double T = double(t);
  double t_2 = pow(T, 2);
  double t_3 = pow(T, 3);
  double t_4 = pow(T, 4);
  double t_5 = pow(T, 5);
  Eigen::Matrix3d A;
  A << t_3,   t_4,    t_5,
       3*t_2, 4*t_3,  5*t_4,
       6*t,   12*t_2, 20*t_3;
  
  double b_0 = start[0] + start[1] * t + 0.5 * start[2] * t_2;
  double b_1 = start[1] + start[2] * t;
  double b_2 = start[2];
  Eigen::MatrixXd b(3,1);
  b << goal[0] - b_0, goal[1] - b_1, goal[2] - b_2;
  
  Eigen::MatrixXd c = A.inverse() * b;
  vector<double> coeff = {start[0], start[1], 0.5*start[2], c.data()[0], c.data()[1], c.data()[2]};
 
  return coeff;
}



int PTG::d_to_lane(double d)
{
  return ((int)(round(round(d - 2.0) / 4.0)));
}

int PTG::find_closest_ahead_vehicle_for_lane(int lane,
                                  vector<double> const &start,
                                  vector<Vehicle> const &vehicles)
{
  int closest = -1;
  double min_dist = 9999;
  for (int i = 0; i < vehicles.size(); i++) {
    vector <double> other_car_d = vehicles[i].get_d();
    int other_car_lane = d_to_lane(other_car_d[0]);
    if (other_car_lane == lane)
      {
        vector <double> other_car_s = vehicles[i].get_s();
        double dist = other_car_s[0] - start[0];
        if (dist > 0.0 && dist < min_dist) {
          min_dist = dist;
          closest = i;
        }
      }
  }
  //cout << "    Lane " << lane << " Vehicle Idx " << closest << " at dist " << min_dist <<  endl;
  return closest;
}


int PTG::find_closest_behind_vehicle_for_lane(int lane,
                                  vector<double> const &start,
                                  vector<Vehicle> const &vehicles)
{
  int closest = -1;
  double min_dist = -9999;
  for (int i = 0; i < vehicles.size(); i++) {
    vector <double> other_car_d = vehicles[i].get_d();
    int other_car_lane = d_to_lane(other_car_d[0]);
    if (other_car_lane == lane)
      {
        vector <double> other_car_s = vehicles[i].get_s();
        double dist = other_car_s[0] - start[0];
        if (dist < 0.0 && dist > min_dist) {
          min_dist = dist;
          closest = i;
        }
      }
  }
  //cout << "    Behind Car: Lane " << lane << " Vehicle Idx " << closest << " at dist " << min_dist <<  endl;
  return closest;
}

bool PTG::check_behind_vehicle_approach(vector<double> const &start,
                                        Vehicle const &from_back)
{
  bool ok_room = false;
  vector <double> behind_car_s = from_back.get_s();
  double behind_dist = abs(start[0] - behind_car_s[0]);
  double behind_vel = behind_car_s[1];
  if (((behind_vel < start[1]) && (behind_dist > (0.5* kCollisionBufferLength))) ||
      ((behind_vel > start[1]) && (behind_dist > kCollisionBufferLength)))
    {
      ok_room = true;
    }

  // if (ok_room == false)
  // {
      //cout << "    Warning: Behind Car in adjacency lane dist: " << behind_dist << " vel " << behind_vel <<  "safe: " << ok_room << endl;
  //}
  return ok_room;
}

vector<int> PTG::find_closest_ahead_vehicle_in_each_lane(vector<double> const &start,
                                                   vector<Vehicle> const &vehicles)
{
  cout << "    Car s: " << start[0] << endl;
  vector<int> closest_vehicles(3);
  for (int i = 0; i < 3; i++) {
    closest_vehicles[i] = find_closest_ahead_vehicle_for_lane(i, start, vehicles);
    if (closest_vehicles[i] != -1)
      {
        vector <double> cr_s = vehicles[closest_vehicles[i]].get_s();
        cout << "    Lane " << i << " Vehicle Idx " << closest_vehicles[i]
             << " at dist " << (cr_s[0] - start[0]) << " w vel " << cr_s[1] <<  endl;
      }
  }
  
  
  return closest_vehicles;
}

double PTG::exceeds_speed_limit_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  for (int i = 0; i < _num_steps_ahead; i++) {
    double curr_vel_s = poly_differentiate_and_eval(traj.first, i);
    double curr_vel_d = poly_differentiate_and_eval(traj.second, i);
    if (curr_vel_s + curr_vel_d > kMaxVelPerUnitTimeStep) {
      //cout << " Exceeded velocity w " << (curr_vel_s + curr_vel_d) << " hard max is " << kMaxVelPerUnitTimeStep << endl;
     return 1.0;
    }
  }
  return 0.0;
}

double PTG::exceeds_accel_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  for (int i = 0; i < _num_steps_ahead; i++) {
    double curr_accel_s = poly_differentiate_twice_and_eval(traj.first, i);
    double curr_accel_d = poly_differentiate_twice_and_eval(traj.second, i);
    if ( curr_accel_s + curr_accel_d > kMaxAccPerUnitTimeStep) {
      //cout << " Exceeded accel w " << (curr_accel_s + curr_accel_d) << " hard max is " << kMaxAccPerUnitTimeStep << endl;
      return 1.0;
    }
  }
  return 0.0;
}

double PTG::exceeds_jerk_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  for (int i = 0; i < _num_steps_ahead; i++) {
    double curr_jerk_s = poly_differentiate_thrice_and_eval(traj.first, i);
    double curr_jerk_d = poly_differentiate_thrice_and_eval(traj.second, i);
    if (curr_jerk_s + curr_jerk_d > kMaxJerkPerUnitTimeStep) {
      //cout << " Exceeded jerk w " << (curr_jerk_s + curr_jerk_d) << " hard max is " << kMaxJerkPerUnitTimeStep << endl;
      return 1.0;
    }
  }
  return 0.0;
}

double PTG::collision_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  for (int t = 0; t < _num_steps_ahead; t++) {
    for (int i = 0; i < vehicles.size(); i++) {
      double self_s = poly_eval(traj.first, t);
      double self_d = poly_eval(traj.second, t);
      vector<double> other_car = vehicles[i].state_at(t); // {s,d}

      double dist_s = abs(other_car[0] - self_s);
      double dist_d = abs(other_car[1] - self_d);
      
      if ((dist_s <= 5*kCarCollisionLength) && (dist_d <= 3*kCarCollisionWidth)) {
        //cout << " Exceeded collision cost s-dist " << dist_s << " d-dist " << dist_d << endl;
        return 1.0;
      }
    }
  }
  return 0.0;
}

// adds cost for getting too close to another vehicle
double PTG::buffer_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  double cost = 0.0;
  for (int i = 0; i < vehicles.size(); i++) {
    for (int t = 0; t < _num_steps_ahead; t++) {
      double self_s = poly_eval(traj.first, t);
      double self_d = poly_eval(traj.second, t);
      vector<double> other_state = vehicles[i].state_at(t); // {s,d}
      
      double dist_s = other_state[0] - self_s;

      // If cars behind are at safe distance
      if (dist_s < 0 && (dist_s < -1.0*(2 * kCarLength)))
        continue;
  
      dist_s = abs(dist_s);
      double dist_d = abs(other_state[1] - self_d);

      // Only if we are in the danger zone
      if ((dist_s <= kCollisionBufferLength) && (dist_d <= kCollisionBufferWidth))
        cost += logistic(1 - (dist_s / kCollisionBufferLength)) / _num_steps_ahead;
    }
  }
  return cost;
}

// Efficiency - If we way overshoot our goal or undershoot our goal, it is inefficient 
double PTG::efficiency_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  double s_dist = goal[0] - poly_eval(traj.first, 0);
  double max_dist = _delta_s_with_max_vel;
  return abs(logistic((max_dist - s_dist) / max_dist)); 
}

double PTG::total_accel_s_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  double cost = 0.0;
  for (int t = 0; t < _num_steps_ahead; t++) {
    cost += abs(poly_differentiate_twice_and_eval(traj.first, t));
  }
  return logistic(cost);
}

double PTG::total_accel_d_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles) {
  double cost = 0.0;
  for (int t = 0; t < _num_steps_ahead; t++) {
    cost += abs(poly_differentiate_twice_and_eval(traj.second, t));
  }
  return logistic(cost);
}

double PTG::total_jerk_cost(pair<vector<double>, vector<double>> const &traj,
                            vector<double> const &goal,
                            vector<Vehicle> const &vehicles) {
  double cost = 0.0;
  for (int t = 0; t < _num_steps_ahead; t++) {
    cost += abs(poly_differentiate_thrice_and_eval(traj.first, t));
    cost += abs(poly_differentiate_thrice_and_eval(traj.second, t));
  }
  return logistic(cost);
}

// Lane Selection Cost: 
double PTG::lane_selection_cost(pair<vector<double>, vector<double>> const &traj,
                                vector<double> const &goal,
                                vector<Vehicle> const &vehicles) {
  
  double ego_s = poly_eval(traj.first, 0);
  double ego_s_vel = poly_differentiate_and_eval(traj.first, 0);

  double ego_d = poly_eval(traj.second, 0);
  double ego_d_tgt = poly_eval(traj.second, _num_steps_ahead);

  int tgt_lane_i = d_to_lane(ego_d_tgt);
  int cur_lane_i = d_to_lane(ego_d);

  
  // If we are already nudging a car to the middle lane, from an outer lane,
  // then we do not want to add this cost
  if (_current_state_category == kNudgeToMiddleLane) {
    cout << "SKIP Lane Selection Cost check when we are trying to nudge to mid lane " << endl;
    return 0;
  }

  vector<double> start = {ego_s, ego_s_vel, 0, ego_d, 0, 0};
  int closest_car_tgt_i = find_closest_ahead_vehicle_for_lane(tgt_lane_i, start,  vehicles);

  // if there is a vehicle in the target lane
  if (closest_car_tgt_i != -1) {
    vector<double> tgt_lane_cars_s = vehicles[closest_car_tgt_i].get_s();
    float dif_s = tgt_lane_cars_s[0] - ego_s;
    if (dif_s < 200) {

      int closest_car_i = find_closest_ahead_vehicle_for_lane(cur_lane_i, start, vehicles);
      
      // Again, remember we are checking only if there is a vehicle in the current lane
      // Also we check for distance < 100
      if ((closest_car_i != -1) && (dif_s < 100)) {
        vector<double> curr_lane_cars_s = vehicles[closest_car_i].get_s();

        // First check if traffic in planned lane slower than in current ?
 
        // cout << "Lane Selection: Closest Vehicle in current lane dist is " << (curr_lane_cars_s[0] - ego_s)  \
        // << " distance to target" << dif_s <<  " curr_lane :  " << cur_lane_i << " target lane " << tgt_lane_i \
        // << " close vel " << traffic_s[1] << " tgt vel " << tgt_lane_cars_s[1] << endl;
        if (tgt_lane_cars_s[1] < 0.99*curr_lane_cars_s[1])
        {
          // cout << " SLOW car in target lane" << endl;
          return 1000;
        }
        // If this a planned lane change, check if there is a vehicle coming from behind
        if (cur_lane_i != tgt_lane_i)
          {
            int behind_car_idx = find_closest_behind_vehicle_for_lane(tgt_lane_i, start, vehicles);
            if (behind_car_idx != -1)
              {
                bool car_in_target_lane_w_behind_room = check_behind_vehicle_approach(start,
                                                                                      vehicles[behind_car_idx]);
                if (car_in_target_lane_w_behind_room == false)
                  {
                    //cout << " RISKY FROM BEHIND - JUST DO NOT DO IT, NOT WORTH IT!!!!!" << endl;
                    return 1000;
                  }
              }
          }
      }
      return logistic(1 - (dif_s / 200));
    }
  }
  return 0.0;
}


// returns a value between 0 and 1 for x in the range [0, infinity]
// and -1 to 1 for x in the range [-infinity, infinity].
// approaches 1 at an input of around 5
double PTG::logistic(double x) {
    return (2.0 / (1 + exp(-x)) - 1.0);
}

double PTG::calculate_cost(pair<vector<double>, vector<double>> const &traj, vector<double> const &goal, vector<Vehicle> const &vehicles, vector<vector<double>> &all_costs) {
  
  double cost = 0.0;
  /*
   * First the Infeasible costs
   */
  double ex_sp_lim_cost = exceeds_speed_limit_cost(traj, goal, vehicles);
  double ex_acc_lim_cost = exceeds_accel_cost(traj, goal, vehicles);
  double ex_jerk_lim_cost = exceeds_jerk_cost(traj, goal, vehicles);
  double col_cost = collision_cost(traj, goal, vehicles);
  
  double infeasible_costs = ex_sp_lim_cost + ex_acc_lim_cost + ex_jerk_lim_cost + col_cost;
  if (infeasible_costs > 0.0) {
    //cout << "Infeasible Costs exceeded " << infeasible_costs << endl;
    //cout << "vals are " << ex_sp_lim_cost << " " <<  ex_acc_lim_cost << " " << ex_jerk_lim_cost << " " << col_cost << endl;
    all_costs.push_back({999999});
    return 999999;
  }  

  double ln_selection_cost = lane_selection_cost(traj, goal, vehicles) * _cost_weights["lane_selection_cost"];
  double buf_cost          = buffer_cost(traj, goal, vehicles) * _cost_weights["buffer_cost"];
  double eff_cost          = efficiency_cost(traj, goal, vehicles) * _cost_weights["eff_cost"];
  double tot_acc_s_cost    = total_accel_s_cost(traj, goal, vehicles) * _cost_weights["tot_acc_s_cost"];
  double tot_acc_d_cost    = total_accel_d_cost(traj, goal, vehicles) * _cost_weights["tot_acc_d_cost"];
  double tot_jerk_cost     = total_jerk_cost(traj, goal, vehicles) * _cost_weights["tot_jerk_cost"];
  
  
  vector<double> cost_vec = {ln_selection_cost, buf_cost, eff_cost, tot_acc_s_cost, tot_acc_d_cost, tot_jerk_cost};
  all_costs.push_back(cost_vec);
  
  cost = ln_selection_cost + buf_cost + eff_cost + tot_acc_s_cost + tot_acc_d_cost + tot_jerk_cost;
  return cost;
}

void PTG::realize_keep_lane_state_trajectories(int cur_lane_i,
                                               bool follow_leader,
                                               vector<double> const &start_s,
                                               vector<double> const &lead_s,
                                               vector<Vehicle> const &vehicles,
                                               vector<vector<double>>& all_goals)

{
  double goal_s_pos;
  double goal_s_vel;
  double goal_s_acc;
  double goal_d_pos;
  double goal_d_vel;
  double goal_d_acc;
  vector<double> goal_vec;

  
  if ((follow_leader == true) && (lead_s.size() > 0)) {
    goal_s_pos = start_s[0] + lead_s[1] * _num_steps_ahead;
    goal_s_vel = lead_s[1];
  } else {
    goal_s_pos = start_s[0] + _delta_s_with_max_vel;
    goal_s_vel = (_max_vel_in_ts);
  }
  goal_s_acc = 0.0;
  goal_d_pos = (4 * cur_lane_i) + 2;
  goal_d_vel = 0.0;
  goal_d_acc = 0.0;
  goal_vec = {goal_s_pos, goal_s_vel, goal_s_acc, goal_d_pos, goal_d_vel, goal_d_acc};
  //cout << " Eval goal for keep lane :" << goal_s_pos << " " << goal_s_vel << " " << goal_s_acc << " " << goal_d_pos << " " << goal_d_vel << " " << goal_d_acc;
  vector<vector<double>> all_goals_straight = {goal_vec};
  perturb_goal(goal_vec, all_goals_straight);
  // add to goal points
  all_goals.reserve(all_goals.size() + all_goals_straight.size());
  all_goals.insert(all_goals.end(), all_goals_straight.begin(), all_goals_straight.end());
}

void PTG::realize_hard_stop_state_trajectories(int cur_lane_i,
                                               bool follow_leader,
                                               vector<double> const &start_s,
                                               vector<double> const &lead_s,
                                               vector<Vehicle> const &vehicles,
                                               vector<vector<double>>& all_goals)
{
  double goal_s_pos;
  double goal_s_vel;
  double goal_s_acc;
  double goal_d_pos;
  double goal_d_vel;
  double goal_d_acc;
  vector<double> goal_vec;
  
  if ((follow_leader == true) && (lead_s.size() > 0)) {
    _num_steps_ahead = 120;
    goal_s_pos = start_s[0] + (lead_s[0] - start_s[0]);
    goal_s_vel = lead_s[1]*0.9;
    goal_s_acc = 0.0;
    goal_d_pos = ((4 * cur_lane_i) + 2);
    goal_d_vel = 0.0;
    goal_d_acc = 0.0;
    goal_vec = {goal_s_pos, goal_s_vel, goal_s_acc, goal_d_pos, goal_d_vel, goal_d_acc};
    vector<vector<double>> all_goals_hard_stop_ahead = {goal_vec};
    cout << " Eval goal for Hard Stop ahead :" << goal_s_pos << " " << goal_s_vel << " " << goal_s_acc << " " << goal_d_pos << " " << goal_d_vel << " " << goal_d_acc;
    all_goals.reserve(all_goals.size() + all_goals_hard_stop_ahead.size());
    all_goals.insert(all_goals.end(),all_goals_hard_stop_ahead.begin(),all_goals_hard_stop_ahead.end());
  }
}  

void PTG::realize_lane_change_left_trajectories(int cur_lane_i,
                                                bool follow_leader,
                                                vector<double> const &start_s,
                                                vector<double> const &lead_s,
                                                vector<Vehicle> const &vehicles,
                                                vector<vector<double>>& all_goals)
{
  double goal_s_pos;
  double goal_s_vel;
  double goal_s_acc;
  double goal_d_pos;
  double goal_d_vel;
  double goal_d_acc;
  vector<double> goal_vec;

  // Lane Change Left
  goal_s_pos = start_s[0] + _delta_s_with_max_vel;
  goal_s_vel = _max_vel_in_ts;

    if ((follow_leader == true) && (lead_s.size() > 0)) {
      if ((lead_s[0] - start_s[0]) < kCollisionBufferLength *0.5) {
        goal_s_pos = start_s[0] + lead_s[1] * _num_steps_ahead;
        goal_s_vel = lead_s[1];
      }
    }
  goal_s_acc = 0.0;
  goal_d_pos = ((4 * cur_lane_i) + 2) - 4;
  goal_d_vel = 0.0;
  goal_d_acc = 0.0;
  goal_vec = {goal_s_pos, goal_s_vel, goal_s_acc, goal_d_pos, goal_d_vel, goal_d_acc};
  vector<vector<double>> all_goals_left = {goal_vec};
  perturb_goal(goal_vec, all_goals_left);
  // add to goal points
  //cout << " Eval goal for Left :" << goal_s_pos << " " << goal_s_vel << " " << goal_s_acc << \
  //" " << goal_d_pos << " " << goal_d_vel << " " << goal_d_acc;

  all_goals.reserve(all_goals.size() + all_goals_left.size());
  all_goals.insert(all_goals.end(), all_goals_left.begin(), all_goals_left.end());
}

void PTG::realize_lane_change_right_trajectories(int cur_lane_i,
                                                 bool follow_leader,
                                                 vector<double> const &start_s,
                                                 vector<double> const &lead_s,
                                                 vector<Vehicle> const &vehicles,
                                                 vector<vector<double>>& all_goals)
{
  double goal_s_pos;
  double goal_s_vel;
  double goal_s_acc;
  double goal_d_pos;
  double goal_d_vel;
  double goal_d_acc;
  vector<double> goal_vec;

  // Lane Change Right
  goal_s_pos = start_s[0] + _delta_s_with_max_vel;
  goal_s_vel = _max_vel_in_ts;
  
    if ((follow_leader == true) && (lead_s.size() > 0)) {
      if ((lead_s[0] - start_s[0]) < kCollisionBufferLength * 0.5) {
        goal_s_pos = start_s[0] + lead_s[1] * _num_steps_ahead;
        goal_s_vel = lead_s[1];
      }
  }
  goal_s_acc = 0.0;
  goal_d_pos = ((4 * cur_lane_i) + 2) + 4;
  goal_d_vel = 0.0;
  goal_d_acc = 0.0;
  goal_vec = {goal_s_pos, goal_s_vel, goal_s_acc, goal_d_pos, goal_d_vel, goal_d_acc};
  vector<vector<double>> all_goals_right = {goal_vec};
  perturb_goal(goal_vec, all_goals_right);
  // add to goal points
  //cout << " Eval goal for Right :" << goal_s_pos << " " << goal_s_vel << " " << goal_s_acc << \
  //" " << goal_d_pos << " " << goal_d_vel << " " << goal_d_acc;
   
  all_goals.reserve(all_goals.size() + all_goals_right.size());
  all_goals.insert(all_goals.end(),all_goals_right.begin(),all_goals_right.end());
}

void PTG::generate_next_state_trajectories(vector<double> const &start,
                                           vector<Vehicle> const &vehicles,
                                           vector<vector<double>>& all_goals)
{
  vector <State> next_states;
  bool follow_leader = false;
  bool nudge_to_mid_lane = false;

  // 1. Build up the next states to evaluate based on cars ahead/behind in lane.
  // ----------------------------------------------------------------------------

  _current_state_category = kKeepLane;
  
  const vector<double> start_s = {start[0], start[1], start[2]};
  const vector<double> start_d = {start[3], start[4], start[5]};
  int cur_lane_i = d_to_lane(start_d[0]);
  cout << "  Car is in Lane: [ " << cur_lane_i << " ]" << endl;

  vector<int> closest_vehicles;
  closest_vehicles = find_closest_ahead_vehicle_in_each_lane(start,
                                                             vehicles);

  // There is no car ahead in this lane.
  // So keep going on lane and do not consider changing lanes
  if (closest_vehicles[cur_lane_i] == -1) {
    next_states.push_back(KEEP_LANE);
    cout << "!!! No Closest Cars in current Lane " << endl;
  } else {

    // There is a car ahead, but it further than our lookahead metric
    // So just keep going on lane, do not consider changing lanes
    vector<double> closest_car_s = vehicles[closest_vehicles[cur_lane_i]].get_s();
    double dist_to_car = closest_car_s[0] - start[0];
    if (dist_to_car > kMaxLookAhead/2) {
      cout << "     KEEP LANE: Ahead car at dist " << dist_to_car << " > 100 gap. " << endl;
      next_states.push_back(KEEP_LANE);

    } else {

      // Ahead car is dangerously close, and it is slowing down. Too risky to
      // change lines now. Just keep this lane and slow down to follow leader.
      // Reduce the trajectory to slow down.
      if ((dist_to_car < 0.5 * kCollisionBufferLength) &&
          (closest_car_s[1] < (start_s[1] * 0.8))) {
        follow_leader = true;
        _current_state_category = kKeepLaneHardStopAhead;
        cout << "   HARD STOP AHEAD: Ahead car at dist" << dist_to_car << " our vel " << start_s[1]  \
             << " is larger than their vel " << closest_car_s[1] << endl;
        next_states.push_back(KEEP_LANE_HARD_STOP_AHEAD);

      } else {

        // Ahead car is close but greater than our designate collision buffer. Potentially match
        // speed of leader car ahead.
        if (dist_to_car < kCollisionBufferLength) {
          follow_leader = true;
          cout << "    Adjust speed to follow car ahead: at dist " << dist_to_car << " our vel " << start_s[1] << \
            " their vel " << closest_car_s[1] << endl;
        }

        // Nudge to Middle Lane
        // ---------------------
        // We are on an outer lane with a car ahead. While the middle lane is populated, the opposite
        // lane has lots of room. We want to nudge the car to middle lane *if it is safe*, so that it will transition
        // to the outer lane. This nudging is required because normally we would not change lanes given there are
        // cars already in middle lane.
        //
        // We are on outer lane
        if ((cur_lane_i == 0) || (cur_lane_i == 2)) {
          // and middle lane has a car ahead
          if (closest_vehicles[1] != -1)
            {
              vector<double> closest_car_mid_s = vehicles[closest_vehicles[1]].get_s();
              vector<double> closest_car_opp_s;

              // our current lane has a car ahead close that warrants a lane nudge consideration
              if (abs(closest_car_s[0] - start_s[0]) <  2 * kCollisionBufferLength) {

                // car in middle lane also has a car ahead close but
                // there is a safety buffer both front and back
                bool car_in_mid_lane_w_front_room = ((abs(closest_car_mid_s[0] - start_s[0]) <  2 * kCollisionBufferLength) &&
                                                     (abs(closest_car_mid_s[0] - start_s[0]) >  kCollisionBufferLength));
                bool car_in_mid_lane_w_behind_room = true;
                int behind_car_idx = find_closest_behind_vehicle_for_lane(1, start, vehicles);
                if (behind_car_idx != -1) {
                  car_in_mid_lane_w_behind_room = check_behind_vehicle_approach(start, vehicles[behind_car_idx]);
                }
                bool car_in_mid_lane_w_room = car_in_mid_lane_w_front_room && car_in_mid_lane_w_behind_room;
                if (car_in_mid_lane_w_room) {
             
                  // Now check the opposite lane. If opposite lane is empty or opposite
                  // lane has plenty of room to transition to, then do a lane change right or left
                  // depending on which outer lane the car is in.
                  bool no_car_in_opp_lane = (closest_vehicles[abs(cur_lane_i - 2)] == -1);
                  bool car_in_opp_lane_w_room = false;
                  if (closest_vehicles[abs(cur_lane_i - 2)] != -1) {
                    vector<double> closest_car_opp_s = vehicles[closest_vehicles[abs(cur_lane_i - 2)]].get_s();
                    car_in_opp_lane_w_room = (abs(closest_car_opp_s[0] - start_s[0]) >  4 * kCollisionBufferLength);
                  }
                  if (no_car_in_opp_lane || car_in_opp_lane_w_room) {
                    nudge_to_mid_lane = true;
                    cout << "   Nudging to middle lane ... " << endl;
                    _current_state_category = kNudgeToMiddleLane;
                    if (cur_lane_i == 0)
                      {
                        cout << "     Nudge Right to Middle Lane"  << endl;
                        next_states.push_back(LANE_CHANGE_RIGHT);
                      }
                    else
                      {
                        cout << "     Nudge Left to Middle Lane"  << endl;  
                        next_states.push_back(LANE_CHANGE_LEFT);
                      }
                  }
                }
              }
            }
        }

        // No Nudging to middle lane business
        // ----------------------------------
        if (nudge_to_mid_lane == false) {

          // Evaluate if its just worth it to keep going straight
          cout << "     Evaluate KEEP LANE:(regular)"  << endl;
          next_states.push_back(KEEP_LANE);

          // If you are in the middle lane (with a car ahead of course),
          // and there are no cars in either left lane or right lane, choose
          // to pass on the left.
          bool no_left_lane_car = (closest_vehicles[0] == -1);
          bool no_right_lane_car = (closest_vehicles[2] == -1);
          if ((cur_lane_i == 1) &&  no_left_lane_car &&  no_right_lane_car)
            {
              cout << "     In middle lane, choose to pass on left.. "  << endl;
              next_states.push_back(LANE_CHANGE_LEFT);
            }
          else
            {
              // Evaluate a lane change left, unless of course you
              // are in the left-most lane.
              if (cur_lane_i != 0)
                {
                  cout << "     Evaluate LANE CHANGE LEFT "  << endl;
                  next_states.push_back(LANE_CHANGE_LEFT);
                }

              // Evaluate a lane change right, unless of course you
              // are in the right-most lane.
              if (cur_lane_i != 2)
                {
                  cout << "     Evaluate LANE CHANGE RIGHT "  << endl;
                  next_states.push_back(LANE_CHANGE_RIGHT);
                }
            }
        }
      }
    }
  }

  //2. Compute Goals for Next State Trajectories 
  //-------------------------------------------- 
  assert(next_states.size() != 0);
  vector<double> lead_s;
  if (closest_vehicles[cur_lane_i] != -1) {
    lead_s = vehicles[closest_vehicles[cur_lane_i]].get_s();
  }
  for (auto state : next_states)  {
    switch (state) {
    case KEEP_LANE:
      realize_keep_lane_state_trajectories(cur_lane_i, follow_leader, start_s, lead_s, vehicles, all_goals);
      break;      
    case KEEP_LANE_HARD_STOP_AHEAD:
      realize_hard_stop_state_trajectories(cur_lane_i, follow_leader, start_s, lead_s, vehicles, all_goals);
      break;      
    case LANE_CHANGE_LEFT:
      realize_lane_change_left_trajectories(cur_lane_i, follow_leader, start_s, lead_s, vehicles, all_goals);
      break;
    case LANE_CHANGE_RIGHT:
      realize_lane_change_right_trajectories(cur_lane_i, follow_leader, start_s, lead_s, vehicles, all_goals);
      break;
    default:;
    }
  }
}


string PTG::get_current_state_category()
{
  return _current_state_category;
}


//  Make the car accelerate more smoothly from a complete stop
void PTG::start_gracefully_from_complete_stop(double start_velocity)
{
  
  if (start_velocity < _max_vel_in_ts / 1.5) {
    cout << "START: Trying to make the car accelerate more smoothly " << endl;
    double v_diff = _max_vel_in_ts - start_velocity;
    _max_vel_in_ts -= v_diff * 0.7;
    _delta_s_with_max_vel = _num_steps_ahead * _max_vel_in_ts;
  }
}

vector<vector<double>> PTG::generate_sd_path(vector<double> const &start,
                                             double max_speed,
                                             double horizon,
                                             vector<Vehicle> const &vehicles)
{
  const vector<double> start_s = {start[0], start[1], start[2]};
  const vector<double> start_d = {start[3], start[4], start[5]};
  
  cout << " GENERATE SD PATH " << endl;

  // Set up parameters
  _num_steps_ahead = horizon;

  // (miles_to_metres * timestep_in_sec /(hr_to_seconds)) or (1609.34*0.02/3600)
  _max_vel_in_ts = 0.00894 * max_speed * 0.95;
   _delta_s_with_max_vel = _num_steps_ahead * _max_vel_in_ts;

  vector<vector<double>> all_goals;
  vector<pair<vector<double>, vector<double>>> traj_coeffs;
  vector<vector<double>> traj_goals;
  vector<double> traj_costs;

  cout << "start s: " << start_s[0] << " s_vel: " << start_s[1] <<  " s_acc: " << start_s[2] << " d: " << start_d[0] << " d_vec: " << start_d[1]  << " d_acc: " << start_d[2] << endl;

  // Avoid Max Acceleration errors
  start_gracefully_from_complete_stop(start_s[1]);

  // Generate goals for next state trajectories
  generate_next_state_trajectories(start, vehicles, all_goals);

  // Compute the trajectory coefficients for each goal
  for (vector<double> goal : all_goals) {
    vector<double> goal_s = {goal[0], goal[1], goal[2]};
    vector<double> goal_d = {goal[3], goal[4], goal[5]};
    if ((goal[3] > 1.0) && (goal[3] < 11.0)) {
      vector<double> traj_s_coeffs, traj_d_coeffs;  
      traj_s_coeffs = jmt(start_s, goal_s, _num_steps_ahead);
      traj_d_coeffs = jmt(start_d, goal_d, _num_steps_ahead);
      traj_coeffs.push_back(std::make_pair(traj_s_coeffs, traj_d_coeffs));
      traj_goals.push_back({goal[0], goal[1], goal[2], goal[3], goal[4], goal[5]});
    }
  }

  // Evaluate the cost for each trajectory
  vector<vector<double>> all_costs;
  //cout << endl;
  for (int i = 0; i < traj_goals.size(); i++) {
   
      double cost = calculate_cost(traj_coeffs[i], traj_goals[i], vehicles, all_costs);
      //  cout << "        ---> Cost : "<< cost << " : for i=" << i << " " << traj_goals[i][0] << " " << traj_goals[i][1] << " " << traj_goals[i][2] << " " << traj_goals[i][3]<< " " << traj_goals[i][4] << " " << traj_goals[i][5] << endl;
      traj_costs.push_back(cost);
  }

  // Compute the Minimum Cost trajectory
  double min_cost = traj_costs[0];
  int min_cost_idx = 0;
  for (int i = 0; i < traj_costs.size(); i++)
    {
      if (traj_costs[i] < min_cost) {
        min_cost = traj_costs[i];
        min_cost_idx = i;
      }
    }
  cout << "The minimum cost is " << min_cost << " with index " << min_cost_idx << endl;

  // Save the current statate category
  if ((_current_state_category !=  kKeepLaneHardStopAhead) &&
      (_current_state_category !=  kNudgeToMiddleLane))
    {
      if (min_cost_idx > kNumSamplesForPerturb)
        {
          _current_state_category = kLaneChange;
        }
    }
  
  // assert (min_cost != 999999);
  // Sometimes when we are boxed in and there is no optimal trajectory
  // we will get infeasible cost. Still continue on keep lane trajectory
  if (min_cost == 999999) {
    int sample_count  = (traj_costs.size() < kNumSamplesForPerturb)? traj_costs.size() : kNumSamplesForPerturb;
    if (sample_count != 10) {
      cout << " Only one sample for hard stop " << endl;
    }
    double min_s = poly_eval(traj_coeffs[0].first, _num_steps_ahead);
    int min_s_i = 0;
    // find trajectory going straight with minimum s
    for (int i = 1; i < sample_count; i++) {
      if (poly_eval(traj_coeffs[i].first, _num_steps_ahead) < min_s) {
        min_s = poly_eval(traj_coeffs[i].first, _num_steps_ahead);
        min_s_i = i;
      }
    }
    min_cost_idx = min_s_i;
  }
  cout << "Goals Target: "<< traj_goals[min_cost_idx][0] << " " << traj_goals[min_cost_idx][1] << " " << traj_goals[min_cost_idx][2] << " " << traj_goals[min_cost_idx][3]<< " " << traj_goals[min_cost_idx][4] << endl;

  cout << "Cost: " << traj_costs[min_cost_idx] << " minimum index: " << min_cost_idx << endl;
  cout << "------" << endl;

  cout << " Lane Selection Cost: " << all_costs[min_cost_idx][0] << endl;
  cout << " Traffic buffer Cost: " << all_costs[min_cost_idx][1] << endl;
  cout << " Efficiency cost:     " << all_costs[min_cost_idx][2] << endl;
  cout << " Total Acc s cost:    " << all_costs[min_cost_idx][3] << endl;
  cout << " Total acc d cost:    " << all_costs[min_cost_idx][4] << endl;

  // Compute the trajectories in s and d based on the optimal coefficients
  vector<double> traj_s(_num_steps_ahead);
  vector<double> traj_d(_num_steps_ahead);
  
  //cout << "Best Traj s and d are " << endl;
  for (int t = 0; t < _num_steps_ahead; t++) {
    traj_s[t] = poly_eval(traj_coeffs[min_cost_idx].first, t);
    traj_d[t] = poly_eval(traj_coeffs[min_cost_idx].second, t);
    //cout << " (" << traj_s[t] << " , " << traj_d[t] << "),  ";
  }
  //cout << endl;

  vector<vector<double>> new_sd_traj(2);
  new_sd_traj[0] = traj_s;
  new_sd_traj[1] = traj_d;
  return new_sd_traj;
}
