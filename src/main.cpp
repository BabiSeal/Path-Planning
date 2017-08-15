#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "PTG.h"
#include "spline.h"
#include "Vehicle.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Constants

static const int kLocalSplineSegPrev = 9;
static const int kNumPreviousPointToReuse = 15;
static const int kLocalSplineSegAfter = 15;
static const int kMaxNumStepsAhead = 175;
static const int kUpdateInterval = 40;
static const double max_s=6945.554;                    // The max s value before wrapping around the track back to 0

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transforms From Spline smoothed Frenet coord s, d to global Cartesian x, y
vector <double> getXY_splines(double s, double d,
                              tk::spline const &spline_s_to_x,
                              tk::spline const &spline_s_to_y,
                              tk::spline const &spline_s_to_dx,
                              tk::spline const &spline_s_to_dy)
{
  double x_edge = spline_s_to_x(s);
  double y_edge = spline_s_to_y(s);
  double dx = spline_s_to_dx(s);
  double dy = spline_s_to_dy(s);

  double x = x_edge + dx * d;
  double y = y_edge + dy * d;

  return {x, y};
}

                                 
// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while((prev_wp < (int)(maps_s.size()-1) ) && (s > maps_s[prev_wp+1]))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


// Compute the indices of the map to use for the local spline segment
// Uses few points behind given s and points ahead of the given s
// It gets tricky since the indices wrap around.
void getSplineSegmentIndices(double s, vector<double> map_s, vector<double> &indices)
{
  int prev_wp = -1;

  while((prev_wp < (int)(map_s.size()-1)) && s > map_s[prev_wp+1])
      prev_wp++;
 
  // kLocalSplineSegPrev indices before the closest index
  for (int idx = kLocalSplineSegPrev; idx > 0; idx--)
    {
      if (prev_wp - idx < 0)
        indices.push_back(map_s.size() + (prev_wp - idx));
      else
        indices.push_back((prev_wp - idx) % map_s.size());
    }

  // the closest point
  indices.push_back(prev_wp);

  // kLocalSplineSegAfter after the closest index
  for (int idx = 1; idx < kLocalSplineSegAfter; idx++)
    indices.push_back((prev_wp + idx) % map_s.size());
}

//
// Returns the localized spline segments that map s to x,y,dx and dy
//
// First we compute the map-indices around the car's location (few points behind,
// and some ahead). Note the indices might wrap around.
//
// Then we compute local s's - (the offsets from the first s of the global waypoints).
// Again it gets tricky since the points wrap around.
//
// Use the spine function then to fit the local s's to corresponding x,y,dx,dy
//
// Note we need not only to fit x and y but also dx and dy. 
//
void getLocalizedSplineSegments(double s, vector<double> const &map_s,
                                vector<double> const &map_x, vector<double> const &map_y,
                                vector<double> const &map_dx, vector<double> const &map_dy,
                                vector<double> &wp_seg_s_local, vector<double> &wp_seg_s_global,
                                tk::spline &spline_s_to_x,
                                tk::spline &spline_s_to_y,
                                tk::spline &spline_s_to_dx,
                                tk::spline &spline_s_to_dy)
{
  vector<double> wp_indices;
  getSplineSegmentIndices(s, map_s, wp_indices);
  
  vector<double> wp_seg_x, wp_seg_y, wp_seg_dx, wp_seg_dy;
  double start_s = map_s[wp_indices[0]];
  bool passed_zero = false;
  for (int idx = 0; idx < wp_indices.size(); idx++) {
      int wp_i = wp_indices[idx];
      wp_seg_x.push_back(map_x[wp_i]);
      wp_seg_y.push_back(map_y[wp_i]);
      wp_seg_dx.push_back(map_dx[wp_i]);
      wp_seg_dy.push_back(map_dy[wp_i]);
      wp_seg_s_global.push_back(map_s[wp_i]);

      if ((idx > 0) && (wp_i < wp_indices[idx - 1])) {
        passed_zero = true;
      }
      if (passed_zero)
        wp_seg_s_local.push_back(abs(max_s - start_s) + map_s[wp_i]);
      else
        wp_seg_s_local.push_back(map_s[wp_i] - start_s);                         
    }

  //cout << "Indices are " << endl;
  //for (int idx = 0; idx < wp_indices.size(); idx++) {
  //  cout << " " << wp_indices[idx];
  //}
  //cout << endl;
  //cout << "Seg s locals  are " << endl;
  //for (int idx = 0; idx < wp_indices.size(); idx++) {
  //  cout << " " << wp_seg_s_local[idx];
  //}
  //cout << endl;
  //cout << "Map s globals are " << endl;
  //for (int idx = 0; idx < wp_indices.size(); idx++) {
  //  cout << " " << map_s[wp_indices[idx]];
  //}
  //cout << endl;
    
  // Fit the respective splines for x, y, dx and dy
  // to the local s (s offset from beginning of local spline segment)
  spline_s_to_x.set_points(wp_seg_s_local, wp_seg_x);
  spline_s_to_y.set_points(wp_seg_s_local, wp_seg_y);
  spline_s_to_dx.set_points(wp_seg_s_local, wp_seg_dx);
  spline_s_to_dy.set_points(wp_seg_s_local, wp_seg_dy);
}

// Returns a path-vector of x,y coordinates, from a path of s,d coordinates
//
// The inputs are the new-path's s,d coordinates and the previous path's x,y
// in addition to the fitted splines of the reference map.
//
// Instead of directly using the spline function to compute corresponding x,y
// introduce additional smoothing
//
// Based on guidance by Marcus Erbar on discussions.udacity.com and on slack channel #p-path-planning
// https://discussions.udacity.com/t/latency-handling/322156/18
//
// Basically first re-use previous path points (around 15) to populate the new-path.
// Start the new-path from the first element in previous_path. Then subsequent points
// add the difference the spline computed points to the previous point
//
// new_path[0] = previous_path[0]
// new_path[1] = new_path[0] + (spline_planned_path[1] - spline_planned_path[0])
// new_path[2] = new_path[1] + (spline_planned_path[2] - spline_planned_path[1])
//
//
vector<vector<double>>
generate_xy_from_sd(vector<vector<double>> const &new_path, int previous_path_size,
                    vector<double> const & previous_path_x,
                    vector<double> const & previous_path_y,
                    tk::spline const &spline_s_to_x,
                    tk::spline const &spline_s_to_y,
                    tk::spline const &spline_s_to_dx,
                    tk::spline const &spline_s_to_dy)
{

  bool first_cycle = (previous_path_size == 0);
  int  use_previous_points = kNumPreviousPointToReuse;
  vector<double> x_vals;
  vector<double> y_vals;
  vector <double> prev_xy;
  vector <double> xy;
  double new_x, new_y;

  // Reuse previous path if it is there
  for (int i = 0; i < use_previous_points; i++) {
    prev_xy = getXY_splines(new_path[0][i], new_path[1][i],
                            spline_s_to_x, spline_s_to_y,
                            spline_s_to_dx, spline_s_to_dy);
    
    // No previous points generate paths
    if (first_cycle) {
      x_vals.push_back(prev_xy[0]);
      y_vals.push_back(prev_xy[1]);
    } else {
      // Reuse previous points
      new_x = previous_path_x[i];
      new_y = previous_path_y[i];
      x_vals.push_back(new_x);
      y_vals.push_back(new_y);
    }
  }

  // Generate new points 
  for (int i = use_previous_points; i < new_path[0].size(); i++) {
    xy = getXY_splines(new_path[0][i], new_path[1][i],
                       spline_s_to_x, spline_s_to_y,
                       spline_s_to_dx, spline_s_to_dy);
    if (first_cycle) {
      x_vals.push_back(xy[0]);
      y_vals.push_back(xy[1]);
    } else {
      double x_dif =  xy[0] - prev_xy[0];
      double y_dif =  xy[1] - prev_xy[1];
      new_x = new_x + x_dif;
      new_y = new_y + y_dif;
      x_vals.push_back(new_x);
      y_vals.push_back(new_y);
      prev_xy = xy;
    }
  }
  
  vector<vector<double>> new_xy_vals(2);
  new_xy_vals[0] = x_vals;
  new_xy_vals[1] = y_vals;
  return new_xy_vals;
}

  // converts world space s coordinate to local space based on provided mapping
double get_local_s(double world_s,
                   vector<double> const &waypoints_segment_s_worldSpace,
                   vector<double> const &waypoints_segment_s)
{
  int prev_wp = 0;
  // special case: first wp in list is larger than s. Meaning we are crossing over 0 somewhere.
  // go to index with value zero first and search from there.
  if (waypoints_segment_s_worldSpace[0] > world_s) {
      while (waypoints_segment_s_worldSpace[prev_wp] != 0.0)
          prev_wp += 1;
  }
  while ((waypoints_segment_s_worldSpace[prev_wp+1] < world_s) &&
         (waypoints_segment_s_worldSpace[prev_wp+1] != 0))
      prev_wp += 1;
  double diff_world = world_s - waypoints_segment_s_worldSpace[prev_wp];
  double s_local = waypoints_segment_s[prev_wp] + diff_world;
  return s_local;
}



int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  Vehicle ego; //Our Car
  PTG ptg;
  
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &ego, &ptg](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
                
                
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
                int num_steps_ahead = kMaxNumStepsAhead;
                int update_interval = kUpdateInterval;

                // Just use the previous path's to update trajectory
                if (previous_path_x.size() >= (num_steps_ahead - update_interval)) {
                  //cout << "USE OLD TRAJECTORY prev path size: " <<  previous_path_x.size() \
                  //      << " car_s : " << car_s << " car_d : " << car_d << endl;
                  for(int i = 0; i < previous_path_x.size(); i++) {
                    next_x_vals.push_back(previous_path_x[i]);
                    next_y_vals.push_back(previous_path_y[i]);
                  }
                } else {
                  // Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
                  
                  ego.set_s(car_s, 0, 0);
                  ego.set_d(car_d, 0, 0);

                  // Localization and Spline Fitting
                  //---------------------------------
                  // We are generating a set of local splines for x,y,dx,dy trajectories.
                  //
                  // We first select some map waypoints close to the location of the car. As per Marcus E's
                  // recommendation we select 9 waypoints behind and 15 waypoints ahead.
                  // We compute a vector of local s correponding to the waypoints selected. This local "s"
                  // starts at 0 for the first point s waypoint selected previously. Subsequent local s's
                  // are the offsets of the other s waypoints from this first s waypoint.
                  //
                  // The splines for corresponding x,y,dx,dy are are then computed as fitting them on the basis
                  // of these "offset" local s's.
                  //
                  // Note the local_s the car is no zero since it is not the reference point, but its s is also
                  // localized on the basis of the first local s. The sensor fusion s parameters are correspondingly
                  // localized relative to that first local s point too.
                  tk::spline spline_s_to_x;
                  tk::spline spline_s_to_y;
                  tk::spline spline_s_to_dx;
                  tk::spline spline_s_to_dy;
                  vector<double> wp_seg_s_local, wp_seg_s_global;
                  getLocalizedSplineSegments(car_s, map_waypoints_s,
                                             map_waypoints_x, map_waypoints_y,
                                             map_waypoints_dx, map_waypoints_dy,
                                             wp_seg_s_local, wp_seg_s_global,
                                             spline_s_to_x, spline_s_to_y,
                                             spline_s_to_dx, spline_s_to_dy);
                

                  // Localization: convert current car_s into our local Frenet space
                  double car_local_s = get_local_s(car_s, wp_seg_s_global, wp_seg_s_local);
                  
                  // Localization: convert sensor fusion data into local Frenet space
                  for (int i = 0; i < sensor_fusion.size(); i++) {
                    sensor_fusion[i][5] = get_local_s(sensor_fusion[i][5], wp_seg_s_global, wp_seg_s_local);
                  }
                
                  // get previous forecasted car state from previous path computations
                  // guidance from marcus on path-planning channel
                  // 
                  int lag = num_steps_ahead - update_interval - previous_path_x.size();
                  cout << "PATH UPDATE" << endl;
                  cout << " prev path size: " <<  previous_path_x.size() << " Update Interval " << update_interval  \
                       << " Global car_s : " << car_s << " Global car_d : " << car_d << " lag : " << lag << endl;
                  if (lag >= 10) lag = 0;  // this is the first run of simulation

                  // get last known car forecasted state
                  // collect best guess at current car state. S position in local segment space
                  double est_car_s_vel = ego.get_forecasted_param(lag, 0);
                  double est_car_s_acc = ego.get_forecasted_param(lag, 1);
                  double est_car_d_vel = ego.get_forecasted_param(lag, 2);
                  double est_car_d_acc = 0.0;
                  vector<double> car_state = {car_local_s, est_car_s_vel , est_car_s_acc, car_d, est_car_d_vel, est_car_d_acc};
                  cout << " Predicted car state vel: " <<  est_car_s_vel << " " << est_car_s_acc << " " << est_car_d_vel << " " << est_car_d_acc << " " << endl;
                  
                  // turn sensor fusion data into Vehicle objects
                  vector<Vehicle> other_cars(sensor_fusion.size());
                  for (int i = 0; i < sensor_fusion.size(); i++) {

                    // Other car has constant velocity and no acceleration in s
                    double other_car_s = sensor_fusion[i][5];
                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];

                    // Note the deliberate localization to unit-time being 0.02 sec and not 1 sec
                    // This helps us in our calculation of velocity.
                    double other_car_s_vel_per_ts = sqrt(pow(vx, 2) + pow(vy, 2)) / 50.0;
                    other_cars[i].set_s(other_car_s, other_car_s_vel_per_ts, 0.0);
                  
                    // Other car has no velocity in d and no acceleration in d
                    double other_car_d = sensor_fusion[i][6];
                    other_cars[i].set_d(other_car_d, 0.0, 0.0);
                  }

                  // Poor hack - car in left-most lane can go fastest, outer lanes with more
                  // curvature - need to reduce the speed limit to not exceed the speed limit
                  double speed_limit;
                  int car_lane = ptg.d_to_lane(car_d);
                  if (car_lane == 0) 
                    speed_limit = 47.35;
                  else if (car_lane == 1)
                    speed_limit = 47.25;
                  else if (car_lane == 2)
                    speed_limit = 47.15;

                  // -----------------------------------------------------
                  // MAIN - Path Generation Function invocation here
                  // This will invoke the Polynomial Trajectory generator
                  // to evaluate the optimal next states, generate goal
                  // points for prospective next next state trajectories,
                  // evaluate them using the cost functions and then return
                  // the minimum cost trajectory.
                  // ---------------------------------------------------
                  vector<vector<double>> new_path = ptg.generate_sd_path(car_state,
                                                                         speed_limit,
                                                                         num_steps_ahead,
                                                                         other_cars);
                   if (ptg.get_current_state_category() == ptg.kKeepLaneHardStopAhead) {
                     cout << "Result: SLOW DOWN !!!" << endl;
                     // Need to reduce the num_steps_ahead and update interval to enable fast stopping
                     num_steps_ahead = 120;
                     update_interval = num_steps_ahead - 80;
                   }
                   if (ptg.get_current_state_category() == ptg.kLaneChange) {
                     cout << "Result: Do a Lane Change..." << endl;
                   }
                   if (ptg.get_current_state_category() == ptg.kNudgeToMiddleLane) {
                      cout << "Result: Nudging to Middle Lane..." << endl;
                   } else {
                     if (ptg.get_current_state_category() == ptg.kKeepLane) {
                       cout << "Result: Stay in Lane" << endl;
                     }
                   }
                
                   // Save forecasted points for next round usage
                   for (int i = 0; i < 10; i++) {
                    double s_start = new_path[0][i + update_interval];
                    double s_first = new_path[0][i + update_interval + 1];
                    double s_second = new_path[0][i + update_interval + 2];
                    double d_start = new_path[1][i + update_interval];
                    double d_first = new_path[1][i + update_interval + 1];
                    double d_second = new_path[1][i + update_interval + 2];
                    double s_vel_first = s_first - s_start;
                    double s_vel_second = s_second - s_first;
                    double s_acc = s_vel_second - s_vel_first;
                    double d_vel_first = d_first - d_start;
                    double d_vel_second = d_second- d_first;
                    double d_acc = d_vel_second - d_vel_first;
                    ego.set_forecasted_param(i, s_vel_first, s_acc, d_vel_first, d_acc);
                    //cout << " For next round forecast (i) : " << i+update_interval <<  " " << s_vel_first << " " << s_acc << " " << d_vel_first << " " << d_acc << endl;
                    
                  }
                   // Convert the s,d coordinates of the new path to xy coordinates.
                   vector<vector<double>> new_xy_vals = generate_xy_from_sd(new_path, previous_path_x.size(),
                                                                            previous_path_x, previous_path_y,
                                                                            spline_s_to_x, spline_s_to_y,
                                                                            spline_s_to_dx, spline_s_to_dy);
                  next_x_vals = new_xy_vals[0];
                  next_y_vals = new_xy_vals[1];
                }
                msgJson["next_x"] = next_x_vals;
                msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































