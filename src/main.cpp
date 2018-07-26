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
#include "spline.h" // Added to use spline.h

using namespace std;

// for convenience
using json = nlohmann::json;

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
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
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

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

  double angle = abs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
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
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

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

  // Reference velocity in mph 
  double ref_vel = 0.0; 

  // Starting the drive in the middle lane. Starting value: 
  int lane = 1;

  h.onMessage([&ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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

          int prev_path_size = previous_path_x.size();

          if (prev_path_size > 0) {
            car_s = end_path_s; // This aids in preventing collisions 
          }

          // Prediction: Let's predict ! We begin by setting these boolean values to false 
          bool car_on_front = false;
          bool car_on_left = false;
          bool car_on_right = false;
          
          // Starting values for lane start and end. We begin with multiples of 4 
          int start_of_first_lane = 0;
          int end_of_first_lane = 4;
          int end_of_second_lane = 8;
          int end_of_third_lane = 12;

          int thirty_ahead = 30;
          
          for ( int itr = 0; itr < sensor_fusion.size(); itr++ ) {
           
            // Find out if car is in my lane 
            float f = sensor_fusion[itr][6];
            int car_lane = -1;

            if ( f > start_of_first_lane && f < end_of_first_lane ) {
              car_lane = 0;
            } 
            else if ( f > end_of_first_lane && f < end_of_second_lane ) {
              car_lane = 1;
            } 
            else if ( f > end_of_second_lane && f < end_of_third_lane ) {
              car_lane = 2;
            }

            // Car is not driving on the road
            if (car_lane < 0) {
              continue;
            }

            // Find out car's speed !
            double vx = sensor_fusion[itr][3];
            double vy = sensor_fusion[itr][4];
            double check_car_speed = sqrt(vx * vx + vy * vy);
            double check_car_s = sensor_fusion[itr][5];

            check_car_s += ((double)prev_path_size*0.02*check_car_speed);

            if ( car_lane == lane ) { // This implies the car is in the same lane
              car_on_front |= check_car_s > car_s && check_car_s - car_s < thirty_ahead;
            } else if ( car_lane - lane == -1 ) { // This implies the car is on the left
              car_on_left |= car_s - thirty_ahead < check_car_s && car_s + thirty_ahead > check_car_s;
            } else if ( car_lane - lane == 1 ) { // This implies the car is on the right
              car_on_right |= car_s - thirty_ahead < check_car_s && car_s + thirty_ahead > check_car_s;
            }
          }
          // We stop predicting here 

          // Code for Behavior: How our car behaves in a situation ?

          double diff_speed = 0;
          const double max_speed = 49.9;
          const double max_acceleration = .225;
          if ( car_on_front ) { // If there is a car ahead of us
            if ( !car_on_left && lane > 0 ) {
              // This implies that there is no car is on the left and there exists a left lane.
              lane--; // Shift lanes, move to the left lane
            } else if ( !car_on_right && lane != 2 ){
              // This implies that there is no car is on the right and there exists a right lane.
              lane++; // // Shift lanes, move to the right lane
            } else {
              diff_speed -= max_acceleration;
            }
          } else {
            if ( lane != 1 ) { 
              if ( ( lane == 0 && !car_on_right ) || ( lane == 2 && !car_on_left ) ) {
                lane = 1; 
                // We make the car drive in the center lane if it isn't already doing that
              }
            }
            if ( ref_vel < max_speed ) {
              diff_speed += max_acceleration;
            }
          }

          // Car Trajectory

          // (x,y) waypoints evenly spaced at 30m. 
          // Later, these waypoints are interpolated with spine to fill more points that control speed 
          vector<double> ptsx;
          vector<double> ptsy;

          // Reference for x,y and yaw states
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          if ( prev_path_size < 2 ) {
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);
            // These two points above are tangent to the car

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);

            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          } else {
            // Redefining reference state as previous path end point
            ref_x = previous_path_x[prev_path_size-1];
            ref_y = previous_path_y[prev_path_size-1];

            double ref_x_prev = previous_path_x[prev_path_size-2];
            double ref_y_prev = previous_path_y[prev_path_size-2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            // Using the 2 points which make the path tangent to the previous path's end point
            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }

          // Add in frenet, evenly 30m spaced points ahead of the starting reference
          vector<double> next_wp0 = getXY(car_s + 30, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s + 60, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s + 90, 2 + 4*lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);
          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          // Converting coordinates to car's local coordinates.
          for ( int i = 0; i < ptsx.size(); i++ ) {

            // shifting car's reference angle to ZERO degrees
            double shift_x = ptsx[i] - ref_x;
            double shift_y = ptsy[i] - ref_y;

            ptsx[i] = shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw);
            ptsy[i] = shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw);
          }

          // Creating a spline.
          tk::spline s;

          // Define the actual (x,y) points we will use for the planner
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          // set (x,y) points to the spline
          s.set_points(ptsx, ptsy);
          
          // Start with all of the previous path points from last time
          for ( int i = 0; i < prev_path_size; i++ ) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // Break up spline points so that we have our desired reference velocity
          double target_x = 30.0;
          double target_y = s(target_x);
          double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y));

          double x_add_on = 0;

          // Fill up rest of path planner after filling it with previous points always outputting the 50 points
          for( int i = 1; i < 50 - prev_path_size; i++ ) {
            ref_vel += diff_speed;
            if ( ref_vel > max_speed ) {
              ref_vel = max_speed;
            } else if ( ref_vel < max_acceleration ) {
              ref_vel = max_acceleration;
            }
            double N = target_dist/(0.02 * ref_vel / 2.24);
            double x_point_val = x_add_on + target_x/N;
            double y_point_val = s(x_point_val);

            x_add_on = x_point_val;

            double x_ref_val = x_point_val;
            double y_ref_val = y_point_val;

            // Rotate back to normal after rotating it earlier
            x_point_val = (x_ref_val*cos(ref_yaw)-y_ref_val*sin(ref_yaw));
            y_point_val = (x_ref_val*sin(ref_yaw)+y_ref_val*cos(ref_yaw));

            x_point_val += ref_x;
            y_point_val += ref_y;

            next_x_vals.push_back(x_point_val);
            next_y_vals.push_back(y_point_val);
          }

          json msgJson;

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
