/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// - Set the number of particles.
  // - Initialize all particles to first position including uncertainties from GPS
	// - Set all weights to 1. 
	// - Add random Gaussian noise to each particle.
  
  num_particles = 100;

  particles.resize(num_particles);
  weights.resize(num_particles);
  
  random_device rd;
  default_random_engine gen(rd());
  normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

  for (int i=0 ; i<num_particles ; ++i ) {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0;
  }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // - Add measurements to each particle and
  // - Add random Gaussian noise - std::normal_distribution, std::default_random_engine useful.
	//    http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//    http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;

  // Normal distributions for sensor noise
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  for (int i=0 ; i<num_particles ; ++i ) {
    double theta = particles[i].theta;
    // Update state
    if (fabs(yaw_rate) < 0.00001) {
      particles[i].x += (velocity * delta_t * cos(theta));
      particles[i].y += (velocity * delta_t * sin(theta));
    }
    else {
      particles[i].x += (velocity/yaw_rate) * (sin(theta + yaw_rate*delta_t) - sin(theta));
      particles[i].y += (velocity/yaw_rate) * (cos(theta) - cos(theta + yaw_rate*delta_t));
      particles[i].theta += yaw_rate*delta_t;
    }
    // Add random Gaussian noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// - Find the predicted measurement that is closest to each observed measurement
  // - Assign the observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  double min_distance;
  double distance;
  int id;

  LandmarkObs obs, pred;

  for (int i=0 ; i<observations.size(); ++i) {
    obs = observations[i];
    min_distance = 0.0;
    for (int j=0 ; j<predicted.size() ; ++j) {
      pred = predicted[j];
      distance = dist(obs.x, obs.y, pred.x, pred.y);
      // Update the closest predicted measurement
      if ( (j==0) || (distance < min_distance) ) {
        min_distance = distance;
        id = pred.id;
      }
    }
    // Decide the closest predicted measurement and store its id
    observations[i].id = id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
  // - Transform the observations in the VEHICLE's coordinate system to MAP's coordinate
  // - Associate transformed observation with a landmark
  // - Update the weights of each particle using a multivariate Gaussian distributtion.

  double gauss_norm = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);

  for (int i=0 ; i<num_particles ; ++i ) {  // Update the weight for each particles
    double x_p = particles[i].x;
    double y_p = particles[i].y;
    double theta_p = particles[i].theta;

    // Transforming landmark observations
    vector<LandmarkObs> observations_t;

    for (int j=0 ; j<observations.size() ; ++j) {
      double x_map = (cos(theta_p) * observations[j].x) - (sin(theta_p) * observations[j].y) + x_p;
      double y_map = (sin(theta_p) * observations[j].x) + (cos(theta_p) * observations[j].y) + y_p;
      observations_t.push_back(LandmarkObs{observations[j].id, x_map, y_map});
    }

    vector<LandmarkObs> predicted;

    // Concern only the observation within the range of sensor
    for (int j=0 ; j<map_landmarks.landmark_list.size(); ++j) {
      double x_l = map_landmarks.landmark_list[j].x_f;
      double y_l = map_landmarks.landmark_list[j].y_f;
      int id_l = map_landmarks.landmark_list[j].id_i;

      if ( dist(x_l, y_l, x_p, y_p) <= sensor_range ) 
        predicted.push_back(LandmarkObs{id_l, x_l, y_l});
    }

    dataAssociation(predicted, observations_t);

    particles[i].weight = 1.0;

    for (int j=0 ; j<observations_t.size(); ++j) {
      double x_obs = observations_t[j].x;
      double y_obs = observations_t[j].y;

      int id_associated = observations_t[j].id; 

      double x_pred;
      double y_pred;
      for (int k=0 ; k<predicted.size() ; ++k) {
        if (predicted[k].id == id_associated) {
          x_pred = predicted[k].x;
          y_pred = predicted[k].y;
        }
      }

      double exponent = pow((x_obs-x_pred),2)/(2*pow(std_landmark[0],2)) + pow((y_obs-y_pred),2)/(2*pow(std_landmark[1],2));
      particles[i].weight *= (gauss_norm * exp(-exponent));
    }
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> resampled_particles;

  // get all the weights for generating distribution
  vector<double> weights;
  for (int i=0 ; i<num_particles ; ++i) {
    weights.push_back(particles[i].weight);
  }

  random_device rd;
  default_random_engine gen(rd());
  for (int i=0 ; i<num_particles ; ++i) {
    discrete_distribution<int> resampling(weights.begin(), weights.end());
    resampled_particles.push_back(particles[resampling(gen)]);
  }

  // Update the resampled particles
  particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
