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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	default_random_engine gen;
	// particles number
        num_particles = 100;
	// Add random Gaussian noise
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	particles.resize(num_particles) ;
	weights.resize(num_particles);
	// Initialize particles
	for(int i =0; i < num_particles; i++){		
		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		particles[i].weight = 1.0f;
		weights[i] = 1.0f;
		//cout << "init particle" << particles[i].x << "," << particles[i].y << "," << particles[i].theta << endl;
	}
	// Set status
	is_initialized = true;
	
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	// Prediction
	for(int i =0; i < num_particles; i++){	
		Particle &p = particles[i];

		double theta_ = p.theta;
		double x_ = p.x;
		double y_ = p.y;
		// avoid division by zero 
		if(fabs(yaw_rate) > 0.000001) {
			theta_ += yaw_rate*delta_t;
			x_ += velocity/yaw_rate*(sin(theta_)-sin(p.theta));
    			y_ -= (velocity/yaw_rate*(cos(theta_)-cos(p.theta)));
		}
		else{
			x_ += velocity * delta_t * cos(p.theta) ;
			y_ += velocity * delta_t * sin(p.theta) ;
		}
		
		// Add random Gaussian noise
		normal_distribution<double> dist_x(x_, std_pos[0]);
		normal_distribution<double> dist_y(y_, std_pos[1]);
		normal_distribution<double> dist_theta(theta_, std_pos[2]);
 		
		// predict all particles
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
		//cout << "predict particle"<< particles[i].x << "," << particles[i].y << "," << particles[i].theta << endl;
		
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
	
	// updateWeights
	weights.clear();
	for(int i =0; i < num_particles; i++){	
	
		//vector<LandmarkObs> particle_obs_list ;

		double particle_theta = particles[i].theta;
      		double cos_theta      = cos(particle_theta) ;
      		double sin_theta      = sin(particle_theta) ; 
		
		// partile weight
		double w_p = 1.0f;
		for(int j = 0; j < observations.size(); j++){
			// observations should be in sensor range
			if (dist(observations[j].x, observations[j].y, 0, 0) < sensor_range) {
				LandmarkObs obs_map;
				//obs_map.id = j;
      				// p_x + (obs_x * cos(theta)) - (obs_y * sin(theta))
      				obs_map.x =  particles[i].x   + (observations[j].x * cos_theta) - (observations[j].y * sin_theta) ;
      				// p_y + (obs_x * sin(theta)) + (obs_y * cos(theta))             
      				obs_map.y =  particles[i].y   + (observations[j].x * sin_theta) + (observations[j].y * cos_theta) ;
				//particle_obs_list.push_back(obs_map);
				//cout << "obs" << obs_map.x << "," << obs_map.y << endl;
				
				// to simply associate the closest landmark
				double dist_ass_landmark = sensor_range;
				int ass_landmark_id = -1;
				for (int l = 0; l < map_landmarks.landmark_list.size(); l++) {
					double dst = dist(obs_map.x,obs_map.y,map_landmarks.landmark_list[l].x_f,map_landmarks.landmark_list[l].y_f);
					if (dst < dist_ass_landmark) {
						dist_ass_landmark = dst;
						ass_landmark_id = l;
					}
				}
				// to get particle each observation updated weight w_o_obs and multipy to particle weight w_p
				if (ass_landmark_id > 0) {
					//cout << "landmark" << map_landmarks.landmark_list[ass_landmark_id].x_f << "," << map_landmarks.landmark_list[ass_landmark_id].y_f << endl;
					// Multivariate-Gaussian probability of observations[j]
					double c1 = ((obs_map.x-map_landmarks.landmark_list[ass_landmark_id].x_f)*(obs_map.x-map_landmarks.landmark_list[ass_landmark_id].x_f))/(2*std_landmark[0]*std_landmark[0]);
					double c2 = ((obs_map.y-map_landmarks.landmark_list[ass_landmark_id].y_f)*(obs_map.y-map_landmarks.landmark_list[ass_landmark_id].y_f))/(2*std_landmark[1]*std_landmark[1]);
					double w_p_obs = (1./(2.*3.1415926*std_landmark[0]*std_landmark[1]))*exp(-1.*(c1+c2));
					
					w_p = w_p*w_p_obs;
					//particles[i].weight = w_p_obs;
					
				}
				//cout << "obs_map" << obs_map.x << "," << obs_map.y << ", weight" << w_p << endl;
			}
			
		}
		particles[i].weight = w_p;
		weights.push_back(w_p);	
		//cout << "i" << i << "weight" << particles[i].weight << endl;	
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	default_random_engine gen;

	// Creates a discrete distribution for weight.
    	std::discrete_distribution<int> distribution_wts(weights.begin(), weights.end());
    	std::vector<Particle> resampled_particles;
    	// Resample
    	for(int i=0;i<num_particles;i++){
		Particle particles_i = particles[distribution_wts(gen)];
        	resampled_particles.push_back(particles_i);
    	}
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
