#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include "math.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  is_initialized_ = false;
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_ ;

  // initial state vector
  x_ = VectorXd(n_x_);
  Xsig_pred_ = Eigen::MatrixXd(n_aug_,2*n_x_ +1);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  R_radar = MatrixXd::Zero(3,3);
  R_lidar = MatrixXd::Zero(2,2);
  R_radar(0,0) = std_radr_;
  R_radar(1,1) = std_radphi_;
  R_radar(2,2) = std_radrd_;
  R_lidar(0,0) = std_laspx_;
  R_lidar(1,1) = std_laspy_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
 
  if(! is_initialized_){
      if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true){
          
          float rho = meas_package.raw_measurements_[0];
          float phi = meas_package.raw_measurements_[1];
          float rhodot = meas_package.raw_measurements_[2];

          x_(0) = rho * cos(phi);
          x_(1) = rho * sin(phi);
          x_(2) = rhodot;
          x_(3) = phi;
          x_(4) = 0;
      }

      else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true){
          
          float px = meas_package.raw_measurements_[0];
          float py = meas_package.raw_measurements_[1];
          
          if(px == 0)
              px = 1e-5;
          if(py == 0)
              py = 1e-5;
          
          x_(0) = px ;
          x_(1) = py;
          x_(2) = 0;
          x_(3) = 0;
          x_(4) = 0;
      }
      is_initialized_ = true;
      time_us_ = meas_package.timestamp_;
      return;
  }
    
   double dt = meas_package.timestamp_ - time_us_;
   if( use_radar_ || use_laser_ ){
       Prediction(dt);
       
       if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
           UpdateRadar(meas_package);
       }
        if(meas_package.sensor_type_ == MeasurementPackage::LASER){
           UpdateLidar(meas_package);
       }
   }
    
   time_us_ = meas_package.timestamp_;
 
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
     //create augmented state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;
  
  //Generate Sigma points 
  //Predict mean and covariance of sigma points 
  //predict state
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
