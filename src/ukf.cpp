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
UKF::UKF(){
  // if this is false, laser measurements will be ignored (except during init)
  is_initialized_ = false;
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_ ;

  // initial state vector
  x_ = VectorXd(n_x_);
  Xsig_pred_ = Eigen::MatrixXd(n_x_,2*n_aug_ +1);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
 // P_.fill(0);
 // P_(0,0) = 1;
 // P_(1,1) = 1;
 // P_(2,2) = 1;
 // P_(3,3) = 1;
 // P_(4,4) = 1;
P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.78;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.65; 

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
  R_radar(0,0) = std_radr_ * std_radr_;
  R_radar(1,1) = std_radphi_ * std_radphi_;
  R_radar(2,2) = std_radrd_ * std_radrd_;
  R_lidar(0,0) = std_laspx_ * std_laspx_;
  R_lidar(1,1) = std_laspy_ * std_laspy_;

  weights_ = VectorXd(2 * n_aug_ +1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i=1 ; i < 2*n_aug_ + 1; i++){
      weights_(i) = 0.5/(lambda_ + n_aug_);
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
 
  if(! is_initialized_){
      if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true){
          
          double rho = meas_package.raw_measurements_[0];
          double phi = meas_package.raw_measurements_[1];

          x_(0) = rho * cos(phi);
          x_(1) = rho * sin(phi);
          x_(2) = 0;
          x_(3) = 0;
          x_(4) = 0;
      }

      else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true){
          
          double px = meas_package.raw_measurements_[0];
          double py = meas_package.raw_measurements_[1];
          
          if(fabs(px) < 1e-3)
              px = 1e-5;
          if(fabs(py) < 1e-3)
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
    
   double dt = (meas_package.timestamp_ - time_us_) * 1e-6;
   if( use_radar_ || use_laser_ ){
       
       //cout<<"Started Prediction"<<endl;
       Prediction(dt);
       //cout<<"Predicted state<<"<<endl<<x_<<endl;
       
       if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
           UpdateRadar(meas_package);
           //cout<<"After Radar Measurement state<<"<<endl<<x_<<endl;
       }
        if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
           UpdateLidar(meas_package);
             //cout<<"After Lidar Measurement state<<"<<endl<<x_<<endl;
       }
   }
    
   time_us_ = meas_package.timestamp_;
   return;
 
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
  MatrixXd sig_points = MatrixXd(n_aug_, 2*n_aug_ +1);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_ + 1) = 0;
  generateSigmaPoints(x_aug , sig_points);
    
  //cout<<"Genereated Sigma Points<<"<<endl<<sig_points<<endl;    
  predictSigmaPoints(sig_points , Xsig_pred_ , delta_t);
  //calculate the mean
  x_.fill(0);
  for(int i = 0 ; i < weights_.size() ; i++){
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  //calculate co-variance
  P_.fill(0);
  

  for(int i =0; i < weights_.size(); i++){
      VectorXd diff = Xsig_pred_.col(i)  - x_;
      while(diff(3) > M_PI) diff(3) -= 2.0*M_PI;
      while(diff(3) < -M_PI) diff(3) += 2.0*M_PI;
      P_ += weights_(i) * diff * diff.transpose();
  }
    return;
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
    int n_z = 2;
    MatrixXd z_sig = MatrixXd (n_z , 2*n_aug_ + 1);
    VectorXd z_pred = VectorXd(n_z);
    VectorXd z_meas = meas_package.raw_measurements_;
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd T = MatrixXd(n_x_,n_z);
    MatrixXd K = MatrixXd(n_x_,n_z);
    S.fill(0);
    T.fill(0);
    z_pred.fill(0);
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        z_sig(0,i) = px;
        z_sig(1,i) = py;
    }
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        z_pred += weights_(i) * z_sig.col(i);
    }
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        VectorXd diff = VectorXd(n_z);
        diff = z_sig.col(i) - z_pred;
        S += weights_(i) * diff * diff.transpose();
    }
    
    S += R_lidar;
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        VectorXd diff = Xsig_pred_.col(i) - x_;
        while (diff(3)> M_PI) diff(3)-=2.*M_PI;
        while (diff(3)<-M_PI) diff(3)+=2.*M_PI;
        VectorXd diff_meas = z_sig.col(i) - z_pred;
        
        T += weights_(i) * diff * diff_meas.transpose();
      
    }
    
      K = T * S.inverse();
      x_ = x_ + K*(z_meas -z_pred);
      P_ = P_ - K*S*K.transpose();
      
     //Calculate NIS 
     VectorXd diff = (z_meas - z_pred);
     NIS_laser_ = diff.transpose() * S.inverse() * diff;  
    
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
    int n_z = 3;
    MatrixXd z_sig = MatrixXd (n_z , 2*n_aug_ + 1);
    VectorXd z_pred = VectorXd(n_z);
    VectorXd z_meas = meas_package.raw_measurements_;
    MatrixXd S = MatrixXd(n_z,n_z);
    MatrixXd T = MatrixXd(n_x_,n_z);
    MatrixXd K = MatrixXd(n_x_,n_z);
    S.fill(0);
    T.fill(0);
    z_pred.fill(0);
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double phi = Xsig_pred_(3,i);
        
        double rho = sqrt(px*px + py*py);
        if(rho < 0.0001)
            return;
        z_sig(0,i) = rho;
        z_sig(1,i) = atan2(py,px);
        z_sig(2,i) = ( px * v * cos(phi) + py * v * sin(phi) ) / rho;
    }
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        z_pred += weights_(i) * z_sig.col(i);
    }
    
    for(int i = 0; i < 2*n_aug_ + 1; i++){
        VectorXd diff = VectorXd(n_z);
        diff = z_sig.col(i) - z_pred;
        
        while(diff(1) > M_PI) diff(1) -= 2.*M_PI;
        while(diff(1) < -M_PI) diff(1) += 2.*M_PI;
        S += weights_(i) * diff * diff.transpose();
    }
    
    S += R_radar;
    
     for(int i = 0; i < 2*n_aug_ + 1; i++){
        VectorXd diff = Xsig_pred_.col(i) - x_;
        VectorXd diff_meas = z_sig.col(i) - z_pred;

        while(diff(3) > M_PI) diff(3) -= 2.*M_PI;
        while(diff(3) < -M_PI) diff(3) += 2.*M_PI;
         
        while(diff_meas(1) > M_PI) diff_meas(1) -= 2.*M_PI;
        while(diff_meas(1) < -M_PI) diff_meas(1) += 2.*M_PI;
         
        T += weights_(i) * diff * diff_meas.transpose();
         
    }
        K = T * S.inverse();
        VectorXd diff = (z_meas -z_pred);
        while(diff(1) > M_PI) diff(1) -= 2.*M_PI;
        while(diff(1) < -M_PI) diff(1) += 2.*M_PI;
        //cout<<"diff: "<<endl<<diff<<endl;
        x_ = x_ + K*diff;
        P_ = P_ - K*S*K.transpose();
        NIS_radar_ = diff.transpose() * S.inverse() * diff;
    
}
void UKF::generateSigmaPoints(const VectorXd x_aug , MatrixXd &sig_points){
    MatrixXd p_aug = MatrixXd(n_aug_ , n_aug_);
    p_aug.fill(0);
    p_aug.topLeftCorner(n_x_,n_x_) = P_;
    p_aug(n_x_,n_x_) = std_a_ * std_a_ ; 
    p_aug(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_ ; 

    MatrixXd L = p_aug.llt().matrixL();

    sig_points.col(0) = x_aug;
    double scale = sqrt(lambda_ + n_aug_);
    for(int i = 0; i<n_aug_; i++){
        sig_points.col(1 + i) = x_aug + scale * L.col(i);
        sig_points.col(1 + i + n_aug_) = x_aug - scale * L.col(i);
    }


}
void UKF::predictSigmaPoints(const MatrixXd& sig_points , MatrixXd& pred_sig_points, double dt){
   
    //cout<<"predicting sigma points"<<endl;
    
    for(int i = 0; i < 2 * n_aug_ + 1; i++){

        double px = sig_points(0,i); 
        double py = sig_points(1,i); 
        double v = sig_points(2,i);
        double phi = sig_points(3,i);
        double phi_dot = sig_points(4,i);
        double nu_a = sig_points(5,i);
        double nu_yaw = sig_points(6,i);

        if(fabs(phi_dot) > 0.001){
            pred_sig_points(0,i) = px + (v/phi_dot) * (sin(phi + phi_dot *dt) - sin(phi) ) + 0.5*dt*dt*cos(phi) * nu_a ;
            pred_sig_points(1,i) = py + (v/phi_dot) * (cos(phi) - cos(phi + phi_dot * dt)) + 0.5*dt*dt*sin(phi) * nu_a;
        }
        else{

            pred_sig_points(0,i) = px + v * cos(phi) * dt + 0.5 * dt * dt * cos(phi) * nu_a;
            pred_sig_points(1,i) = py + v * sin(phi) * dt + 0.5 * dt * dt * sin(phi) * nu_a;
            
        }

        pred_sig_points(2,i) = v + dt * nu_a;
        pred_sig_points(3,i) = phi + dt * phi_dot + 0.5 * dt * dt * nu_yaw;
        pred_sig_points(4,i) = phi_dot + dt * nu_yaw;
        
        //cout<<"phi: "<<phi<<" phi_dot: "<<phi_dot<<" dt :"<<dt<<endl;
    } 
    //cout<<"predicted sigma points"<<endl<<pred_sig_points<<endl;
}
