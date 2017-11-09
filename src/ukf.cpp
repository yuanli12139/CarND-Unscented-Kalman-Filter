#include "ukf.h"
#include "Eigen/Dense"
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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .5; //30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .15; //30;

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
  is_initialized_ = false;

  // Time when the state is true, in us
  time_us_ = 0;
  
  // State dimension in CTRV model
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2*n_aug_ + 1; i++)
    weights_(i) = 1 / (2 * (lambda_ + n_aug_));

  // NIS for lidar
  NIS_laser_ = 0.;

  // NIS for radar
  NIS_radar_ = 0.;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    // Initializion
    time_us_ = meas_package.timestamp_;

    double p_x = meas_package.raw_measurements_[0];
    double p_y = meas_package.raw_measurements_[1];
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      x_ << p_x, p_y, 0., 0., 0.;
    } 
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      x_ << rho * cos(phi), rho * sin(phi), rho_dot * cos(phi), 0., 0.;
    }

    P_.fill(0.0);
    P_(0, 0) = 1.;
    P_(1, 1) = 1.;
    P_(2, 2) = 1.;
    P_(3, 3) = 1.;
    P_(4, 4) = 1.;

    is_initialized_ = true;
  }
  else {

    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.; //time diff in sec

    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      Prediction(delta_t);
      UpdateLidar(meas_package);
      time_us_ = meas_package.timestamp_;
    } 

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      Prediction(delta_t);
      UpdateRadar(meas_package);

      time_us_ = meas_package.timestamp_;
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;

  int i;
  for (i = 1; i < n_aug_ + 1; i++)
    Xsig_aug.col(i) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i - 1);
  for (; i < 2 * n_aug_ + 1; i++)
    Xsig_aug.col(i) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i - n_aug_ - 1);

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > .001) {
        px_p = p_x + v/yawd * (sin(yaw + yawd * delta_t) - sin(yaw)) + .5 * delta_t * delta_t * cos(yaw) * nu_a;
        py_p = p_y + v/yawd * (-cos(yaw + yawd * delta_t) + cos(yaw)) + .5 * delta_t * delta_t * sin(yaw) * nu_a;
    }
    else {
        px_p = p_x + v * cos(yaw) * delta_t + .5 * delta_t * delta_t * cos(yaw) * nu_a;
        py_p = p_y + v * sin(yaw) * delta_t + .5 * delta_t * delta_t * sin(yaw) * nu_a;
    }

    double v_p = v + delta_t * nu_a;
    double yaw_p = yaw + yawd * delta_t + .5 * delta_t * delta_t * nu_yawdd;
    double yawd_p = yawd + delta_t * nu_yawdd;
    
    //write predicted sigma points into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;

  }

  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);

  //predict state corvariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
 
// cout << "Xsig_pred :" << endl << Xsig_pred_.col(i) << endl;
    //angle normalization
    while(x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
    while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();  
  }

  cout << "Predicted x_: " << x_ << endl;
  cout << "Predicted P_: " << P_ << endl;
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
  MatrixXd R = MatrixXd(2, 2);
  MatrixXd H = MatrixXd(2, n_x_);

  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;

  H << 1, 0, 0, 0, 0,
       0, 1, 0, 0, 0;

  VectorXd z = VectorXd(2);
  z(0) = meas_package.raw_measurements_[0];
  z(1) = meas_package.raw_measurements_[1];

  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R;
  MatrixXd K = P_ * H.transpose() * S.inverse();

  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H) * P_;

  cout << "Laser updated x_: " << x_ << endl;
  cout << "Laser updated P_: " << P_ << endl;

  //calculate NIS
  NIS_laser_ = y.transpose() * S.inverse() * y;
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
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z(0) = meas_package.raw_measurements_[0];
  z(1) = meas_package.raw_measurements_[1];
  z(2) = meas_package.raw_measurements_[2];

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v  = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double r = sqrt(p_x*p_x + p_y*p_y);
    double phi = atan2(p_y, p_x);
    double r_dot = (p_x*cos(yaw)*v + p_y*sin(yaw)*v) / r;

    Zsig(0, i) = r;  
    Zsig(1, i) = phi;
    Zsig(2, i) = r_dot;
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) 
    z_pred = z_pred + weights_(i) * Zsig.col(i);

    //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
    
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add mearsurement noise
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0., 0.,
       0., std_radphi_*std_radphi_, 0.,
       0., 0., std_radrd_*std_radrd_;
  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
      while(x_diff(3) > M_PI) x_diff(3) -= 2 * M_PI;
      while(x_diff(3) < -M_PI) x_diff(3) += 2 * M_PI;
      
      while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
      while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
      
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  while(z_diff(1) > M_PI) z_diff(1) -= 2 * M_PI;
  while(z_diff(1) < -M_PI) z_diff(1) += 2 * M_PI;
  x_ = x_ + K * z_diff;

  P_ = P_ - K * S * K.transpose();

  cout << "Radar updated x_: " << x_ << endl;
  cout << "Radar updated P_: " << P_ << endl;

  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}
