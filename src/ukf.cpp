#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Set state dimension
  n_x_ = 5;

  // Set augmented dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  is_initialized_ = false;

  // Initialize the state covariance matrix as the Identity matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Define spreading parameter
  lambda_ = 3 - n_aug_;

  // Set Sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 *  Angle normalization to [-Pi, Pi]
 */
void UKF::NormalizeAngle(double *ang) {
  *ang = atan2(sin(*ang), cos(*ang));
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  cout << "---\n> Process new measurement (" << meas_package.sensor_type_ << ")" << endl;

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/

  if(!is_initialized_) {
    cout << "> Initialize" << endl;
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
     // Initialize state
     x_.fill(0.0);
     x_(0) = meas_package.raw_measurements_(0); //px
     x_(1) = meas_package.raw_measurements_(1); //py
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
     // Convert radar from polar to cartesian coordinates
     float ro = meas_package.raw_measurements_(0);
     float theta = meas_package.raw_measurements_(1);
     // Initialize state
     x_.fill(0.0);
     x_(0) = ro*cos(theta); //px
     x_(1) = ro*sin(theta); //py
    }

    time_us_ = meas_package.timestamp_;

    cout << "x_:\n" << x_ << endl;
    cout << "P_:\n" << P_ << endl;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // Skip prediction & update if sensor type is ignored
  if((meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) ||
    (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)) {

    /*****************************************************************************
    *  Prediction
    ****************************************************************************/

    // Compute delat_t and update previous measurement time
    double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; // in seconds
    time_us_ = meas_package.timestamp_;

    // Call prediction process
    Prediction(delta_t);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  cout << "> Predict:" << endl;

  /*****************************************************************************
   *  Generate sigma points
   ****************************************************************************/

  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

  // Create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // Create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_)*A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
  }

  /*****************************************************************************
   *  Predict sigma points using unscented transform
   ****************************************************************************/

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double phi = Xsig_aug(3, i);
    double phi_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_phi = Xsig_aug(6, i);

    double px_pred, py_pred;
    if(fabs(phi_dot) > 0.001) {
      px_pred = px + v/phi_dot * (sin(phi + phi_dot*delta_t) - sin(phi));
      py_pred = py + v/phi_dot * (-cos(phi + phi_dot*delta_t) + cos(phi));
    }
    else {
      px_pred = px + v*cos(phi)*delta_t;
      py_pred = py + v*sin(phi)*delta_t;
    }

    double v_pred = v;
    double phi_pred = phi + phi_dot*delta_t;
    double phi_dot_pred = phi_dot;

    // Add noise
    px_pred = px_pred + 0.5*nu_a*delta_t*delta_t * cos(phi);
    py_pred = py_pred + 0.5*nu_a*delta_t*delta_t * sin(phi);
    v_pred = v_pred + nu_a*delta_t;

    phi_pred = phi_pred + 0.5*nu_phi*delta_t*delta_t;
    phi_dot_pred = phi_dot_pred + nu_phi*delta_t;

    Xsig_pred_.col(i) << px_pred, py_pred, v_pred, phi_pred, phi_dot_pred;
  }

  /*****************************************************************************
   *  Predict state mean and state covariance matrix
   ****************************************************************************/

  // Set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for(int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // Predict state mean
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
  }

  // Predict state covariance matrix
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    // angle normalization
    NormalizeAngle(&(x_diff(3)));

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();
  }

  x_ = x_pred;
  P_ = P_pred;

  cout << "x_:\n" << x_ << "\n" << endl;
  cout << "P_:\n" << P_ << "\n" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  cout << "> Update LIDAR" << endl;

  // Set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 2;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  // Transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    Zsig.col(i) << px, py;
  }

  // Calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Calculate innovation covariance matrix S
  S.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z_,n_z_);
  R(0, 0) = std_laspx_*std_laspx_;
  R(1, 1) = std_laspy_*std_laspy_;

  S += R;

  // Calculate Cross-correlation Matrix
  MatrixXd Tc = MatrixXd(n_x_,n_z_);
  Tc.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  }

  // Calculate Kalman gain
  MatrixXd K = MatrixXd(n_x_,n_z_);
  K = Tc * S.inverse();

  // Update state mean and covariance matrix
  x_ = x_ +  K * (meas_package.raw_measurements_ - z_pred);
  cout << "x_:\n" << x_ << "\n" << endl;
  P_ = P_ - K * S * K.transpose();
  cout << "P_:\n" << P_ << "\n" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  cout << "> Update RADAR" << endl;

  // Set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);

  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  // Transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double phi = Xsig_pred_(3, i);

    double ro = sqrt(px*px + py*py);
    double theta = atan2(py, px);
    double ro_dot;
    if(ro == 0) {
      ro_dot = 0;
    }
    else {
      ro_dot = (px*cos(phi)*v + py*sin(phi)*v) / ro;
    }

    Zsig.col(i) << ro, theta, ro_dot;
  }

  // Calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Calculate innovation covariance matrix S
  S.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    NormalizeAngle(&(z_diff(1)));

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z_,n_z_);
  R(0, 0) = std_radr_*std_radr_;
  R(1, 1) = std_radphi_*std_radphi_;
  R(2, 2) = std_radrd_*std_radrd_;

  S = S + R;

  // Calculate Cross-correlation Matrix
  MatrixXd Tc = MatrixXd(n_x_,n_z_);
  Tc.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    NormalizeAngle(&(z_diff(1)));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    NormalizeAngle(&(x_diff(3)));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Calculate Kalman gain
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  // angle normalization
  NormalizeAngle(&(z_diff(1)));

  // Update state mean and covariance matrix
  x_ = x_ +  K * z_diff;
  P_ = P_ - K * S * K.transpose();

  cout << "x_:\n" << x_ << "\n" << endl;
  cout << "P_:\n" << P_ << "\n" << endl;
}
