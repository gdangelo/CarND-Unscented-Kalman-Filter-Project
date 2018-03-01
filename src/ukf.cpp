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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Set state dimension
  n_x_ = 5;

  // Set augmented dimension
  n_aug_ = 7;

  // Set sigma points dimension
  n_sig_ = 2 * n_aug_ + 1;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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

  // Set augmented Sigma points matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sig_);

  // Set Sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  // Set measurement noise covariance matrix for Radar
  R_Radar_ = MatrixXd(3, 3);
  R_Radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;

  // Set measurement noise covariance matrix for Lidar
  R_Laser_ = MatrixXd(2, 2);
  R_Laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  // Set weights
  weights_ = VectorXd(n_sig_);

  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for(int i = 1; i < n_sig_; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
}

UKF::~UKF() {}

/**
 *  Compute NIS
 *  @param:
 *    - z, measurement vector at time k+1
 *    - sensor_type, sensor type (RADAR or LASER)
 */
void UKF::ComputeNIS(VectorXd z, MeasurementPackage::SensorType sensor_type) {
  double epsilon;

  VectorXd z_diff = z - z_pred_;
  epsilon = z_diff.transpose() * S_.inverse() * z_diff;

  if(sensor_type == MeasurementPackage::LASER) {
    NIS_Laser_ = epsilon;
  }
  if(sensor_type == MeasurementPackage::RADAR) {
    NIS_Radar_ = epsilon;
  }
}

/**
 *  Angle normalization to [-Pi, Pi]
 */
void UKF::NormalizeAngle(double *angle) {
  *angle = atan2(sin(*angle), cos(*angle));
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/

  if(!is_initialized_) {

    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      cout << "**********\nINITIALIZE with LASER\n**********" << endl;

      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      // Handle small px, py
      if(fabs(px) < 0.001){
       px = 0.001;
      }
      if(fabs(py) < 0.001){
       py = 0.001;
      }

      // Initialize state
      x_ << px, py, 0, 0, 0;

      // done initializing, no need to predict or update
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;

      PrintStateMeanAndCovariance();
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      cout << "**********\nINITIALIZE with RADAR\n**********" << endl;

      // Convert radar from polar to cartesian coordinates
      double ro = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      double rho_dot = meas_package.raw_measurements_(1);
      double px = ro*cos(theta);
      double py = ro*sin(theta);
      double vx = rho_dot*cos(theta);
      double vy = rho_dot*sin(theta);

      // Initialize state
      x_ << px, py, sqrt(vx*vx + vy*vy), 0, 0;

      // done initializing, no need to predict or update
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;

      PrintStateMeanAndCovariance();
    }

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
 * GenerateAugmentedSigmaPoints Generates the sigma points for prediction
 * @param None
 */
void UKF::GenerateAugmentedSigmaPoints() {
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

  // Set augmented sigma point matrix
  Xsig_aug_.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
}

/**
 * PredictSigmaPoints Predicts sigma points
 * @param delta_t Time between k and k+1 in s
 */
void UKF::PredictSigmaPoints(double delta_t) {
  // Iterate over each sigma points
  for(int i = 0; i < n_sig_; i++) {
    double px = Xsig_aug_(0, i);
    double py = Xsig_aug_(1, i);
    double v = Xsig_aug_(2, i);
    double phi = Xsig_aug_(3, i);
    double phi_dot = Xsig_aug_(4, i);
    double nu_a = Xsig_aug_(5, i);
    double nu_phi = Xsig_aug_(6, i);

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
}

/**
 * PredictStateMeanAndCovariance Predicts the state mean and covariance
 * @param None
 */
void UKF::PredictStateMeanAndCovariance() {
  // Predict state mean
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for(int i = 0; i < n_sig_; i++) {
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
  }

  // Predict state covariance matrix
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  for(int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    // angle normalization
    NormalizeAngle(&(x_diff(3)));

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose();
  }

  // Update with predictions
  x_ = x_pred;
  P_ = P_pred;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  cout << "**********\nPREDICTION\n**********" << endl;

  GenerateAugmentedSigmaPoints();
  PredictSigmaPoints(delta_t);
  PredictStateMeanAndCovariance();
  PrintStateMeanAndCovariance();
}

/**
 * Transform LIDAR sigma points into measurement space, calculate the mean
 * predicted measurement vector, and calculate the innovation covariance matrix
 * @param n_z_ measurement dimension
 */
void UKF::PredictLidarMeasurement(int n_z_) {
  // Set matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, n_sig_);

  // Set mean predicted measurement
  z_pred_ = VectorXd(n_z_);

  // Set easurement covariance matrix
  S_ = MatrixXd(n_z_, n_z_);

  // Transform sigma points into measurement space
  for(int i = 0; i < n_sig_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    Zsig_.col(i) << px, py;
  }

  // Calculate mean predicted measurement
  z_pred_.fill(0.0);
  for(int i = 0; i < n_sig_; i++) {
    z_pred_ += weights_(i) * Zsig_.col(i);
  }

  // Calculate innovation covariance matrix
  S_.fill(0.0);
  for(int i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    S_ += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add noise
  S_ = S_ + R_Laser_;
}

/**
 * Transform RADAR sigma points into measurement space, calculate the mean
 * predicted measurement vector, and calculate the innovation covariance matrix
 * @param n_z_ measurement dimension
 */
void UKF::PredictRadarMeasurement(int n_z_) {
  // Set matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, n_sig_);

  // Set mean predicted measurement
  z_pred_ = VectorXd(n_z_);

  // Set easurement covariance matrix
  S_ = MatrixXd(n_z_, n_z_);

  // Transform sigma points into measurement space
  for(int i = 0; i < n_sig_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    if (fabs(px) > 0.001 || fabs(py) > 0.001) {
      Zsig_(0,i) = sqrt(px*px + py*py);
      Zsig_(1,i) = atan2(py, px);
      Zsig_(2,i) = (px*cos(yaw)*v + py*sin(yaw)*v) / sqrt(px*px + py*py);
    }
    else {
      Zsig_(0,i) = 0.0;
      Zsig_(1,i) = 0.0;
      Zsig_(2,i) = 0.0;
    }
  }

  // Calculate mean predicted measurement
  z_pred_.fill(0.0);
  for(int i = 0; i < n_sig_; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  // Calculate innovation covariance matrix S_
  S_.fill(0.0);
  for(int i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    // angle normalization
    NormalizeAngle(&(z_diff(1)));

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add noise
  S_ = S_ + R_Radar_;
}

/**
 * Update state mean and covariance matrix and compute NIS
 * @param:
 *  - n_z_, measurement dimension
 *  - z, measurement vector
 *  - sensor_type, sensor type (RADAR or LASER)
 */
void UKF::UpdateUKF(int n_z_, VectorXd z, MeasurementPackage::SensorType sensor_type) {
  // Calculate Cross-correlation Matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  for(int i = 0; i < n_sig_; i++) {
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    if(sensor_type == MeasurementPackage::RADAR) {
      // angle normalization
      NormalizeAngle(&(z_diff(1)));
      NormalizeAngle(&(x_diff(3)));
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Calculate Kalman gain
  MatrixXd K = Tc * S_.inverse();

  // residual
  VectorXd z_diff = z - z_pred_;

  if(sensor_type == MeasurementPackage::RADAR) {
    // angle normalization
    NormalizeAngle(&(z_diff(1)));
  }

  // Update state mean and covariance matrix
  x_ = x_ +  K * z_diff;
  P_ = P_ - K * S_ * K.transpose();

  // Compute and store NIS for this step
  ComputeNIS(z, sensor_type);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  cout << "***********\nUPDATE LIDAR\n***********" << endl;

  // Set measurement dimension, lidar can measure px, py
  int n_z_ = 2;

  PredictLidarMeasurement(n_z_);
  UpdateUKF(n_z_, meas_package.raw_measurements_, meas_package.sensor_type_);
  PrintStateMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  cout << "***********\nUPDATE RADAR\n***********" << endl;

  // Set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  PredictRadarMeasurement(n_z_);
  UpdateUKF(n_z_, meas_package.raw_measurements_, meas_package.sensor_type_);
  PrintStateMeanAndCovariance();
}

/**
 * Print state mean vector and state covariance matrix
 * @param None
 */
void UKF::PrintStateMeanAndCovariance() {
  cout << "x_:\n" << x_ << "\n" << endl;
  cout << "P_:\n" << P_ << "\n" << endl;
}
