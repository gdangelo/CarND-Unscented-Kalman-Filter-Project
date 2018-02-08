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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  // Initialize the state covariance matrix as the Identity matrix
  P_ = MatrixXd::Identity(5, 5);

  // Set state dimension
  int n_x_ = 5;

  // Set augmented dimension
  int n_aug_ = 7;

  // Define spreading parameter
  double lambda_ = 3 - n_aug_;

  // Set Sigma points matrix
  Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if(!is_initialized_) {

    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state
      x_.fill(0.0);
      x_(0) = meas_package.raw_measurements_(0); //px
      x_(0) = meas_package.raw_measurements_(1); //py
    }

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      float ro = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      // Initialize state
      x_.fill(0.0);
      x_(0) = ro*cos(theta); //px
      x_(0) = ro*sin(theta); //py
    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

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

   // TODO
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /*
   * 1. Generates sigma points
  */

  // Create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  // Create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_aug_-2,n_aug_-2) = std_a_*std_a_;
  P_aug(n_aug_-1,n_aug_-1) = std_yawdd_*std_yawdd_;

  // Create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  MatrixXd A = MatrixXd(n_aug_, n_aug_);
  A = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_x_)*A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_x_)*A.col(i);
  }

  /*
   * 2. Predict sigma points using unscented transform
  */

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double phi = Xsig_aug(3, i);
    double phi_dot = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_phi = Xsig_aug(6, i);

    double px_pred, py_pred;
    if(phi_dot == 0) {
      px_pred = px + v*cos(phi)*delta_t;
      py_pred = py + v*sin(phi)*delta_t;
    }

    px_pred += 0.5 * delta_t*delta_t * cos(phi) * nu_a;
    py_pred += 0.5 * delta_t*delta_t * sin(phi) * nu_a;

    double v_pred = v + delta_t * nu_a;
    double phi_pred = phi + phi_dot*delta_t + 0.5 * delta_t*delta_t*nu_phi;
    double phi_dot_pred = phi_dot + delta_t*nu_phi;

    Xsig_pred_.col(i) << px_pred, py_pred, v_pred, phi_pred, phi_dot_pred;
  }

  /*
   * 3. Predict state mean and state covariance matrix
  */

  // Set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 1 / (2*(lambda_ + n_aug_));
  }

  // Predict state mean
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  // Predict state covariance matrix
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    P_ += weights_(i) * (Xsig_pred_.col(i) - x_) * (Xsig_pred_.col(i) - x_).transpose();
  }
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
