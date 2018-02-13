#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* sigma points in measurement space matrix
  MatrixXd Zsig_;

  ///* mean predicted measurement vector
  VectorXd z_pred_;

  ///* measurement covariance matrix
  MatrixXd S_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* NIS values for each processed radar measurement
  std::vector<double> NIS_Radar_;

  ///* NIS values for each processed laser measurement
  std::vector<double> NIS_Laser_;

  ///* Measurement noise covariance matrix for Lidar
  MatrixXd R_Laser_;

  ///* Measurement noise covariance matrix for Radar
  MatrixXd R_Radar_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   *  Compute NIS
   *  @param:
   *    - z, measurement vector at time k+1
   *    - sensor_type, sensor type (RADAR or LASER)
   */
   void ComputeNIS(VectorXd z, MeasurementPackage::SensorType sensor_type);

  /**
   *  Angle normalization to [-Pi, Pi]
   *  @param ang The angle to normalize
   */
  void NormalizeAngle(double *ang);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * GenerateAugmentedSigmaPoints Generates the sigma points for prediction
   * @param None
   */
  void GenerateAugmentedSigmaPoints();

  /**
   * PredictSigmaPoints Predicts sigma points
   * @param delta_t Time between k and k+1 in s
   */
  void PredictSigmaPoints(double delta_t);

  /**
   * PredictStateMeanAndCovariance Predicts the state mean and covariance
   * @param None
   */
  void PredictStateMeanAndCovariance();

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Transform RADAR sigma points into measurement space, calculate the mean
   * predicted measurement vector, and calculate the innovation covariance matrix
   * @param n_z_ measurement dimension
   */
  void PredictRadarMeasurement(int n_z_);

  /**
   * Transform LIDAR sigma points into measurement space, calculate the mean
   * predicted measurement vector, and calculate the innovation covariance matrix
   * @param n_z_ measurement dimension
   */
  void PredictLidarMeasurement(int n_z_);

  /**
   * Update state mean and covariance matrix and compute NIS
   * @param:
   *  - n_z_, measurement dimension
   *  - z, measurement vector
   *  - sensor_type, sensor type (RADAR or LASER)
   */
  void UpdateUKF(int n_z_, VectorXd z, MeasurementPackage::SensorType sensor_type);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
