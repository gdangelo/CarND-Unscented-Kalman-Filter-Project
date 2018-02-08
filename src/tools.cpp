#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Initialize rmse vector
  VectorXd rmse = VectorXd(4);
  rmse.fill(0.0);

  // check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    std::cout << "CalculateRMSE() - ERROR - Invalid inputs for calcaluting RMSE." << std::endl;
    return rmse;
  }

  // Accumulate squared residuals
  for(int i = 0; i < estimations.size(); i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array(); //coefficient-wise multiplication
    rmse += residual;
  }

  // Compute the mean
  rmse = rmse/estimations.size();

  // Compute the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}
