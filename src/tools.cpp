#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
  	cout << "Invalid estimation or ground truth data" << endl;
  	return rmse;
  }

  //accumulate squared residuals
  VectorXd residuals(4);
  residuals << 0, 0, 0, 0;
  for (int i = 0; i < estimations.size(); i++) {
  	VectorXd residual = (estimations[i] - ground_truth[i]).array() * (estimations[i] - ground_truth[i]).array();
  	residuals += residual;
  }

  //calculate the mean
  VectorXd mean = residuals / estimations.size();

  //calculate the squared root
  rmse = mean.array().sqrt();

  //return the result
  return rmse;
}