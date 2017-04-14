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
    Eigen::VectorXd sum(estimations[0].size());
    sum.fill(0);
    for(int i=0; i < estimations.size(); i++){
        Eigen::VectorXd err = (estimations[i] - ground_truth[i]).array().pow(2);
        sum += err;
    }
    return (sum / estimations.size()).array().sqrt();
}
