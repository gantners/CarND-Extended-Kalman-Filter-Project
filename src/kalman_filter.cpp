#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  CommonUpdate(z - H_ * x_, H_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //map from cartesian to polar coordinates
  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double phi = atan(x_(1) / x_(0));
  double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
  VectorXd hj = VectorXd(3);
  hj << rho, phi, rho_dot;
  CommonUpdate(z - hj, hj);
}

void KalmanFilter::CommonUpdate(const VectorXd &y, const MatrixXd &H){
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H) * P_;
}
