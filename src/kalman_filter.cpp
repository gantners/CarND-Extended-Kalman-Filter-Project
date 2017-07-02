#include "kalman_filter.h"
#include "tools.h"
#include <cmath>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    //cout << "Predict" << endl;
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft;
    P_ = P_ + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    //cout << "Update" << endl;
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    //cout << "y: " << y(1) << endl;
    CommonUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    //cout << "Update EKF" << endl;
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    //Do not let rho be smaller than 1e-5 avoiding bad values
    double rho = std::max(1e-5, sqrt(px * px + py * py));
    double phi = atan2(py, px);
    double rho_dot = (px * vx + py * vy) / (rho);

    VectorXd z_pred = VectorXd(3);
    z_pred << rho, phi, rho_dot;
    VectorXd y = z - z_pred;
    //Correcting angles as suggested on tips and tricks
    //cout << "y EKF: " << y(1) << endl;
    if (y(1) < -M_PI || y(1) > M_PI) {
        //cout << "y(1) out of range: " << y(1) << endl;
        y(1) = fmod(y(1),M_PI);
        //cout << "y(1) corrected: " << y(1) << endl;
    }
    CommonUpdate(y);
}

void KalmanFilter::CommonUpdate(const VectorXd &y) {
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht;
    S = S + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
