#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  noise_ax = 9;
  noise_ay = 9;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    Init(measurement_pack);
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  //1. Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //2. Set the process covariance matrix Q
  double t2 = pow(dt, 2);
  double t3 = pow(dt, 3) / 2;
  double t4 = pow(dt, 4) / 4;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << t4 * noise_ax, 0, t3 * noise_ax, 0,
          0, t4 * noise_ay, 0, t3 * noise_ay,
          t3 * noise_ax, 0, t2 * noise_ax, 0,
          0, t3 * noise_ay, 0, t2 * noise_ay;

  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  //3. Call the Kalman Filter predict() function
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  //4. Call the Kalman Filter update() function
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    //Modify the F matrix so that the time is integrated
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

void FusionEKF::Init(const MeasurementPackage &measurement_pack) {
  //create a 4D state vector, we don't know yet the values of the x state
  VectorXd x = VectorXd(4);

  //Init state covariance matrix P
  MatrixXd P = MatrixXd(4, 4);
  P << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  MatrixXd Q = MatrixXd(4, 4);

  //the initial transition matrix F_
  MatrixXd F = MatrixXd(4, 4);
  F << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;

  MatrixXd R;
  MatrixXd H;

  // first measurement
  cout << "EKF: " << endl;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
     ro, phi, ro_dot
    */
    cout << "Initializing radar sensor." << endl;

    double rho = measurement_pack.raw_measurements_[0]; // range
    double phi = measurement_pack.raw_measurements_[1]; // bearing
    double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho

    // Coordinates conversion from polar to cartesian
    double px = rho * cos(phi);
    double py = rho * sin(phi);
    double vx = rho_dot * cos(phi);
    double vy = rho_dot * sin(phi);
    x << px, py, vx, vy;
    R = R_radar_;

  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    /**
    Initialize state.
     x, y
    */
    cout << "Initializing laser sensor." << endl;

    //set the state with the initial location and zero velocity
    x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

    H = H_laser_;
    R = R_laser_;
  } else {
    cout << "Unknown sensor type" << endl;
  }

  if (fabs(x(0)) < 0.0001 && fabs(x(1)) < 0.0001) {
    x(0) = 0.0001;
    x(1) = 0.0001;
  }

  //timestamp update
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.Init(x, P, F, H, R, Q);

  is_initialized_ = true;
}
