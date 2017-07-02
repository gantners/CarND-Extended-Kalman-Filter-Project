#include "FusionEKF.h"
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
    H_radar = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0.0,
            0.0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0.0, 0.0,
            0.0, 0.0009, 0.0,
            0.0, 0.0, 0.09;

    H_laser_ << 1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0;
    /**
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
        Init(measurement_pack);
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
       * Update the state transition matrix F according to the new elapsed time.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    if (!Update(measurement_pack))
        return;

    //3. Call the Kalman Filter predict() function
    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     * Update the state and covariance matrices.
     */
    //4. Call the Kalman Filter update() function
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        H_radar = Tools::CalculateJacobian(ekf_.x_);
        ekf_.H_ = H_radar;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    //cout << "x_: " << ekf_.x_ << endl;
    //cout << "P_: " << ekf_.P_ << endl;
}

void FusionEKF::Init(const MeasurementPackage &measurement_pack) {
    //create a 4D state vector, we don't know yet the values of the x state
    VectorXd x = VectorXd(4);
    x << 1, 1, 1, 1; // Important for RSME

    //Init state covariance matrix P
    MatrixXd P = MatrixXd(4, 4);
    P << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    MatrixXd Q(4, 4);


    // first measurement

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
         ro, phi, ro_dot
        */
        double rho = measurement_pack.raw_measurements_[0]; // range
        double phi = measurement_pack.raw_measurements_[1]; // bearing

        // Coordinates conversion from polar to cartesian
        double px = rho * cos(phi);
        double py = rho * sin(phi);
        /**
         * Although radar gives velocity data in the form of the range rate ​ρ,
         * a radar measurement does not contain enough information to determine
         * the state variable velocities vx and vy.
         * You can, however, use the radar measurements ρ and ϕ to initialize
         * the state variable locations p​x and p​y.
         */
        double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
        double vx = rho_dot * cos(phi);
        double vy = rho_dot * sin(phi);
        x << px, py, vx, vy;

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        //set the state with the initial location and zero velocity
        x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    } else {
        throw "Unknown sensor type";
    }
    /*
    float threshold = 0.001;
    if (fabs(x(0)) < threshold && fabs(x(1)) < threshold) {
        cout << "px,py below threshold" << endl;
        x(0) = threshold;
        x(1) = threshold;
    }
    */
    //the initial transition matrix F_
    MatrixXd F = MatrixXd(4, 4);
    F << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    ekf_.Init(x, P, F, Q);

    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
}

bool FusionEKF::Update(const MeasurementPackage &measurement_pack) {
    //compute the time elapsed between the current and previous measurements
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (dt <= 0.005) {
        cout << "Discarding second measure within short timeframe" << endl;
        return false;
    }

    //1. Set the process covariance matrix Q
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;

    //2. Modify the F matrix so that the time is integrated
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
            0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
            dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
            0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;
    return true;
}
