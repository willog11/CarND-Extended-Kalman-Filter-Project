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

	//measurement matrix - laser
	H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;

	//measurement matrix - radar - Set to 0's for now
	Hj_ << 0, 0, 0, 0,
		0, 0, 0, 0,
		0, 0, 0, 0;

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;

	//initial transition matrix with dt=0
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;

	//inital measurement matrix
	ekf_.H_ = MatrixXd(4, 4);
	ekf_.H_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1;
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
	// first measurement
	cout << "FusionEKF: " << endl;
	ekf_.x_ = VectorXd(4);
	ekf_.x_ << 1, 1, 1, 1;

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		/**
		Convert radar from polar to cartesian coordinates and initialize state.
		*/
		cout << "FusionEKF: Initializing  radar" << endl;
		float rho = measurement_pack.raw_measurements_[0]; // Range - radial distance
		float phi = measurement_pack.raw_measurements_[1]; // Bearing - angel betwee p and x
		float rhodot = measurement_pack.raw_measurements_[2]; // Radial velocity - rate of change of np

		cout << "FusionEKF: Radar measurements read" << endl;
		
		// Polar -> cartesian: x = r * cos(angle), y = r * sin(angle)
		float p1 = rho * cos(phi);
		float p2 = rho * sin(phi);
		float v1 = rhodot * cos(phi);
		float v2 = rhodot * sin(phi);

		cout << "FusionEKF: Updating ekf_x" << endl;
		ekf_.x_ << p1, p2, v1, v2;

		ekf_.R_ << R_radar_;
		ekf_.H_ << Hj_;
		cout << "FusionEKF: Initialized complete" << endl;
	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		/**
		Initialize state.
		*/
		cout << "FusionEKF: Initializing  laser" << endl;
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;

		ekf_.R_ << R_laser_;
		ekf_.H_ << H_laser_;
		cout << "FusionEKF: Initialized complete" << endl;
	}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		cout << "FusionEKF: Initialized" << endl;
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

	//set the acceleration noise components
	float noise_ax = 9;
	float noise_ay = 9;

	//compute the time elapsed between the current and previous measurements
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;

	float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;

	//the initial transition matrix F_
	ekf_.F_ << 1, 0, dt, 0,
				0, 1, 0, dt,
				0, 0, 1, 0,
				0, 0, 0, 1;

	//set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
			0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
			dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
			0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

	cout << "EKF: Prediction state" << endl;
	ekf_.Predict();

	/*****************************************************************************
	*  Update
	****************************************************************************/

	/**
	TODO:
		* Use the sensor type to perform the update step.
		* Update the state and covariance matrices.
	*/

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		Tools tool;
		ekf_.H_ << tool.CalculateJacobian(measurement_pack.raw_measurements_);
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	} else {
		// Laser updates
		ekf_.Update(measurement_pack.raw_measurements_);
	}

	// print the output
	cout << "x_ = " << ekf_.x_ << endl;
	cout << "P_ = " << ekf_.P_ << endl;
}
