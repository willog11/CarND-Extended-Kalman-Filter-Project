#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

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
  /**
  TODO:
    * predict the state
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
	cout << "KalmanFilter: Prediction state complete" << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
	cout << "KalmanFilter: Update state complete" << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	cout << "KalmanFilter::UpdateEKF() - started" << endl;
	float px = z[0];
	float py = z[1];
	float vx = z[2];
	float vy = z[3];

	float rho = sqrt((px*px) + (py*py));
	float phi = atan(py / px);

	//check division by zero
	if (fabs(rho) < 0.0001) {
		cout << "KalmanFilter::UpdateEKF() - Error - Division by Zero" << endl;
		px = 0.0001;
		py = 0.0001;
		rho = sqrt((px*px) + (py*py));
	}

	float rhodot = (px*vx + py*vy) / rho;

	while (phi < -M_PI || phi > M_PI)
	{
		phi += 2 * M_PI;
	}
	cout << "KalmanFilter::UpdateEKF() - phi ="<< phi << endl;

	VectorXd hx(3);
	hx << rho, phi,	rhodot;
	cout << "KalmanFilter::UpdateEKF() - hx updated" << endl;

	VectorXd y = z - hx;

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
	cout << "KalmanFilter::UpdateEKF() state complete" << endl;
}
