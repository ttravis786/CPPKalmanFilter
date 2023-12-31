//
// Created by tomtr on 02/12/2023.
//

#include "MathFunctions.h"


int Sum(int a, int b)
{
    return a + b;
}

double linear_func(std::vector<double> array, std::vector<double> params)
{
    int n = params.size();

    double total_sum =  0;
    for (double x : array) {
        for (int i = 0; i < n; i++) {
            total_sum += params[i] * pow(x, (n - i - 1));
        }
    }
    return total_sum;
}

void linear_regression(Eigen::VectorXd x, Eigen::VectorXd y, int order, double error,
                       Eigen::VectorXd * B, Eigen::MatrixXd * BCovMat, double * var)
{
    int num_points = x.size();
    auto W = error * Eigen::MatrixXd::Identity(num_points,num_points);
    Eigen::MatrixXd X(num_points, order + 1);
    for (int i = 0; i <=order; i++) {
        X.col(i) = x.array().pow(order - i);
    }
    Eigen::MatrixXd X_T = X.transpose();
    // compute coeffiencients, B
    *B = (X_T * W * X).inverse() * X_T * W * y;
    Eigen::VectorXd res = y - X * *B;
    *var = (res.array().pow(2) * (1 / (num_points - order - 1))).sum();
    // compute B_Cov
    *BCovMat = *var * (X_T * X).inverse();
}