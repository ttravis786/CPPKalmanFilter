//
// Created by tomtr on 02/12/2023.
//
# include <cmath>
# include <vector>
# include <Eigen/Dense>
# include <iostream>
#include <unordered_set>
# include <array>


using Eigen::MatrixXd;

int Sum(int a, int b);
double linear_func(std::vector<double> array, std::vector<double> params);
void linear_regression(Eigen::VectorXd x, Eigen::VectorXd y, int order,
                       double error, Eigen::VectorXd * B, Eigen::MatrixXd * BCovMat, double * var);

struct ArrayHash;

int countCommonElements(std::vector<std::array<double, 2>> *arr1,
                        std::vector<std::array<double, 2>> *arr2);