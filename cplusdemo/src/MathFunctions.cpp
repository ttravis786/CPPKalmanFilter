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

// Custom hash function for std::array<double, 2>
struct ArrayHash {
    size_t operator()(const std::array<double, 2>& arr) const {
        // Use a simple hash combining method; you can customize it based on your needs
        return std::hash<double>()(arr[0]) ^ (std::hash<double>()(arr[1]) << 1);
    }
};

int countCommonElements(std::vector<std::array<double, 2>> *arr1,
                        std::vector<std::array<double, 2>> *arr2) {
    std::unordered_set<std::array<double, 2>, ArrayHash> set1(arr1->begin(), arr1->end());
    int commonCount = 0;
    for (size_t i = 0; i < arr2->size(); ++i) {
        if (set1.count((*arr2)[i]) > 0) {
            // Element from the second array found in the first array
            ++commonCount;
        }
    }
    return commonCount;
}