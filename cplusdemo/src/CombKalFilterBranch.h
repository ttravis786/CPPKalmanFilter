//
// Created by tomtr on 02/12/2023.
//
#ifndef CPLUSDEMO_COMBKALFILTERBRANCH_H
#define CPLUSDEMO_COMBKALFILTERBRANCH_H
# include <cmath>
# include <vector>
# include <Eigen/Dense>
#include "CombKalFilter.h"
#include <memory>

struct nextPointData {
    double x;
    double y;
    int binX;
    int binY;
    Eigen::MatrixXd X;
    Eigen::MatrixXd P;

    nextPointData(double x, double y, int binX, int binY, Eigen::MatrixXd X, Eigen::MatrixXd P);
};

class CombKalFilter;

class CombKalFilterBranch
{
private:
public:


    CombKalFilterBranch(Eigen::VectorXd X, Eigen::MatrixXd P, int currBinX, int currBinY, int prevBinX, int prevBinY,
                        std::vector<double> x, std::vector<double> y, CombKalFilterBranch *master, int totalPoints,
                        std::shared_ptr<CombKalFilter> ckf, int iteration);

//    CombKalFilterBranch(Eigen::VectorXd X, Eigen::MatrixXd P, int currBinX, int currBinY, int prevBinX, int prevBinY,
//                        std::vector<double> x, std::vector<double> y, int totalPoints, CombKalFilter *ckf,
//                        int iteration);

    Eigen::VectorXd X;
    Eigen::MatrixXd P;
    int currBinX;
    int currBinY;
    int prevBinX;
    int prevBinY;
    std::vector<double> xPoints;
    std::vector<double> yPoints;
    CombKalFilterBranch  * master;
    std::vector<CombKalFilterBranch> children;
    int totalPoints;
    std::shared_ptr<CombKalFilter> ckf;
    int iteration;
    std::vector<nextPointData> possibleNext;



    CombKalFilterBranch()=default;
    // ~Person(); destructor
    void propogate();

    void residFunc(double x, double *maxResPoint);

    void goodnessMeasure(double *measure);

    void fullPoints(std::vector<double> *fullXPoints, std::vector<double> *fullYPoints);

    void getBestChild(std::shared_ptr<CombKalFilterBranch> bestchild, double *measure, int *childPoints);

    void updateKF(double xP, double yP, Eigen::VectorXd *nX, Eigen::MatrixXd *nP, double *res);
};


#endif //CPLUSDEMO_COMBKALFILTERBRANCH_H