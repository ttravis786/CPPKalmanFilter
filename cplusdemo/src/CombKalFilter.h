//
// Created by tomtr on 02/12/2023.
//

#ifndef CPLUSDEMO_COMBKALFILTER_H
#define CPLUSDEMO_COMBKALFILTER_H

#include <string>
#include <tuple>
#include "cmath"
#include <Eigen/Dense>
#include <utility>
#include <map>
#include <vector>
#include <iostream>
#include <array>
#include "CombKalFilterBranch.h"
#include "MathFunctions.h"
#include "memory"


class CombKalFilterBranch;

class CombKalFilter
{
private:
public:
    CombKalFilter(
            int maxRange, double binWidth, int minPoints,
            int fitOrder, double RK,
            std::array<std::array<int, 2>, 2> startRange,
            double maxRes, int maxIteration);

    int maxRange = 5;
    double binWidth = 0.5;
    int minPoints = 8;
    int fitOrder=  2;
    std::array<std::array<int, 2>, 2> startRange;
    double RK=pow(0.5,2);
    // unc on y
    double maxRes = 2.0;
    int maxIteration = 10;
    Eigen::MatrixXd F;
    Eigen::MatrixXd H;
    std::array<std::array<std::vector<std::array<double, 2>>, 112>,  112> pointsArray;
    std::array<std::vector<std::vector<std::array<int, 2>>>, 4> nextCellVec;

    Eigen::MatrixXd getH(double x);
    int binSelector(double x) const;
    void binFinder(double x, double y, int *binX, int *binY);
    void addEvent(Eigen::VectorXd x, Eigen::VectorXd y);

    std::array<std::array<std::vector<std::array<double, 2>>, 112>, 112>  getPointsArray();

    CombKalFilter();

    // ~Person(); destructor
    void onStart();

    void generateNextCellVec();

    std::vector<std::array<int, 2>> generateNextCell1D(int depth, int const_val, bool inclusive,
                                                       bool doubleLength, bool append_left, bool doubleFixed);

    std::array<std::vector<std::vector<std::array<int, 2>>>, 4> getNextCellVec() const;

    std::vector<CombKalFilterBranch> find_tracks();


    void findSeeds(std::array<double, 2>, int i, int j, std::vector<std::array<double, 2>> *usedPoints,
                   std::vector<CombKalFilterBranch> *possibleTracks);

    std::vector<CombKalFilterBranch>
    processSeed(std::vector<double> x, std::vector<double> y, int currBinX, int currBinY, int prevBinX, int prevBinY);

    std::vector<CombKalFilterBranch>
    addPoints(std::vector<double> x, std::vector<double> y, int currBinX, int currBinY, int prevBinX, int prevBinY,
              int currentLev);

    void getCellDeltas(int currBin, int prevBin, int *magX);

    std::vector<std::array<int, 2>> getNextCells(int i, int magX, int magY);

    void getNextBin(int mag, int relCell, int currBin, int *nextBin);

    bool isin_coords(std::array<double, 2> coord, std::vector<std::array<double, 2>> coordList);

    bool isin_coords(double x, double y, std::vector<double> xV, std::vector<double> yV);
};

#endif //CPLUSDEMO_COMBKALFILTER_H