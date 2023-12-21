//
// Created by tomtr on 02/12/2023.
//
#include "../src/CombKalFilter.h"
#include "../src/CSVFunctions.h"
#include "../src/MathFunctions.h"

#include <iostream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

void writeCSVTest(){
    std::vector<double> vec1(100, 1.0);
    std::vector<double> vec2(100, 2.0);
    std::vector<double> vec3(100, 3.0);

    // Wrap into a vector
    std::vector<std::pair<std::string, std::vector<double>>> table = readCSV("twoTracks.csv");

    // Write the vector to CSV
    writeCSV("three_cols.csv", table);
}

void readCSVTest(){
    std::vector<std::pair<std::string, std::vector<double>>> table = readCSV("twoTracks.csv");
    std::cout << table[0].first;
}

int main()
{
    //writeCSVTest();
    //readCSVTest();
    std::array<int, 2> xRange = {6, 8};
    std::array<int, 2> yRange = {6, 8};
    std::array<std::array<int, 2>, 2> startRange{xRange, yRange};
    CombKalFilter p3(
            4, 0.5, 5,
            2, 0.5, startRange,
            0.5, 10);
    Eigen::VectorXd x {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}};
    Eigen::VectorXd y {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}};
    p3.addEvent(x, y);
    p3.find_tracks();

    // orientation tests
//    std::array<std::array<std::vector<std::array<double, 2>>, 112>, 112> points_array = p3.getPointsArray();
//    std::array<std::vector<std::vector<std::array<int, 2>>>, 4> nextCells = p3.getNextCellVec();
//    for (int i = 0; i <nextCells.size(); ++i){
//        for (int j = 0; j <nextCells[i].size(); ++j){
//            for (int k = 0; k <nextCells[i][j].size(); ++k){
//                printf("orientation %d depth %d Coord: (%d, %d)\n", i, j, nextCells[i][j][k][0], nextCells[i][j][k][1]);
//            }
//        }
//    }
//    for (int i = 0; i <points_array.size(); ++i) {
//        for (int j = 0; j <points_array.size(); ++j) {
//            for (int k = 0; k<points_array[i][j].size(); ++k){
//                std::cout << points_array[i][j][k][0] << " ";
//            }
//        }
//    }
//    std::string name = p1.getName();
//    std::cout << name;
//    std::cout << Sum(1,3);
//    int order = 3;
//    double error = 0.5;
//    linear_regression(x, y, order, error);
//    return 0;
}
