//
// Created by tomtr on 02/12/2023.
//
#include "../src/CombKalFilter.h"
#include "../src/CSVFunctions.h"
#include "../src/MathFunctions.h"

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
//#include <single_include/nlohmann/json.hpp>
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
//using json = nlohmann::json;

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
    std::vector<std::pair<std::string, std::vector<double>>> table = readCSV("fourTracks.csv");
    std::cout << table[0].first;
}

void writeBranchToCSV(std::vector<double> &x, std::vector<double> &y, std::string &name){
    std::vector<std::pair<std::string, std::vector<double>>> table = {{"x", x}, {"y", y}};
    writeCSV(name, table);
}

void saveBranches(std::vector<CombKalFilterBranch> &branches, std::string &name){
    for (int i=0; i<branches.size(); ++i){
        std::vector<double> x;
        std::vector<double> y;
        std::string table_name;
        branches[i].fullPoints(&x, &y);
        table_name = name + "_" + std::to_string(i) + ".csv";
        writeBranchToCSV(x, y, table_name);
    }
}

void assignTrackID(std::vector<CombKalFilterBranch> &branches, int &event,
                   std::vector<int> *trackIDVec , std::vector<float> *momData, std::vector<int> *particleVec){
    double sf = 10;
    std::ostringstream oss;
    oss << "../simulation_data/tenThousandTracks/event" << event << "TrueHits.csv";
    std::string hitsPath = oss.str();
    std::ostringstream ossb;
    ossb << "../simulation_data/tenThousandTracks/event" << event << "TrueStats.csv";
    std::string statsPath = ossb.str();
    std::vector<std::pair<std::string, std::vector<double>>> hitsTable = readCSV(hitsPath);
    std::vector<std::pair<std::string, std::vector<double>>> statsTable = readCSV(statsPath);
    std::vector<double> xTruth = hitsTable[1].second;
    std::vector<double> yTruth = hitsTable[2].second;
    std::vector<double> IDTruth = hitsTable[0].second;
    for (int i=0; i < branches.size(); ++i){
        std::unordered_map<int, int> IDFrequencyMap;
        CombKalFilterBranch branch = branches[i];
        std::vector<double> x;
        std::vector<double> y;
        std::string table_name;
        branches[i].fullPoints(&x, &y);
        for (int j=0; j <x.size(); ++j){
            for (int k=0; k<xTruth.size(); ++k){
                if ((std::fabs(x[j] - xTruth[k]*sf) < 1.5)
                    and (std::fabs(y[j] - yTruth[k]*sf) < 1.5)){
                    IDFrequencyMap[IDTruth[k]]++;
                }
            }
        }
        int branchID = -1;
        int maxFre = 0;
        for (const auto& entry : IDFrequencyMap) {
            if (entry.second > maxFre) {
                branchID = entry.first;
                maxFre = entry.second;
            }
        }
        trackIDVec->emplace_back(branchID);
    }
    for (int i=0; i<trackIDVec->size(); ++i){
        float mom = -1;
        int particle = -1;
        for (int j=0; j<statsTable[0].second.size(); ++j){
            if ((*trackIDVec)[i] == statsTable[0].second[j]) {
                mom = std::sqrt(pow(statsTable[3].second[j],2)
                                      + pow(statsTable[4].second[j],2)
                                      + pow(statsTable[5].second[j],2));
                particle = statsTable[1].second[j];
                break;
            }
        }
        momData->emplace_back(mom);
        particleVec->emplace_back(particle);
    }
    int a =1;
}

void saveEvent(std::vector<CombKalFilterBranch> &branches, std::vector<double> &metaData, std::string &name, int &event){
    std::vector<std::pair<std::string, std::vector<double>>> overallTable;
    std::vector<std::pair<std::string, std::vector<double>>> momentumTable;
    std::string table_name = name + ".csv";
    std::vector<int> trackIDVec;
    std::vector<float> momData;
    std::vector<int> particleVec;
    assignTrackID(branches, event, &trackIDVec, &momData, &particleVec);
    for (int i=0; i<branches.size(); ++i) {
        std::vector<double> x;
        std::vector<double> y;
        branches[i].fullPoints(&x, &y);
        std::vector<double> trackID(x.size(), i);
        std::vector<double> eventID(trackID.size(), event);
        std::vector<double> a = std::vector<double>(1, branches[i].X[0,0]);
        std::vector<double> b = std::vector<double>(1, branches[i].X[1,1]);
        std::vector<double> c = std::vector<double>(1, branches[i].X[2,2]);
        std::vector<double> var_a = std::vector<double>(1, branches[i].P(0,0));
        std::vector<double> var_b = std::vector<double>(1, branches[i].P(1,1));
        std::vector<double> var_c = std::vector<double>(1, branches[i].P(2,2));
        std::vector<double> e = std::vector<double>(1, event);
        std::vector<double> t = std::vector<double>(1, i);
        std::vector<double> nPoints = std::vector<double>(1, x.size());
        std::vector<double> trueID = std::vector<double>(1, trackIDVec[i]);
        std::vector<double> mom = std::vector<double>(1, momData[i]);
        std::vector<double> particle = std::vector<double>(1, particleVec[i]);
        if (i == 0) {
            overallTable = {
                    {"event", eventID},
                    {"track", trackID},
                    {"x",     x},
                    {"y",     y}};
            momentumTable = {
                    {"event", e},
                    {"track", t},
                    {"a", a},
                    {"b", b},
                    {"c", c},
                    {"aVar", var_a},
                    {"bVar", var_b},
                    {"cVar", var_c},
                    {"TrueMom", mom},
                    {"TrueID", trueID},
                    {"TrueParticle", particle},
                    {"numPoints", nPoints}};
        }
        else{
            overallTable[0].second.insert(overallTable[0].second.end(), eventID.begin(), eventID.end());
            overallTable[1].second.insert(overallTable[1].second.end(), trackID.begin(), trackID.end());
            overallTable[2].second.insert(overallTable[2].second.end(), x.begin(), x.end());
            overallTable[3].second.insert(overallTable[3].second.end(), y.begin(), y.end());
            momentumTable[0].second.insert(momentumTable[0].second.end(), e.begin(), e.end());
            momentumTable[1].second.insert(momentumTable[1].second.end(), t.begin(), t.end());
            momentumTable[2].second.insert(momentumTable[2].second.end(), a.begin(), a.end());
            momentumTable[3].second.insert(momentumTable[3].second.end(), b.begin(), b.end());
            momentumTable[4].second.insert(momentumTable[4].second.end(), c.begin(), c.end());
            momentumTable[5].second.insert(momentumTable[5].second.end(), var_a.begin(), var_a.end());
            momentumTable[6].second.insert(momentumTable[6].second.end(), var_b.begin(), var_b.end());
            momentumTable[7].second.insert(momentumTable[7].second.end(), var_c.begin(), var_c.end());
            momentumTable[8].second.insert(momentumTable[8].second.end(), mom.begin(), mom.end());
            momentumTable[9].second.insert(momentumTable[9].second.end(), trueID.begin(), trueID.end());
            momentumTable[10].second.insert(momentumTable[10].second.end(), particle.begin(), particle.end());
            momentumTable[11].second.insert(momentumTable[11].second.end(), nPoints.begin(), nPoints.end());
        }
    }
    if (overallTable.empty()){
        return;
    }
    writeCSV(table_name, overallTable);
    writeCSV(name + "MomentumData.csv", momentumTable);
    std::string metaDataName = name + "MetaData.csv";
    std::vector<std::pair<std::string, std::vector<double>>> metaTable = {
            {"event", std::vector<double>(1, event)},
            {"time", std::vector<double>(1, metaData[0])},
            {"predictedTrackNum", std::vector<double>(1, metaData[1])}};
    writeCSV(metaDataName, metaTable);
}

//void benchMarkModel(){
//    std::ifstream f("simulation_data/ten_track.json");
//    json data = json::parse(f);
//    int b = 5;
//}


void runEvents(){

    bool fixedView = false;
    std::vector<std::array<std::array<int, 2>, 2>> startRange;
    int size = 4;
    double factor = 0.7;
    startRange.reserve(size);
    //Assign values to elements
    for (int i = 0; i < size; ++i) {
        // never start near the end hence factor
        std::array<int, 2> xRange = {static_cast<int>(28 - factor * i * 56 / size) -3,
                                     static_cast<int>(28 - factor * i * 56 / size)};
        std::array<int, 2> yRange = {-25, 25};
        startRange.push_back({xRange, yRange});
    }
//i = 50 bugs, i = 374, i = 408, i = 874
    for (int i=1; i < 3; ++i){
        if (i == 50 or i == 374 or i == 408 or i == 874 or i==1042 or i==1260 or i==1392 or i==1558 or i==1703 or i==1889 or i==1926){
            continue;
        }
        int event = i;
        CombKalFilter p3 = CombKalFilter(
                4, 0.5, 8,
                2, 0.5, startRange,
                0.3, 10, fixedView,
                false, false);
        std::ostringstream oss;
        oss << "../simulation_data/tenThousandTracks/event" << event << ".csv";
        std::string path = oss.str();
        std::vector<std::pair<std::string, std::vector<double>>> table = readCSV(path);
        std::vector<double> x = table[0].second;
        std::vector<double> y = table[1].second;
        Eigen::Map<Eigen::VectorXd> xV(x.data(), x.size());
        Eigen::Map<Eigen::VectorXd> yV(y.data(), y.size());
        p3.addEvent(xV, yV);
        std::pair<std::vector<CombKalFilterBranch>, std::vector<double>> data;
        std::vector<CombKalFilterBranch> branches;
        printf("Event %d started", i);
        data = p3.find_tracks();
        branches = data.first;
        std::vector<double> metaData = data.second;
        std::string branchesName = "tenThousandTrack";
        saveEvent(branches, metaData, branchesName, i);
        printf("Event %d saved", i);
        p3.deleteBranches();
        printf("Event %d deleted", i);
    }

}

int main()
{
    runEvents();

//    //writeCSVTest();
//    //readCSVTest();
//    //benchMarkModel();
//    // only look left!
//    //runEvents();
//    bool fixedView = false;
//    std::array<int, 2> xRange = {25, 28};
//    std::array<int, 2> yRange = {-25, 25};
//    std::array<int, 2> xRangeVertex = {3, 6};
//    std::array<int, 2> xRangeVertex2 = {-6, -3};
//    std::array<int, 2> yRangeVertex2 = {-10,-8};
//    std::vector<std::array<std::array<int, 2>, 2>> startRange;
////    startRange.push_back({xRange, yRange});
////    startRange.push_back({xRangeVertex,
////                          yRange});
////    startRange.push_back({xRangeVertex2,
////                          yRangeVertex2});
//    int size = 4;
//    double factor = 0.7;
//    startRange.reserve(size);
//     //Assign values to elements
//    for (int i = 0; i < size; ++i) {
//        // never start near the end hence factor
//        std::array<int, 2> xRange = {static_cast<int>(28 - factor * i * 56 / size) -3,
//                                     static_cast<int>(28 - factor * i * 56 / size)};
//        std::array<int, 2> yRange = {-25, 25};
//        startRange.push_back({xRange, yRange});
//    }
//// last two are remove and multiple
//    CombKalFilter * p3 = new CombKalFilter(
//            4, 0.5, 7,
//            2, 0.5, startRange,
//            0.3, 10, fixedView,
//            false, false);
//    int event = 0;
//    std::ostringstream oss;
//    oss << "../simulation_data/tenThousandTracks/event" << event << ".csv";
//    std::string path = oss.str();
//    std::vector<std::pair<std::string, std::vector<double>>> table = readCSV(path);
//    //Eigen::VectorXd x {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}};
//    //Eigen::VectorXd y {{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0}};
//    std::vector<double> x = table[0].second;
//    std::vector<double> y = table[1].second;
//    Eigen::Map<Eigen::VectorXd> xV(x.data(), x.size());
//    Eigen::Map<Eigen::VectorXd> yV(y.data(), y.size());
//    p3->addEvent(xV, yV);
//    std::pair<std::vector<CombKalFilterBranch>, std::vector<double>> data;
//    std::vector<CombKalFilterBranch> branches;
//    data = p3->find_tracks();
//    branches = data.first;
//    std::vector<double> metaData = data.second;
////stdout (data.second[1]);
//    std::string branchesName = "sullyEGBestChildTrack";
//    saveEvent(branches, metaData, branchesName, event);
//    p3->deleteBranches();





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
