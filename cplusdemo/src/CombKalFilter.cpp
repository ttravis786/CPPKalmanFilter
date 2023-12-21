//
// Created by tomtr on 02/12/2023.
//

#include "CombKalFilter.h"


CombKalFilter::CombKalFilter(
        int maxRange, double binWidth, int minPoints,
        int fitOrder, double RK,
        std::array<std::array<int, 2>, 2> startRange,
        double maxRes, int maxIteration):
        maxRange(maxRange),
        binWidth(binWidth),
        minPoints(minPoints),
        fitOrder(fitOrder),
        RK(RK),
        startRange(startRange),
        maxRes(maxRes),
        maxIteration(maxIteration)
{
    onStart();
}

CombKalFilter::CombKalFilter(){
    onStart();
}

std::vector<std::array<int, 2>> CombKalFilter::generateNextCell1D(int depth, int const_val, bool inclusive,
                                                                  bool doubleLength, bool append_left, bool doubleFixed){
    std::vector<std::array<int, 2>> cells = {};
    int end = depth;
    int start;
    if (doubleLength){
        start = -depth;
        if (!inclusive){
            start +=1;
        }
    }
    else{
        start = 0;
    }
    if (inclusive){
        end += 1;
    }
    for (int i=start; i < end; ++i){
        if (append_left){
            cells.emplace_back(std::array<int, 2> {i, const_val});
            if (doubleFixed){
                cells.emplace_back(std::array<int, 2> {i, -const_val});
            }
        }
        else {
            cells.emplace_back(std::array<int, 2> {const_val, i});
            if (doubleFixed){
                cells.emplace_back(std::array<int, 2> {-const_val, i});
            }
        }
    }
    return cells;
}

void CombKalFilter::generateNextCellVec(){
    for (int i=1; i <= maxRange; ++i){

        // quarter top right
        std::vector<std::array<int, 2>> xInc = generateNextCell1D(i, i, true,
                                                                  false, false, false);
        std::vector<std::array<int, 2>> yInc = generateNextCell1D(i, i, false, false,
                                                                  true, false);
        xInc.insert(xInc.end(), yInc.begin(), yInc.end());
        nextCellVec[0].push_back(xInc);
        //printf("My Integer: %d\n", i);
        // half up
        // first one is half up, second is half right
        for (int j=1; j<3; ++j){
            bool append_left_setting;
            if (j == 1){
                append_left_setting = true;
            }
            else{
                append_left_setting = false;
            }
            std::vector<std::array<int, 2>> yhuIncL = generateNextCell1D(i, i, true, false, append_left_setting, true);
            std::vector<std::array<int, 2>> xhuIncD = generateNextCell1D(i, i, false, true, not append_left_setting, false);
            yhuIncL.insert(yhuIncL.end(), xhuIncD.begin(), xhuIncD.end());
            nextCellVec[j].push_back(yhuIncL);
        }
        std::vector<std::array<int, 2>> yfullInc = generateNextCell1D(i, i, true, true, false, true);
        std::vector<std::array<int, 2>> xfullInc = generateNextCell1D(i, i, false, true, true, true);
        yfullInc.insert(yfullInc.end(), xfullInc.begin(), xfullInc.end());
        nextCellVec[3].push_back(yfullInc);
    }
}

void CombKalFilter::onStart(){
    F = Eigen::MatrixXd::Identity(fitOrder + 1, fitOrder + 1);
    H = Eigen::MatrixXd::Zero(fitOrder + 1, fitOrder + 1);
    generateNextCellVec();
}


std::array<std::array<std::vector<std::array<double, 2>>, 112>,  112> CombKalFilter::getPointsArray()
{
    return pointsArray;
}

std::array<std::vector<std::vector<std::array<int, 2>>>, 4> CombKalFilter::getNextCellVec() const {
    return nextCellVec;
};


void CombKalFilter::getCellDeltas(int currBin, int prevBin, int *mag){
    *mag = currBin - prevBin;
    if (*mag !=0) {
        if (*mag > 0) {
            *mag = 1;
        } else {
            *mag = -1;
        }
    }
}


std::vector<std::array<int, 2>> CombKalFilter::getNextCells(int i, int magX, int magY){
    if ((magX !=0) & (magY !=0)){
        printf("a");
        return nextCellVec[0][i];
    }
    // u/d
    else if (magX !=0){
        printf("b");
        return nextCellVec[1][i];
    }
    // r/l
    else if (magY !=0){
        printf("c");
        return nextCellVec[2][i];
    }
    // full
    else{
        printf("d");
        return nextCellVec[3][i];
    }

}

Eigen::MatrixXd CombKalFilter::getH(double x)
{
    // may want to change H at some point
    for (int i = 0; i <=fitOrder; i++) {
        H(0,i) = pow(x, fitOrder - i);
    }
    return H;
}
// 56 added to remove negatives
int CombKalFilter::binSelector(double x) const {
    return std::floor(x / binWidth) + 56;
}

void CombKalFilter::binFinder(double x, double y, int *binX, int *binY){
    * binX = binSelector(x);
    * binY = binSelector(y);
}

void CombKalFilter::addEvent(Eigen::VectorXd xVec, Eigen::VectorXd yVec){
//    for (int i = 0; i <pointsArray.size(); ++i) {
//        for (int j = 0; j <pointsArray.size(); ++j) {
//            pointsArray[i][j] = std::vector<std::array<double, 2>>;
//        }
//    }
    for(int i=0; i<xVec.size(); i++){
        double x = xVec(i);
        double y = yVec(i);
        int binX;
        int binY;
        binFinder(x, y, &binX, &binY);
        std::array<double, 2> point = {x, y};
        pointsArray[binX][binY].emplace_back(point);
    }
}

bool CombKalFilter::isin_coords(std::array<double,2> coord, std::vector<std::array<double,2>> coordList){
    for (int i=0; i<coordList.size(); ++i){
        if ((coord[0] == coordList[i][0])
        & (coord[1] == coordList[i][1])) {
            return true;
        }
    }
    return false;
}

bool CombKalFilter::isin_coords(double x, double y, std::vector<double> xV, std::vector<double> yV){
    for (int i=0; i<xV.size(); ++i){
        if ((x == xV[i])
            & (y == yV[1])) {
            return true;
        }
    }
    return false;
}

void CombKalFilter::getNextBin(int mag, int relCell, int currBin, int *nextBin){
    if (mag==0){
        *nextBin = relCell + currBin;
    }
    else{
        *nextBin = mag * relCell + currBin;
    }
}

std::vector<CombKalFilterBranch> CombKalFilter::addPoints(
        std::vector<double> x, std::vector<double> y, int currBinX,
        int currBinY, int prevBinX, int prevBinY,
        int currentLev){
    if (currentLev == fitOrder + 2) {
        return processSeed(x, y, currBinX, currBinY, prevBinX, prevBinY);
    }
    std::vector<CombKalFilterBranch> possible_seeds;
    int magX;
    int magY;
    getCellDeltas(currBinX, prevBinX, &magX);
    getCellDeltas(currBinY, prevBinY, &magY);
    // note here zeroth elemnt is 1 out!
    for (int i=0; i<maxRange; ++i){
        std::vector<std::array<int, 2>> relativeCells = getNextCells(i, magX, magY);
        for (int j=0; j<relativeCells.size(); ++j){
            // get cell
            int nextBinX;
            int nextBinY;
            getNextBin(magX, relativeCells[j][0], currBinX, &nextBinX);
            getNextBin(magY, relativeCells[j][1], currBinY, &nextBinY);
            for (int k=0; k < pointsArray[nextBinX][nextBinY].size(); ++k){
                double x_p = pointsArray[nextBinX][nextBinY][k][0];
                double y_p = pointsArray[nextBinX][nextBinY][k][1];
                printf("n %f ", x_p);
                std::vector<double> newX = x;
                std::vector<double> newY = y;
                newX.push_back(x_p);
                newY.push_back(y_p);
                std::vector<CombKalFilterBranch> newSeeds = addPoints(newX, newY, nextBinX,
                                                      nextBinY, currBinX, currBinY,
                                                      currentLev +1);
                possible_seeds.insert(possible_seeds.end(), newSeeds.begin(), newSeeds.end());
            }
        }
        if (possible_seeds.size() != 0 ){
            printf("  empty  ");
            return possible_seeds;
        }
    }
    return possible_seeds;

}

//void CombKalFilter::passToBranch(CombKalFilterBranch* branch) {
//    // Pass 'this' pointer to another object
//    branch->getBranch(this);
//}

void CombKalFilter::findSeeds(std::array<double, 2> point, int i , int j,
                                std::vector<std::array<double, 2>> * usedPoints,
                                std::vector<CombKalFilterBranch> * possibleTracks){
    * possibleTracks  = addPoints(std::vector<double>{point[0]}, std::vector<double>{point[1]}, i,
            j, i+1, i+1,1);
    printf("yo");
    CombKalFilterBranch track;
    std::vector<double> fullXPoints;
    std::vector<double> fullYPoints;
    // first reserve
    usedPoints ->reserve(possibleTracks->size() * 5);
    for (int k=0; k< possibleTracks->size(); ++k){
        track = (*possibleTracks)[k];
        track.fullPoints(&fullXPoints, &fullYPoints);
        for(int l=0; l<5; ++l){
            std::array<double, 2> point = {fullXPoints[l], fullYPoints[l]};
            usedPoints->push_back(std::array<double,2>
                    {track.xPoints[l], track.yPoints[l]});
        }
    }

}

std::vector<CombKalFilterBranch> CombKalFilter::processSeed(std::vector<double> x, std::vector<double> y, int currBinX,
                                int currBinY, int prevBinX, int prevBinY) {
    double var;
    Eigen::VectorXd B;
    Eigen::MatrixXd BCovMat;
    linear_regression(Eigen::Map<Eigen::VectorXd> (x.data(), x.size()),
                      Eigen::Map<Eigen::VectorXd> (y.data(), y.size()),
                      fitOrder, 0.5, &B, &BCovMat, &var);
    std::shared_ptr<CombKalFilter> ckf_pointer(this);
    CombKalFilterBranch masterBranch = CombKalFilterBranch(
            B, BCovMat, currBinX, currBinY,
            prevBinX, prevBinY,
            x,  y, nullptr, x.size(),
            ckf_pointer, 1);

    masterBranch.propogate();
    // will eventually return all the suitable children not just best
    std::shared_ptr<CombKalFilterBranch> bestChild;
    double measure;
    int childPoints;
    masterBranch.getBestChild(bestChild, &measure, &childPoints);
    int gh = 1;
    return std::vector<CombKalFilterBranch> {};
    //return std::vector<std::shared_ptr<CombKalFilterBranch>> {bestChild};
}



std::vector<CombKalFilterBranch> CombKalFilter::find_tracks() {
    printf("yo");
    std::vector<CombKalFilterBranch> possibleTrackVec;
    std::vector<std::array<double, 2>> trackSeedsUsed;
    int minBinY;
    int minBinX;
    int maxBinY;
    int maxBinX;
    binFinder(    startRange[0][0],     startRange[0][1], &minBinX, &maxBinX);
    binFinder(    startRange[1][0],     startRange[1][1], &minBinY, &maxBinY);
    for (int i=maxBinX; i >= minBinX; --i){
        for (int j=maxBinY; j>=minBinY; --j){
            for (int k=0; k<pointsArray[i][j].size(); ++k){
                if (isin_coords(pointsArray[i][j][k], trackSeedsUsed)){
                    continue;
                }
                std::vector<std::array<double, 2>> usedPoints;
                std::vector<CombKalFilterBranch> possibleTracks;
                findSeeds(pointsArray[i][j][k], i , j, &usedPoints, &possibleTracks);
                possibleTrackVec.insert(possibleTrackVec.end(), possibleTracks.begin(), possibleTracks.end());
                trackSeedsUsed.insert(trackSeedsUsed.end(), usedPoints.begin(), usedPoints.end());
            }
        }
    }
    printf("yo");
    return possibleTrackVec;

}



