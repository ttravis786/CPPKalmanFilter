//
// Created by tomtr on 02/12/2023.
//

#include "CombKalFilter.h"


CombKalFilter::CombKalFilter(
        int maxRange, double binWidth, int minPoints,
        int fitOrder, double RK,
        std::vector<std::array<std::array<int, 2>, 2>>  startRange,
        double maxRes, int maxIteration, bool fixedView,
        bool remove, bool multiple):
        maxRange(maxRange),
        binWidth(binWidth),
        minPoints(minPoints),
        fitOrder(fitOrder),
        RK(RK),
        startRange(startRange),
        maxRes(maxRes),
        maxIteration(maxIteration),
        fixedView(fixedView),
        remove(remove),
        multiple(multiple)

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

        return nextCellVec[0][i];
    }
    // u/d
    else if (magX !=0){
        //printf("b");
        return nextCellVec[1][i];
    }
    // r/l
    else if (magY !=0){
        //printf("c");
        return nextCellVec[2][i];
    }
    // full
    else{
        //printf("d");
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
    int num =  std::floor(x / binWidth) + 56;
    if (num >= 111){
        return 111;
    }
    if (num < 0){
        return 0;
    }
    return num;
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
        double x = xVec[i];
        double y = yVec[i];
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
        if ((std::fabs(x - xV[i]) < 0.01)
            & (std::fabs(y - yV[i]) <0.01)) {
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
    if (*nextBin >= 112){
        *nextBin = 111;
    }
    else if (*nextBin < 0){
        *nextBin = 0;
    }
}

void CombKalFilter::addPoints(
        std::vector<double> x, std::vector<double> y, int currBinX,
        int currBinY, int prevBinX, int prevBinY,
        int currentLev, std::vector<childData> *possibleSeedsData,
        std::vector<CombKalFilterBranch> *possibleSeeds){;
    if (currentLev == fitOrder + 2) {
        processSeed(x, y, currBinX, currBinY, prevBinX, prevBinY, possibleSeedsData, possibleSeeds);
        if (possibleSeeds->size() > 0) {
            CombKalFilterBranch branch = (*possibleSeeds)[0];
            int a = 1;
        }
        return;
    }
    int magX;
    int magY;
    if (fixedView){
        magX = -1;
        magY = 0;
    }
    else {
        getCellDeltas(currBinX, prevBinX, &magX);
        getCellDeltas(currBinY, prevBinY, &magY);
    }
    // note here zeroth element is 1 out!
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
                //printf("n %f ", x_p);
                std::vector<double> newX = x;
                std::vector<double> newY = y;
                newX.push_back(x_p);
                newY.push_back(y_p);
                std::vector<childData> newSeedsData;
                std::vector<CombKalFilterBranch> newSeeds;
                addPoints(newX, newY, nextBinX, nextBinY,
                          currBinX, currBinY, currentLev +1, &newSeedsData, &newSeeds);
                possibleSeeds->insert(possibleSeeds->end(), newSeeds.begin(), newSeeds.end());
                possibleSeedsData->insert(possibleSeedsData->end(), newSeedsData.begin(), newSeedsData.end());
            }
        }
        if (possibleSeeds->size() != 0 ){
            return;
        }
    }

}

//void CombKalFilter::passToBranch(CombKalFilterBranch* branch) {
//    // Pass 'this' pointer to another object
//    branch->getBranch(this);
//}

void CombKalFilter::findSeeds(std::array<double, 2> point, int i , int j,
                                std::vector<std::array<double, 2>> * usedPoints,
                                std::vector<childData> * possibleTrackData,
                                std::vector<CombKalFilterBranch> * possibleTracks){
     addPoints(std::vector<double>{point[0]}, std::vector<double>{point[1]}, i,
            j, i+1, i+1,1, possibleTrackData, possibleTracks);
    CombKalFilterBranch track;
    // first reserve
    usedPoints ->reserve(possibleTracks->size() * 5);
    for (int k=0; k< possibleTracks->size(); ++k) {
        for (int l = 0; l < 5; ++l) {
            usedPoints->push_back((*possibleTrackData)[k].points[l]);
        }
    }
}

void CombKalFilter::processSeed(std::vector<double> x, std::vector<double> y, int currBinX,
                           int currBinY, int prevBinX, int prevBinY, std::vector<childData> *bestChildData,
                           std::vector<CombKalFilterBranch> *bestChild) {
    double var;
    Eigen::VectorXd B;
    Eigen::MatrixXd BCovMat;
    Eigen::VectorXd B_2;
    Eigen::MatrixXd BCovMat_2;
    linear_regression(Eigen::Map<Eigen::VectorXd> (x.data(), x.size()),
                      Eigen::Map<Eigen::VectorXd> (y.data(), y.size()),
                      fitOrder, 0.5, &B, &BCovMat, &var);
//    linear_regression(Eigen::Map<Eigen::VectorXd> (x.data(), x.size()),
//                      Eigen::Map<Eigen::VectorXd> (y.data(), y.size()),
//                      fitOrder-1, 0.5, &B_2, &BCovMat_2, &var);
//    B(0,0) = 0.0;
//    for (int i=1; i <=2; ++i){
//        for (int j=1; j <=2; ++j){
//            B(i, j) = B_2(i-1,j-1);
//            BCovMat(i,j) = BCovMat_2(i-1, j-1);
//        }
//    }

    //std::shared_ptr<CombKalFilter> ckf_pointer(this);
    CombKalFilterBranch * masterBranch = new CombKalFilterBranch(
            B, BCovMat, currBinX, currBinY,
            prevBinX, prevBinY,
            x,  y, nullptr, x.size(),
            this, 1);

    (*masterBranch).propogate();
    // will eventually return all the suitable children not just best
    //std::shared_ptr<CombKalFilterBranch> bestChild;
    (*masterBranch).getBestChild(bestChildData);
    masterBranches.emplace_back(masterBranch);
    for (int i=0; i<bestChildData->size(); ++i){
        bestChild->emplace_back(*(*bestChildData)[i].branch);
    }

    //return std::vector<std::shared_ptr<CombKalFilterBranch>> {bestChild};
}

void CombKalFilter::deleteBranches(){
    for (int i=0; i<masterBranches.size(); ++i){
        (*masterBranches[i]).deleteBranches();
        delete masterBranches[i];
    }
}

void CombKalFilter::finalMerge(std::vector<childData> *possibleTrackData,
                               std::vector<CombKalFilterBranch> *possibleTrackVec,
                               std::vector<CombKalFilterBranch> *trackVec){
    std::vector<int> removeIndexes;
    for (int i=0; i <possibleTrackData->size(); i++){
        for (int j=i+1; j <possibleTrackData->size(); j++){
            if (similarity_check((*possibleTrackData)[i].points, (*possibleTrackData)[j].points)){
                removeIndexes.emplace_back(j);
                // try measure instead
                double mes1;
                double mes2;
                (*possibleTrackVec)[i].goodnessMeasure(&mes1);
                (*possibleTrackVec)[j].goodnessMeasure(&mes2);
                //if (mes1 > mes2){
                if ((*possibleTrackData)[i].points.size() < (*possibleTrackData)[j].points.size()) {
                    (*possibleTrackData)[i] = (*possibleTrackData)[j];
                    (*possibleTrackVec)[i] = (*possibleTrackVec)[j];
                }
            }
        }
    }
    std::sort(removeIndexes.begin(), removeIndexes.end());

    // Use std::unique to get unique values
    auto uniqueEnd = std::unique(removeIndexes.begin(), removeIndexes.end());

    // Resize the vector to contain only unique values
    removeIndexes.resize(std::distance(removeIndexes.begin(), uniqueEnd));
    std::vector<int> keepIndexes;
    bool include;
    for (int i=0; i <(*possibleTrackData).size(); ++i){
        include = true;
        for (int j=0; j <removeIndexes.size(); ++j){
            if (removeIndexes[j] == i){
                include = false;
                break;
            }
        }
        // duplicate check
        for (int k=0; k<trackVec->size(); ++k){
            if (((*possibleTrackVec)[i]).X == (*trackVec)[k].X){
                include = false;
                break;
            }
        }
        if (include){
            keepIndexes.push_back(i);
            trackVec->push_back((*possibleTrackVec)[i]);
        }
    }
    // occasionally algo produce duplicates, remove them

    int b = 1;
}

std::pair<std::vector<CombKalFilterBranch>, std::vector<double>> CombKalFilter::find_tracks() {
    //time_t start, end;
    //time(&start);
    auto start = high_resolution_clock::now();
    std::vector<CombKalFilterBranch> possibleTrackVec;
    std::vector<CombKalFilterBranch> trackVec;
    std::vector<childData> possibleTrackData;
    std::vector<std::array<double, 2>> trackSeedsUsed;
    std::vector<double> metaData;
    int minBinY;
    int minBinX;
    int maxBinY;
    int maxBinX;
    for (int p=0; p<startRange.size(); ++p) {
        binFinder(startRange[p][0][0], startRange[p][0][1], &minBinX, &maxBinX);
        binFinder(startRange[p][1][0], startRange[p][1][1], &minBinY, &maxBinY);
        for (int i = maxBinX; i >= minBinX; --i) {
            for (int j = maxBinY; j >= minBinY; --j) {
                for (int k = 0; k < pointsArray[i][j].size(); ++k) {
                    // ensure point is in array.
                    if (isin_coords(pointsArray[i][j][k], trackSeedsUsed)) {
                        continue;
                    }
                    std::vector<std::array<double, 2>> usedPoints;
                    std::vector<childData> newTrackData;
                    std::vector<CombKalFilterBranch> newTrackVec;
                    findSeeds(pointsArray[i][j][k], i, j, &usedPoints, &newTrackData, &newTrackVec);
                    possibleTrackData.insert(possibleTrackData.end(), newTrackData.begin(), newTrackData.end());
                    possibleTrackVec.insert(possibleTrackVec.end(), newTrackVec.begin(), newTrackVec.end());
                    trackSeedsUsed.insert(trackSeedsUsed.end(), usedPoints.begin(), usedPoints.end());
                }
            }
        }
    }
    //time(&end);
    finalMerge(&possibleTrackData, &possibleTrackVec, &trackVec);
    //double time_taken = double(end - start);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    //printf("%ld", duration.count());
    metaData.emplace_back(double (duration.count())/1e6);
    //metaData.emplace_back(time_taken);
    metaData.emplace_back(trackVec.size());
    return std::make_pair(trackVec, metaData);
}

bool CombKalFilter::similarity_check(std::vector<std::array<double, 2>> &child1,
                                     std::vector<std::array<double, 2>> &child2) {
    double nCommon = static_cast<double>(countCommonElements(&child1, &child2));
    double similar_points_threshold = 0.4;
    if ((nCommon/child1.size() > similar_points_threshold) or
        (nCommon/child2.size() > similar_points_threshold)){
        return true;
    }
    return false;
}

