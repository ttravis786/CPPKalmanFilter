//
// Created by tomtr on 02/12/2023.
//


#include "CombKalFilterBranch.h"

nextPointData::nextPointData(double x, double y, int binX, int binY, Eigen::MatrixXd X, Eigen::MatrixXd P):
x(x), y(y), binY(binY), binX(binX), X(X), P(P) {
}

CombKalFilterBranch::CombKalFilterBranch(Eigen::VectorXd X, Eigen::MatrixXd P, int currBinX, int currBinY, int prevBinX,
                                         int prevBinY, std::vector<double> x, std::vector<double> y,
                                         CombKalFilterBranch *master,
                                         int totalPoints, CombKalFilter * ckf, int iteration): X(X), P(P),
                                         currBinX(currBinX), currBinY(currBinY), prevBinX(prevBinX), prevBinY(prevBinY),
                                         xPoints(x), yPoints(y), master(master), totalPoints(totalPoints),
                                         ckf(ckf), iteration(iteration){
}

void CombKalFilterBranch::residFunc(double x, double *maxResPoint){
    *maxResPoint = ckf->maxRes + 0;
    for (int i=0; i< P.rows(); ++i){
        *maxResPoint += abs(pow(P(i,i), 0.5) * pow(x, P.rows() - i -1));
    }
    int b = 1;
}

void CombKalFilterBranch::goodnessMeasure(double *measure){
    *measure = 0;
    for (int i=0; i< P.rows(); ++i){
        *measure += abs(pow(P(i,i),0.5) / X(i));
    }
}

void CombKalFilterBranch::fullPoints(std::vector<double> *fullXPoints,
                                     std::vector<double> *fullYPoints){

    if (master != nullptr){
        master->fullPoints(fullXPoints, fullYPoints);
    }
    fullYPoints->insert(fullYPoints->end(),
                        yPoints.begin(),
                        yPoints.end());
    fullXPoints->insert(fullXPoints->end(),
                        xPoints.begin(),
                        xPoints.end());
}

std::vector<std::array<double, 2>> CombKalFilterBranch::pointConverter(std::vector<double> *x, std::vector<double> *y) {
    std::vector<std::array<double, 2>> pointVector;
    for (int i = 0; i <x->size(); ++i){
        pointVector.emplace_back(std::array<double, 2> {(*x)[i], (*y)[i]});
    }
    return pointVector;
}



void CombKalFilterBranch::startChild(std::vector<childData> * bestchildren){
    double measure_threshold = 0.2;
    //bestchild.reset(this);
    childData bestChild;
    bestChild.branch = this;
    bestChild.measure = 0;
    goodnessMeasure(&bestChild.measure);
    // will get rid of point converter later
    bestChild.points = pointConverter(&xPoints, &yPoints);
    bestChild.totalPoints = totalPoints;
    if ((bestChild.measure < measure_threshold) & (bestChild.totalPoints >= ckf->minPoints)){
        bestchildren->emplace_back(bestChild);
    }
}


void CombKalFilterBranch::getBestChild(std::vector<childData> * bestchildren){
    if (children.empty()){
        startChild(bestchildren);
    }
    // get single best child
    else {
        if (not ckf->multiple){
            double pcMeasure = 1e8;
            childData bestChild;
            std::vector<childData> allChildren;
            for (int i = 0; i < children.size(); ++i) {
                std::vector<childData> pcv;
                (*children[i]).getBestChild(&pcv);
                for (int j=0; j< pcv.size(); ++j){
                    if (pcv[j].measure < pcMeasure){
                        bestChild = pcv[j];
                        pcMeasure = pcv[j].measure;
                    }
                }
            }
            if (pcMeasure != 1e8){
                bestchildren->emplace_back(bestChild);
            }
        }
        // if contain less than 20 % of points in common of smaller seed then they are seperate.
        //std::shared_ptr<CombKalFilterBranch> pc;
        else {
            std::vector<childData> allChildren;
            for (int i = 0; i < children.size(); ++i) {
                std::vector<childData> pc1;
                (*children[i]).getBestChild(&pc1);
                allChildren.insert(allChildren.end(), pc1.begin(), pc1.end());
            }
            // remove those that need to be merged.
            if (ckf->remove) {
                std::vector<int> removeIndexes;
                for (int i = 0; i < allChildren.size(); i++) {
                    for (int j = i + 1; j < allChildren.size(); j++) {
                        if (ckf->similarity_check(allChildren[i].points, allChildren[j].points)) {
                            removeIndexes.emplace_back(j);
                            if (allChildren[i].points.size() < allChildren[j].points.size()) {
                                allChildren[i] = allChildren[j];
                            }
                        }
                    }
                }
                bool include;
                for (int i = 0; i < allChildren.size(); i++) {
                    include = true;
                    for (int j = 0; j < removeIndexes.size(); ++j) {
                        if (removeIndexes[j] == i) {
                            include = false;
                            break;
                        }
                    }
                    if (include) {
                        bestchildren->push_back(allChildren[i]);
                    }
                }
            }
            // can chose to keep them all.
            else {
                *bestchildren = allChildren;
            }
        }
        std::vector<std::array<double, 2>> newPoints = pointConverter(&xPoints, &yPoints);
        if (bestchildren->empty()){
            startChild(bestchildren);
        }
        else{
            for (int i = 0; i < bestchildren->size(); ++i) {
                // add new points, new to change this.
                (*bestchildren)[i].points.insert((*bestchildren)[i].points.begin(), newPoints.begin(), newPoints.end());
            }
        }
//        std::vector<childData> pc1;
//        std::vector<int> removeIndexes;
//        std::vector<childData> pc2;
//        (*children[0]).getBestChild(&pc1);
//        for (int i=1; i < children.size(); ++i){
//            (*children[i]).getBestChild(&pc2);
//            // only check orginal ones
//            // I simplified the algo here to fix a bug, should produce the same result but less efficient.
//            for (int j=0; j < pc1.size(); ++j){
//                for (int k=0; k < pc2.size(); ++k){
//                    // check for similarity between pc1[j] and pc2[j], merge them if there is.
//                    if (ckf->similarity_check(pc1[j].points, pc2[k].points)) {
//                        // currently pick one with more points, maybe make more sophisticated.
//                        if (pc2[j].totalPoints > pc1[k].totalPoints) {
//                            pc1[j] = pc2[k];
//                        }
//                    }
//                    else{
//                        pc1.emplace_back(pc2[k]);
//                    }
//                }
//            }
//        }

        // I simplified the algo here to fix a bug, should produce the same result but less efficient.
    }
}


void CombKalFilterBranch::updateKF(double xP,double yP, Eigen::VectorXd * nX, Eigen::MatrixXd * nP, double *res){
    Eigen::VectorXd X_k_km1 = ckf->F * X;
    Eigen::MatrixXd P_k_km1 = ckf->F * P;
    // H should really be a vector in our context
    Eigen::MatrixXd H = ckf->getH(xP);
    // this stuff will be fixed
    Eigen::VectorXd yVec = Eigen::VectorXd::Zero(P.rows());
    Eigen::MatrixXd RK = Eigen::MatrixXd::Zero(P.rows(), P.rows());
    Eigen::MatrixXd invSK = Eigen::MatrixXd::Zero(P.rows(), P.rows());
    RK(0,0) = ckf->RK / 2;
    yVec[0] = yP;
    //
    Eigen::VectorXd yK = yVec - H * X_k_km1;
    Eigen::MatrixXd SK = (H * P_k_km1) * H.transpose() + RK;
    // fix this to
    invSK(0,0) = 1/SK(0,0);
    //
    Eigen::MatrixXd K = P_k_km1 * H.transpose() * invSK;
    * nX = X_k_km1 + K * yK;
    * nP = (Eigen::MatrixXd::Identity(P.rows(), P.rows()) - K * H) * P_k_km1;
    * res = yK[0];
    //* res = abs((yVec - H**nX)(0));
}

bool CombKalFilterBranch::suitableNextBin(int bin){
    if (bin > 111 or bin < 0){
        return false;
    }
    return true;
}

void CombKalFilterBranch::deleteBranches(){
    if (children.empty()){
        return;
    }
        // get single best child
    else {
        for (int i = 0; i < children.size(); ++i) {
            // delete all it's children
            (*children[i]).deleteBranches();
            // then delete it
            delete children[i];
        }
    }
}

void CombKalFilterBranch::propogate(){
    if (iteration > ckf->maxIteration){
        return;
    }
    //printf("%d \n", iteration);
    int magX;
    int magY;
    if (ckf->fixedView){
        magX = -1;
        magY = 0;
    }
    else {
        ckf->getCellDeltas(currBinX, prevBinX, &magX);
        ckf->getCellDeltas(currBinY, prevBinY, &magY);
    }
    possibleNext.clear();
    std::vector<double> fullXPoints;
    std::vector<double> fullYPoints;
    fullPoints(&fullXPoints,
               &fullYPoints);

//    if (std::fabs((fullXPoints.back() - -12.75) < 0.1) and std::fabs((fullYPoints.back() - -16) < 0.1)){
//        int a = 1;
//    }

    for (int i=0; i<ckf->maxRange; ++i){
        std::vector<std::array<int, 2>> relativeCells = ckf->getNextCells(i, magX, magY);
        for (int j=0; j<relativeCells.size(); ++j) {
            // get cell
            int nextBinX;
            int nextBinY;
            ckf->getNextBin(magX, relativeCells[j][0], currBinX, &nextBinX);
            ckf->getNextBin(magY, relativeCells[j][1], currBinY, &nextBinY);
            if (not (suitableNextBin(nextBinX) and suitableNextBin(nextBinX))){
                continue;
            }
            for (int k = 0; k < ckf->pointsArray[nextBinX][nextBinY].size(); ++k) {
//                if (i==1){
//                    break;
//                }
                //printf("%d \n", k);
                double xP = ckf->pointsArray[nextBinX][nextBinY][k][0];
                double yP = ckf->pointsArray[nextBinX][nextBinY][k][1];

                if (ckf->isin_coords(xP, yP, fullXPoints, fullYPoints)) {
                    continue;
                }

                double res = 0;
                Eigen::VectorXd nX;
                Eigen::MatrixXd nP;
                double maxRes;
                residFunc(xP, &maxRes);

//                printf("%f \n", maxRes);
//                std::cout << "My double value: " << maxRes << std::endl;
//                if (std::fabs(maxRes - 0.447790) < 0.0001){
//                    int a = 1;
//                }
                updateKF(xP, yP, &nX, &nP, &res);
                if (std::fabs(res) < maxRes) {
                    // X is placeholder
                    possibleNext.push_back(nextPointData(xP, yP, nextBinX, nextBinY, nX, nP));
                }
            }
        }
//        fullXPoints.swap(v);
//        // Release the memory
//        fullXPoints.shrink_to_fit();
//        fullYPoints.swap(v);
//        fullYPoints.clear();
//        // Release the memory
//        fullYPoints.shrink_to_fit();
        if (possibleNext.size() == 1){
            X = possibleNext[0].X;
            P = possibleNext[0].P;
            prevBinY = currBinY;
            prevBinX = currBinX;
            currBinX = possibleNext[0].binX;
            currBinY = possibleNext[0].binY;
            xPoints.push_back(possibleNext[0].x);
            yPoints.push_back(possibleNext[0].y);
            totalPoints +=1;
//            iteration += 1;
            std::vector<double>().swap(fullXPoints);
            std::vector<double>().swap(fullYPoints);
            propogate();
            return;
        }
        else if (possibleNext.size() > 1){
            for (int i=0; i < possibleNext.size(); ++i){
                //printf("Obj(%d): %p\n",i,this);
                children.push_back(
                        new CombKalFilterBranch(
                                possibleNext[i].X,
                                possibleNext[i].P,
                                possibleNext[i].binX,
                                possibleNext[i].binY,
                                currBinX,
                                currBinY,
                                std::vector<double> {possibleNext[i].x},
                                std::vector<double> {possibleNext[i].y},
                                this,
                                totalPoints + 1,
                                ckf,
                                iteration + 1
                                ));
                (*children[i]).propogate();
            }
            std::vector<double>().swap(fullXPoints);
            std::vector<double>().swap(fullYPoints);
            return;
        }
    }
    std::vector<double>().swap(fullXPoints);
    std::vector<double>().swap(fullYPoints);
}

//CombKalFilterBranch::CombKalFilterBranch(Eigen::VectorXd X, Eigen::MatrixXd P, int currBinX, int currBinY, int prevBinX, int prevBinY,
//                                         std::vector<double> x, std::vector<double> y, int totalPoints,
//                                         CombKalFilter *ckf, int iteration): X(X), P(P),
//currBinX(currBinX), currBinY(currBinY), prevBinX(prevBinX), prevBinY(prevBinY),
//xPoints(x), yPoints(y), totalPoints(totalPoints),
//ckf(ckf), iteration(iteration), master(this) {
//}
