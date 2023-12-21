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
                                         int totalPoints, std::shared_ptr<CombKalFilter> ckf, int iteration): X(X), P(P),
                                         currBinX(currBinX), currBinY(currBinY), prevBinX(prevBinX), prevBinY(prevBinY),
                                         xPoints(x), yPoints(y), master(master), totalPoints(totalPoints),
                                         ckf(ckf), iteration(iteration){
}

void CombKalFilterBranch::residFunc(double x, double *maxResPoint){
    *maxResPoint = ckf->maxRes + 0;
    for (int i=0; i< P.rows(); ++i){
        *maxResPoint += abs(P(i,i) * pow(x, P.rows() - i -1));
    }
}

void CombKalFilterBranch::goodnessMeasure(double *measure){
    *measure = 0;
    for (int i=0; i< P.rows(); ++i){
        *measure += abs(pow(P(i,i),0.5) / X(i));
    }
}

void CombKalFilterBranch::fullPoints(std::vector<double> *fullXPoints,
                std::vector<double> *fullYPoints){
    int *a;
    int *b;
    a = nullptr;
    if (a == nullptr){
        int c;
        c =  ckf->fitOrder;
    }
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

void CombKalFilterBranch::getBestChild(std::shared_ptr<CombKalFilterBranch> bestchild, double *measure, int *childPoints){
    if (children.empty()){
        bestchild.reset(this);
        goodnessMeasure(measure);
        *childPoints = totalPoints;
    }
    else{
        *measure = 1e8;
        std::shared_ptr<CombKalFilterBranch> pc;
        double pcMeasure;
        int pcTPoints;
        for (int i= 0; children.size(); ++i){
            children[i].getBestChild(pc, &pcMeasure, &pcTPoints);
            if ((pcMeasure > *measure) & (pcTPoints >= ckf->minPoints)){
                *measure = pcMeasure;
                bestchild.reset(pc.get());
                *childPoints = pcTPoints;
            }
        }
        if (*measure  == 1e8){
            bestchild.reset(this);
            goodnessMeasure(measure);
            *childPoints = totalPoints;
        }
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
    RK(0,0) = ckf->RK;
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
    * res = abs((yVec - H*X)(0));
}

void CombKalFilterBranch::propogate(){
    if (iteration > ckf->maxIteration){
        return;
    }
    int magX;
    int magY;
    ckf->getCellDeltas(currBinX, prevBinX, &magX);
    ckf->getCellDeltas(currBinY, prevBinY, &magY);
    possibleNext.clear();

    for (int i=0; i<ckf->maxRange; ++i){
        std::vector<std::array<int, 2>> relativeCells = ckf->getNextCells(i, magX, magY);
        for (int j=0; j<relativeCells.size(); ++j){
            // get cell
            int nextBinX;
            int nextBinY;
            ckf->getNextBin(magX, relativeCells[j][0], currBinX, &nextBinX);
            ckf->getNextBin(magY, relativeCells[j][1], currBinY, &nextBinY);
            for (int k=0; k < ckf->pointsArray[nextBinX][nextBinY].size(); ++k) {
                double xP = ckf->pointsArray[nextBinX][nextBinY][k][0];
                double yP = ckf->pointsArray[nextBinX][nextBinY][k][1];
                if (ckf->isin_coords(xP, yP, xPoints, yPoints)) {
                    continue;
                }
                double res = 0;
                Eigen::VectorXd nX;
                Eigen::MatrixXd nP;
                double maxRes;
                residFunc(xP, &maxRes);
                updateKF(xP, yP, & nX, & nP, &res);
                if (abs(res) < maxRes){
                    // X is placeholder
                    possibleNext.push_back(nextPointData(xP, yP, nextBinX, nextBinY, nX, nP));
                }
            }
        }
        if (possibleNext.size() == 1){
            X = possibleNext[0].X;
            P = possibleNext[0].P;
            prevBinY = currBinX;
            prevBinX = currBinY;
            currBinX = possibleNext[0].binX;
            currBinY = possibleNext[0].binY;
            xPoints.push_back(possibleNext[0].x);
            xPoints.push_back(possibleNext[0].y);
            totalPoints +=1;
            iteration += 1;
            propogate();
        }
        else if (possibleNext.size() > 1){
            for (int i=0; i < possibleNext.size(); ++i){
                children.push_back(
                        CombKalFilterBranch(
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
                children[-1].propogate();
            }
        }
    }

}

//CombKalFilterBranch::CombKalFilterBranch(Eigen::VectorXd X, Eigen::MatrixXd P, int currBinX, int currBinY, int prevBinX, int prevBinY,
//                                         std::vector<double> x, std::vector<double> y, int totalPoints,
//                                         CombKalFilter *ckf, int iteration): X(X), P(P),
//currBinX(currBinX), currBinY(currBinY), prevBinX(prevBinX), prevBinY(prevBinY),
//xPoints(x), yPoints(y), totalPoints(totalPoints),
//ckf(ckf), iteration(iteration), master(this) {
//}
