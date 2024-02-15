//
// Created by tomtr on 20/12/2023.
//

#ifndef CPLUSDEMO_CSVFUNCTIONS_H
#define CPLUSDEMO_CSVFUNCTIONS_H
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <sstream>
#include "iostream"

void writeCSV(std::string filename, std::vector<std::pair<std::string, std::vector<double>>> dataset);

std::vector<std::pair<std::string, std::vector<double>>> readCSV(std::string filename);


#endif //CPLUSDEMO_CSVFUNCTIONS_H
