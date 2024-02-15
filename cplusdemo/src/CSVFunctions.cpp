//
//
// Created by tomtr on 20/12/2023.
//

#include "CSVFunctions.h"


void writeCSV(std::string filename, std::vector<std::pair<std::string, std::vector<double>>> dataset){
    // Make a CSV file with one or more columns of integer values
    // Each column of data is represented by the pair <column name, column data>
    //   as std::pair<std::string, std::vector<int>>
    // The dataset is represented as a vector of these columns
    // Note that all columns should be the same size

    // Create an output filestream object
    std::ofstream myFile(filename, std::ios::app);
    std::streampos currentPosition = myFile.tellp();
    if (currentPosition == 0){
        // Send column names to the stream
        for(int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).first;
            if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }

    // Send data to the stream
    for(int i = 0; i < dataset.at(0).second.size(); ++i)
    {
        for(int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).second.at(i);
            if(j != dataset.size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }

    // Close the file
    myFile.close();
}

std::vector<std::pair<std::string, std::vector<double>>> readCSV(std::string filename){
    // Reads a CSV file into a vector of <string, vector<int>> pairs where
    // each pair represents <column name, column values>

    // Create a vector of <string, int vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<double>>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    double val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
        //line.erase(line.size() - 1);
        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){

            // Initialize and add <colname, int vector> pairs to result
            result.emplace_back(colname, std::vector<double> {});
        }
    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
        // Create a stringstream of the current line
        std::stringstream ss(line);

        // Keep track of the current column index
        int colIdx = 0;

        // Extract each integer
        while(std::getline(ss, colname, ',')){
            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(std::stod(colname));
            // If the next token is a comma, ignore it and move on
            // Increment the column index
            colIdx++;
        }
    }

    // Close file
    myFile.close();

    return result;
}