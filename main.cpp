#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <unordered_map>
#include <algorithm>
#include <numeric>

// Helper function to trim trailing newline characters
bool canConvertToInt(const std::string& str) {
    size_t i = (str[0] == '-' || str[0] == '+') ? 1 : 0; // Skip sign
    return !str.empty() && std::all_of(str.begin() + i, str.end(), ::isdigit);
}

int main() {
    // Open and process the first file (new1.txt)
    std::ifstream dataFile("info.txt");
    if (!dataFile.is_open()) {
        std::cerr << "Error: Could not open info.txt" << std::endl;
        return 1;
    }
    std::ofstream outputFile("Compound_Info.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 1;
    }

    std::vector<std::string> useful;
    std::string dataLine;
    while (std::getline(dataFile, dataLine)) {
        if (dataLine.find("True True") != std::string::npos || dataLine.find("\\MS24GE_") != std::string::npos) {
            useful.push_back(dataLine);
        }
    }
    dataFile.close();

    // Open and process the formula file (formula.txt)
    std::ifstream formulaFile("formula.txt");
    if (!formulaFile.is_open()) {
        std::cerr << "Error: Could not open formula.txt" << std::endl;
        return 1;
    }

    std::unordered_map<std::string, std::string> formula;
    std::string formulaLine;
    while (std::getline(formulaFile, formulaLine)) {
        // Find the position of the last space
        size_t lastSpace = formulaLine.find_last_of(' ');

        std::string key = formulaLine.substr(0, lastSpace); // Everything before the last space
        std::string value = formulaLine.substr(lastSpace + 1); // Everything after the last space

        formula[key] = value;
    }
    formulaFile.close();

    // Extract address
    std::string address = useful[0].substr(useful[0].find_last_of('\\') + 1);
    address = address.substr(0, address.find('.'));

    // Extract names and methods
    std::vector<std::string> names;
    std::vector<std::string> methods;

for (size_t i = 1; i < useful.size(); ++i) {
    // Step 1: Extract the part of the string before "Group"
    std::string str = useful[i].substr(0, useful[i].find(" Group"));

    // Step 2: Find the last space in the extracted string to get the method
    size_t lastSpace = str.find_last_of(' ');
    if (lastSpace == std::string::npos) {
        std::cerr << "Warning: Malformed line, no space found: " << str << std::endl;
        continue;
    }

    // Step 3: Extract the method (last word in the remaining string)
	
    std::string method = str.substr(lastSpace + 1);
    methods.push_back(method);

    // Step 4: Extract the name by removing the method from the string
    std::string name = str.substr(0, lastSpace);
    names.push_back(name);
}
	std::string zapeta="\'";
    int count = 1;
    int count_previous = 2;
outputFile <<   "Filename	Experiment	CH$NAME:	CH$FORMULA:	AC$MASS_SPECTROMETRY: FRAGMENTATION_TYPE" << std::endl;
    for (size_t ind = 0; ind < names.size() - 1; ++ind) {
        count++;
        if ((canConvertToInt(methods[ind]) ? names[ind] : names[ind] + methods[ind]) != 
        (canConvertToInt(methods[ind + 1]) ? names[ind + 1] : names[ind + 1] + methods[ind + 1])) {
        if (count_previous == count) {
                outputFile   << address << "\t" << zapeta << count_previous << "\t" << names[ind] << "\t" << formula[names[ind]] << "\t" << methods[ind] << std::endl;
			} else {
                outputFile   << address << "\t" <<  zapeta << count_previous << "-" << count << "\t" << names[ind] << "\t" << formula[names[ind]] << "\t" << methods[ind] << std::endl;
			}
			 zapeta="";

            count_previous = count + 1;
        }
    }

    size_t ind = names.size() - 1;
    if (count_previous == count + 1) {
        outputFile   << address << "\t" << count_previous << "\t" << names[ind] << "\t" << formula[names[ind]] << "\t" << methods[ind] << std::endl;
    } else {
        outputFile   << address << "\t" << count_previous << "-" << count + 1 << "\t" << names[ind] << "\t" << formula[names[ind]] << "\t" << methods[ind] << std::endl;
    }
    outputFile.close();

    return 0;
}
