#ifndef ENCODING_H
#define ENCODING_H

#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class Encoding {
public:
    vector<double> labelEncode(const vector<string>& labels);
    vector<vector<double>> oneHotEncode(const vector<string>& labels);
    vector<double> mapEncode(const vector<string>& labels, const unordered_map<string, double>& valueMap);
};

// Member function definitions

vector<double> Encoding::labelEncode(const vector<string>& labels) {
    unordered_map<string, double> labelMap;
    vector<double> encoded;
    double index = 0.0;
    for (const string& label : labels) {
        if (labelMap.find(label) == labelMap.end()) {
            labelMap[label] = index;
            index++;
        }
        encoded.push_back(labelMap[label]);
    }

    return encoded;
}

vector<vector<double>> Encoding::oneHotEncode(const vector<string>& labels) {
    unordered_map<string, int> labelMap;
    int index = 0;
    for (const string& label : labels) {
        if (labelMap.find(label) == labelMap.end()) {
            labelMap[label] = index++;
        }
    }

    vector<vector<double>> encoded(labels.size(), vector<double>(labelMap.size(), 0.0));
    for (int i = 0; i < labels.size(); ++i) {
        encoded[i][labelMap[labels[i]]] = 1.0;
    }

    return encoded;
}

vector<double> Encoding::mapEncode(const vector<string>& labels, const unordered_map<string, double>& valueMap) {
    vector<double> encoded;
    for (const string& label : labels) {
        if (valueMap.find(label) != valueMap.end()) {
            encoded.push_back(valueMap.at(label));
        }
        else {
            // Assign -1 if label is not found in the map
            encoded.push_back(-1.0);
        }
    }
    return encoded;
}

#endif // ENCODING_H
