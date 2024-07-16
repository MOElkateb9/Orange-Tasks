#ifndef SCALERS_H
#define SCALERS_H

#include <vector>
#include <iostream>
#include "StatMeasures.h"

using namespace std;

class Scalers {
public:
    vector<double> minMaxScaler(const vector<double>& column);
    vector<double> standardScaler(const vector<double>& column);
};

// Member function definitions

vector<double> Scalers::minMaxScaler(const vector<double>& column) {
    vector<double> scaledValues;
    StatMeasures statMeasures;

    double maxVal = statMeasures.maxValue(column);
    double minVal = statMeasures.minValue(column);
    double range = maxVal - minVal;

    if (range == 0.0) {
        cout << "Cannot perform Min-Max scaling. Range is zero." << endl;
        return scaledValues;
    }

    for (double value : column) {
        double scaledValue = (value - minVal) / range;
        scaledValues.push_back(scaledValue);
    }

    return scaledValues;
}

vector<double> Scalers::standardScaler(const vector<double>& column) {
    vector<double> scaledValues;
    StatMeasures statMeasures;

    double meanVal = statMeasures.mean(column);
    double stdDev = statMeasures.standardDeviation(column);
    if (stdDev == 0.0) {
        cout << "Cannot perform standardization. Standard deviation is zero." << endl;
        return scaledValues;
    }

    for (double value : column) {
        double scaledValue = (value - meanVal) / stdDev;
        scaledValues.push_back(scaledValue);
    }

    return scaledValues;
}

#endif // SCALERS_H
