#ifndef STATMEASURES_H
#define STATMEASURES_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "CSVProcessor.h"

using namespace std;

class StatMeasures {
public:
    double maxValue(const vector<double>& column);
    double minValue(const vector<double>& column);
    double mean(const vector<double>& column);
    double median(vector<double> column);
    double variance(const vector<double>& column);
    double standardDeviation(const vector<double>& column);
    double covarianceBetweenVectors(const vector<double>& column1, const vector<double>& column2);
    double correlationBetweenVectors(const vector<double>& column1, const vector<double>& column2);
    vector<vector<double>> covariance(vector<vector<double>>& dataset);
    vector<vector<double>> correlation(vector<vector<double>>& dataset);
};

// Member function definitions

double StatMeasures::maxValue(const vector<double>& column) {
    if (column.empty()) {
        cout << "Empty column. Cannot find maximum value." << endl;
        return 0.0;
    }

    double maxVal = column[0];
    for (double value : column) {
        if (value > maxVal) {
            maxVal = value;
        }
    }
    return maxVal;
}

double StatMeasures::minValue(const vector<double>& column) {
    if (column.empty()) {
        cout << "Empty column. Cannot find minimum value." << endl;
        return 0.0;
    }

    double minVal = column[0];
    for (double value : column) {
        if (value < minVal) {
            minVal = value;
        }
    }
    return minVal;
}

double StatMeasures::mean(const vector<double>& column) {
    double sum = 0.0;
    for (double value : column) {
        sum += value;
    }
    return sum / column.size();
}

double StatMeasures::median(vector<double> column) {
    sort(column.begin(), column.end());
    size_t col_size = column.size();
    if (col_size % 2 == 0) {
        return (column[col_size / 2 - 1] + column[col_size / 2]) / 2.0;
    }
    else {
        return column[col_size / 2];
    }
}

double StatMeasures::variance(const vector<double>& column) {
    double meanValue = mean(column);
    double sumOfSquares = 0.0;
    for (double value : column) {
        sumOfSquares += pow(value - meanValue, 2);
    }
    return sumOfSquares / column.size();
}

double StatMeasures::standardDeviation(const vector<double>& column) {
    return sqrt(variance(column));
}

double StatMeasures::covarianceBetweenVectors(const vector<double>& column1, const vector<double>& column2) {
    double mean1 = mean(column1);
    double mean2 = mean(column2);
    double sum = 0.0;
    int n = column1.size();
    if (n != column2.size()) {
        cout << "Columns must have the same size for covariance calculation." << endl;
        cout << n << " " << column2.size() << endl;
        return 0.0;
    }
    for (int i = 0; i < n; ++i) {
        sum += (column1[i] - mean1) * (column2[i] - mean2);
    }
    return sum / n;
}

double StatMeasures::correlationBetweenVectors(const vector<double>& column1, const vector<double>& column2) {
    double cov = covarianceBetweenVectors(column1, column2);
    double stdDev1 = standardDeviation(column1);
    double stdDev2 = standardDeviation(column2);
    if (stdDev1 == 0 || stdDev2 == 0) {
        cout << "Standard deviation is zero, cannot compute correlation." << endl;
        return 0.0;
    }
    return cov / (stdDev1 * stdDev2);
}

vector<vector<double>> StatMeasures::covariance(vector<vector<double>>& dataset) {
    int numColumns = dataset[0].size();
    vector<vector<double>> covMatrix(numColumns, vector<double>(numColumns, 0.0));

    CSVProcessor csvProcessor;

    for (int i = 0; i < numColumns; ++i) {
        for (int j = i; j < numColumns; ++j) {
            vector<double> column1 = csvProcessor.extractColumn(dataset, i);
            vector<double> column2 = csvProcessor.extractColumn(dataset, j);
            double cov = covarianceBetweenVectors(column1, column2);
            covMatrix[i][j] = cov;
            covMatrix[j][i] = cov;
        }
    }

    return covMatrix;
}

vector<vector<double>> StatMeasures::correlation(vector<vector<double>>& dataset) {
    int numColumns = dataset[0].size();
    vector<vector<double>> corrMatrix(numColumns, vector<double>(numColumns, 0.0));

    CSVProcessor csvProcessor;

    for (int i = 0; i < numColumns; ++i) {
        for (int j = i; j < numColumns; ++j) {
            vector<double> column1 = csvProcessor.extractColumn(dataset, i);
            vector<double> column2 = csvProcessor.extractColumn(dataset, j);
            double corr = correlationBetweenVectors(column1, column2);
            corrMatrix[i][j] = corr;
            corrMatrix[j][i] = corr;
        }
    }

    return corrMatrix;
}

#endif // STATMEASURES_H
