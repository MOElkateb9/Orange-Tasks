#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include<cmath>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <cstdlib>
using namespace std;

typedef vector<vector<string>> DataFrame;


class dataOperations {
public:
    static DataFrame load(string filename)
    {
        DataFrame df;

        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error " << endl;
            return df;
        }

        string line;
        while (getline(file, line)) { // read line from file
            stringstream ss(line);
            string cell;
            vector<string> row;
            while (getline(ss, cell, ',')) {
                row.push_back(cell);
            }
            df.push_back(row);
        }

        return df;
    }



    static void print(const vector<vector<string>>& df, int numRows = 6, int numCols = -1) {
        if (df.empty()) return;

        if (numCols == -1)  numCols = df[0].size();

        vector<int> col_max_width(numCols, 0); // store max width / column

        // Find the maximum width for each column
        for (const vector<string>& current_row : df) {
            for (size_t i = 0; i < min(current_row.size(), static_cast<size_t>(numCols)); ++i) {
                if (current_row[i].size() > col_max_width[i]) {
                    col_max_width[i] = current_row[i].size();
                }
            }
        }

        int rowCount = 0;
        for (const vector<string>& row : df) {
            for (size_t i = 0; i < min(row.size(), static_cast<size_t>(numCols)); ++i) {
                int padding = col_max_width[i] - row[i].size();
                int left = padding / 2;
                int right = padding - left;

                for (int j = 0; j < left; ++j) {
                    cout << " ";
                }

                cout << row[i];

                for (int j = 0; j < right; ++j) {
                    cout << " ";
                }

                cout << " ";
            }
            cout << endl;
            ++rowCount;
            if (rowCount >= numRows) {
                break;
            }
        }
    }




    static void shape(DataFrame df)
    {
        if (df.empty())
            cout << 0 << " " << 0;
        cout << "rows:" << '(' << df.size() - 1 << ')' << " " << "cols:" << '(' << df[0].size() << ')' << endl; // exclude the header row
    }

    static vector<int> count_missing_values(DataFrame df) {
        vector<int> missing(df[0].size());

        for (int col = 0; col < df[0].size(); col++)
        {
            for (int row = 0; row < df.size(); row++)
            {
                if (df[row][col].empty())
                    missing[col]++;
            }
        }

        return  missing;
    }


    static void info(DataFrame df) {
        cout << "DataFrame Info:" << endl;
        cout << "Number of Rows: " << df.size() - 1 << endl;
        cout << "Number of Columns: " << df[0].size() << endl;


        vector<int> nonEmpty(df[0].size(), 0);
        for (int col = 0; col < df[0].size(); col++)
        {
            for (int row = 1; row < df.size(); row++)
            {
                if (!df[row][col].empty())
                {
                    nonEmpty[col]++;
                }
            }
        }


        cout << "Column Information:" << endl;
        for (int i = 0; i < df[0].size(); ++i)
        {
            cout << df[0][i] << ": ";
            cout << "Non-null Count: " << nonEmpty[i] << ", ";
            cout << endl;
        }
    }

    static void count_Duplicated(DataFrame df)
    {
        int duplicate = 0;
        set<string> seen;

        for (auto row : df)
        {
            string line = accumulate(row.begin(), row.end(), string()); //concate
            if (seen.find(line) != seen.end()) //return iter || seen.end()
            {
                duplicate++;
            }
            else
            {
                seen.insert(line);
            }
        }

        cout << duplicate << " dupliacted ";
    }


    static void remove_duplicated(DataFrame& df)
    {
        set<string> seen;
        df.erase(remove_if(df.begin(), df.end(), [&seen](vector<string> row) {
            string Row = accumulate(row.begin(), row.end(), string());
            return !seen.insert(Row).second;
            }), df.end()); /// using set auto

        count_Duplicated;
        cout << "are deleted" << endl;
    }



    static map<string, int> value_counts(DataFrame df, string column_name)
    {
        map<string, int> counts;
        int columnIndex = 0;
        for (int i = 0; i < df[0].size(); i++)
        {
            if (df[0][i] == column_name)
            {
                columnIndex = i;
                break;
            }
        }

        for (int r = 1; r < df.size(); r++) // Start from row 1 to skip header
        {
            counts[df[r][columnIndex]]++;
        }

        return counts;
    }


    static set<string> unique(DataFrame df, string column_name) {
        set<string> unique_Val;

        int columnIndex = 0;
        for (int i = 0; i < df[0].size(); i++)
        {
            if (df[0][i] == column_name)
            {
                columnIndex = i;
                break;
            }
        }


        for (int row = 1; row < df.size(); row++)
        {
            unique_Val.insert(df[row][columnIndex]);
        }

        return unique_Val;

    }




    static void replace(DataFrame& df, string columnName, string old_value, string new_value) {

        int col_index = 0;
        for (int i = 0; i < df[0].size(); i++) // Find the index of the column
        {
            if (df[0][i] == columnName)
            {
                col_index = i;
                break;
            }
        }

        for (int row = 1; row < df.size(); ++row) { // Start from row 1 to skip header
            if (df[row][col_index] == old_value)
                df[row][col_index] = new_value;

        }
    }

    static void renameColumn(DataFrame& df, const string& old_name, const string& new_name)
    {
        for (int i = 0; i < df[0].size(); i++)
        {
            if (df[0][i] == old_name)
            {
                df[0][i] = new_name;
                break;
            }
        }

    }

    template <typename T>
    static void export_to_csv(const T& df, const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file for writing" << endl;
            return;
        }

        for (const auto& row : df) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i < row.size() - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }

        file.close();
        if (file.fail()) {
            cerr << "Error closing file after writing" << endl;
        }
    }


    template <typename T>
    static void export_to_csv_vector(const vector<T>& vec, const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file for writing" << endl;
            return;
        }

        for (const auto& value : vec) {
            file << value << "\n";
        }

        file.close();
        if (file.fail()) {
            cerr << "Error closing file after writing" << endl;
        }
    }


    static void fillNulls(vector<vector<string>>& csvData, int columnIndex, const string& newValue) {
        if (columnIndex >= csvData[0].size()) {
            cout << "Invalid column index." << endl;
            return;
        }

        for (auto& row : csvData) {
            if (row.size() - 1 < columnIndex) {
                row.resize(columnIndex + 1);
                row[columnIndex] = newValue;
            }
            else if (row[columnIndex].empty()) {
                row[columnIndex] = newValue;
            }
        }
        cout << "Null values in column " << columnIndex << " filled with: " << newValue << endl;
    }

    template<typename T>
    static vector<vector<double>> convertToDouble(const vector<vector<T>>& data) {
        vector<vector<double>> doubleData;

        // Skip the first row assuming it's a header row
        bool firstRow = true;

        for (const auto& row : data) {
            if (firstRow) {
                firstRow = false;
                continue; // Skip the hroweader 
            }

            vector<double> doubleRow;
            for (const auto& cell : row) {
                try {
                    // Convert cell to string and then to double
                    double value = stod(static_cast<string>(cell));
                    doubleRow.push_back(value);
                }
                catch (const invalid_argument& e) {
                    cerr << "Invalid argument: " << cell << " cannot be converted to double" << endl;
                    // Handle invalid data here, like skipping or setting to a default value
                    // For example, you could push back a NaN or 0.0 to indicate an error
                    doubleRow.push_back(0.0); // Example: default to 0.0
                }
                catch (const out_of_range& e) {
                    cerr << "Out of range: " << cell << " is out of the range for a double" << endl;
                    // Handle out of range data here
                    // Example: push back a large or small value depending on context
                    doubleRow.push_back(0.0); // Example: default to 0.0
                }
            }
            doubleData.push_back(doubleRow);
        }
        return doubleData;
    }


    static vector<int> convertToInt(const vector<double>& data) {
        vector<int> intData;

        for (const auto& cell : data) {
            try {

                int value = static_cast<int>(round(cell));
                intData.push_back(value);
            }
            catch (const exception& e) {
                cerr << "Error: " << e.what() << " while converting " << cell << " to int" << endl;
                // Handle error by setting a default value
                intData.push_back(0); // Example: default to 0
            }
        }
        return intData;
    }

    template<typename T>
    static vector<T> extractColumn(vector<vector<T>>& csvData, int columnIndex) {
        vector<T> columnData;
        for (auto& row : csvData) {
            if (row.size() > columnIndex) {
                columnData.push_back(row[columnIndex]);
            }
        }
        return columnData;
    }

    template<typename T>
    static void dropColumn(vector<vector<T>>& csvData, int columnIndex) {
        for (auto& row : csvData) {
            if (row.size() > columnIndex) {
                row.erase(row.begin() + columnIndex); // Remove the column from the row
            }
        }
    }

    template<typename T>
    static void addColumn(vector<vector<T>>& csvData, const vector<T>& newColumn)
    {
        if (csvData.size() != newColumn.size()) {
            cout << "Error: Number of rows in CSV data does not match number of elements in new column." << endl;
            return;
        }

        for (int i = 0; i < csvData.size(); ++i) {
            csvData[i].push_back(newColumn[i]); // Append new element to each row
        }
    }

    static void remove_rows_with_value(DataFrame& df, const string& column_name, const string& target_value) {
        int col_index = -1;
        for (int i = 0; i < df[0].size(); i++) {
            if (df[0][i] == column_name) {
                col_index = i;
                break;
            }
        }

        if (col_index == -1) {
            cerr << "Column not found" << endl;
            return;
        }

        int initial_size = df.size();
        df.erase(remove_if(df.begin() + 1, df.end(), [&](vector<string>& row) {
            return row[col_index] == target_value;
            }), df.end());
        int final_size = df.size();

        cout << "column  (" << column_name << ')' << " with value: (" << target_value << ')' << " removed, " << "Number of rows deleted " << (initial_size - final_size) << endl;

    }

    static void encode_column(DataFrame& df, const string& column_name) {
        // Find column index
        int col_index = -1;
        for (int i = 0; i < df[0].size(); i++) {
            if (df[0][i] == column_name) {
                col_index = i;
                break;
            }
        }

        if (col_index == -1) {
            cerr << "Column not found" << endl;
            return;
        }

        // Create a mapping from categorical values to numerical values
        unordered_map<string, double> labelMap;
        double index = 0.0;

        for (int row = 1; row < df.size(); row++) {
            string& value = df[row][col_index];
            if (labelMap.find(value) == labelMap.end()) {
                labelMap[value] = index;
                index++;
            }
            value = to_string(labelMap[value]);
        }

        // Print the encoding map
        cout << "Encoding for column \"" << column_name << "\":\n";
        for (const auto& pair : labelMap) {
            cout << pair.first << " -> " << pair.second << endl;
        }
        cout << endl;
    }

    static void oneHotEncode(DataFrame& df, const string& column_name) {
        // Find column index
        int col_index = -1;
        for (int i = 0; i < df[0].size(); i++) {
            if (df[0][i] == column_name) {
                col_index = i;
                break;
            }
        }

        if (col_index == -1) {
            cerr << "Column not found" << endl;
            return;
        }

        // Create a mapping from categorical values to numerical indices
        unordered_map<string, int> labelMap;
        int index = 0;
        for (int row = 1; row < df.size(); row++) {
            const string& value = df[row][col_index];
            if (labelMap.find(value) == labelMap.end()) {
                labelMap[value] = index++;
            }
        }

        // Add new columns for one-hot encoding
        int original_col_size = df[0].size();
        for (const auto& pair : labelMap) {
            string new_col_name = column_name + "_" + pair.first;
            df[0].push_back(new_col_name);
        }

        // Initialize new columns with zeros
        for (int row = 1; row < df.size(); row++) {
            for (int i = original_col_size; i < df[0].size(); i++) {
                df[row].push_back("0");
            }
        }

        // Set ones for the corresponding categorical values
        for (int row = 1; row < df.size(); row++) {
            int col_pos = original_col_size + labelMap[df[row][col_index]];
            df[row][col_pos] = "1";
        }

        // Optionally, remove the original column
        // Uncomment the following lines if you want to remove the original column
        /*
        for (int row = 0; row < df.size(); row++) {
            df[row].erase(df[row].begin() + col_index);
        }
        */
    }


};

//=====================================================================================================

class StatisticsMeasures {
public:

    static double max(const vector<double>& column) {
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

    static double min(const vector<double>& column) {
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

    static double mean(const vector<double>& column) {
        double sum = 0.0;
        for (double value : column) {
            sum += value;
        }
        return sum / column.size();
    }

    static double median(vector<double> column) {
        sort(column.begin(), column.end());
        int col_size = column.size();
        if (col_size % 2 == 0) {
            return (column[col_size / 2 - 1] + column[col_size / 2]) / 2.0;
        }
        else {
            return column[col_size / 2];
        }
    }

    static double variance(const vector<double>& column) {
        double meanValue = mean(column);
        double sumOfSquares = 0.0;
        for (double value : column) {
            sumOfSquares += pow(value - meanValue, 2);
        }
        return sumOfSquares / column.size();
    }

    static double standardDeviation(const vector<double>& column) {
        return sqrt(variance(column));
    }

    static double covarianceBetweenVectors(const vector<double>& column1, const vector<double>& column2) {
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

    static double correlationBetweenVectors(const vector<double>& column1, const vector<double>& column2) {
        double cov = covarianceBetweenVectors(column1, column2);
        double stdDev1 = standardDeviation(column1);
        double stdDev2 = standardDeviation(column2);
        if (stdDev1 == 0 || stdDev2 == 0) {
            cout << "Standard deviation is zero, cannot compute correlation." << endl;
            return 0.0;
        }
        return cov / (stdDev1 * stdDev2);
    }

    static vector<vector<double>> covariance(vector<vector<double>>& dataset) {
        int numColumns = dataset[0].size();
        vector<vector<double>> covMatrix(numColumns, vector<double>(numColumns, 0.0));

        // Calculate covariance for each pair of columns
        for (int i = 0; i < numColumns; ++i) {
            for (int j = i; j < numColumns; ++j) {
                vector<double> column1 = dataOperations::extractColumn(dataset, i);
                vector<double> column2 = dataOperations::extractColumn(dataset, j);
                double cov = covarianceBetweenVectors(column1, column2);
                covMatrix[i][j] = cov;
                covMatrix[j][i] = cov;
            }
        }

        return covMatrix;
    }

    static vector<vector<double>> correlation(vector<vector<double>>& dataset) {
        int numColumns = dataset[0].size();
        vector<vector<double>> corrMatrix(numColumns, vector<double>(numColumns, 0.0));

        // Calculate correlation for each pair of columns
        for (int i = 0; i < numColumns; ++i) {
            for (int j = i; j < numColumns; ++j) {
                vector<double> column1 = dataOperations::extractColumn(dataset, i);
                vector<double> column2 = dataOperations::extractColumn(dataset, j);
                double corr = correlationBetweenVectors(column1, column2);
                corrMatrix[i][j] = corr;
                corrMatrix[j][i] = corr;
            }
        }

        return corrMatrix;
    }



};

////============================================================================================
//class Scalers {
//public:
//    static vector<double> minMaxScaler(const vector<double>& column) {
//        vector<double> scaledValues;
//        double maxVal = StatisticsMeasures::max(column);
//        double minVal = StatisticsMeasures::min(column);
//        double range = maxVal - minVal;
//
//        if (range == 0.0) {
//            cout << "Cannot perform Min-Max scaling. Range is zero." << endl;
//            return scaledValues;
//        }
//
//        for (double value : column) {
//            double scaledValue = (value - minVal) / range;
//            scaledValues.push_back(scaledValue);
//        }
//
//        return scaledValues;
//    }
//
//
//    static vector<vector<double>> minMaxScalerdf(const vector<vector<double>>& data) {
//        if (data.empty() || data[0].empty()) {
//            cout << "Input data is empty." << endl;
//            return {};
//        }
//
//        size_t numRows = data.size();
//        size_t numCols = data[0].size();
//        vector<vector<double>> scaledData(numRows, vector<double>(numCols));
//
//        for (size_t col = 0; col < numCols; ++col) {
//            double maxVal = -numeric_limits<double>::infinity();
//            double minVal = numeric_limits<double>::infinity();
//
//            // Find max and min values in the column
//            for (size_t row = 0; row < numRows; ++row) {
//                double value = data[row][col];
//                if (value > maxVal) maxVal = value;
//                if (value < minVal) minVal = value;
//            }
//
//            double range = maxVal - minVal;
//
//            if (range == 0.0) {
//                cout << "Cannot perform Min-Max scaling. Range is zero for column " << col << endl;
//                for (size_t row = 0; row < numRows; ++row) {
//                    scaledData[row][col] = data[row][col]; // Preserve the original values
//                }
//                continue;
//            }
//
//            // Scale the column
//            for (size_t row = 0; row < numRows; ++row) {
//                double value = data[row][col];
//                double scaledValue = (value - minVal) / range;
//                scaledData[row][col] = scaledValue;
//            }
//        }
//
//        return scaledData;
//    }
//
//
//    static vector<double> standardScaler(const vector<double>& column) {
//        vector<double> scaledValues;
//        double meanVal = StatisticsMeasures::mean(column);
//        double stdDev = StatisticsMeasures::standardDeviation(column);
//        if (stdDev == 0.0) {
//            cout << "Cannot perform standardization. Standard deviation is zero." << endl;
//            return scaledValues;
//        }
//
//        for (double value : column) {
//            double scaledValue = (value - meanVal) / stdDev;
//            scaledValues.push_back(scaledValue);
//        }
//
//        return scaledValues;
//    }
//
//    static vector<vector<double>> standardScalerdf(const vector<vector<double>>& data) {
//        vector<vector<double>> scaledData;
//        scaledData.resize(data.size(), vector<double>(data[0].size()));
//
//        for (size_t col = 0; col < data[0].size(); ++col) {
//            vector<double> column;
//            for (size_t row = 0; row < data.size(); ++row) {
//                column.push_back(data[row][col]);
//            }
//
//            double meanVal = StatisticsMeasures::mean(column);
//            double stdDev = StatisticsMeasures::standardDeviation(column);
//
//            if (stdDev == 0.0) {
//                cout << "Cannot perform standardization on column " << col << ". Standard deviation is zero." << endl;
//                continue;
//            }
//
//            for (size_t row = 0; row < data.size(); ++row) {
//                scaledData[row][col] = (data[row][col] - meanVal) / stdDev;
//            }
//        }
//
//        return scaledData;
//    }
//
//
//
//};



//============================================================================================
class outliers
{
public:

    static double percentile(vector<double>& data, double percent) {
        sort(data.begin(), data.end());
        int n = data.size();
        double rank = percent * (n - 1);
        int lower_rank = floor(rank);
        int upper_rank = ceil(rank);

        if (lower_rank == upper_rank)
        {
            return data[lower_rank];
        }
        else {
            double d0 = data[lower_rank] * (upper_rank - rank);
            double d1 = data[upper_rank] * (rank - lower_rank);
            return d0 + d1;
        }
    }

    static vector<double> detect_outliers(DataFrame& df, const string& column_name) {
        // Find column index
        int col_index = -1;
        for (int i = 0; i < df[0].size(); i++) {
            if (df[0][i] == column_name) {
                col_index = i;
                break;
            }
        }

        if (col_index == -1) {
            cerr << "Column not found" << endl;
            return {};
        }

        // Extract numerical values from the column
        vector<double> column_data;
        for (int row = 1; row < df.size(); row++) {
            double value = stod(df[row][col_index]);
            column_data.push_back(value);
        }

        // Calculate quartiles and IQR
        sort(column_data.begin(), column_data.end());
        double Q1 = percentile(column_data, 0.25);
        double Q3 = percentile(column_data, 0.75);
        double IQR = Q3 - Q1;

        // Determine outliers
        double lower_bound = Q1 - 1.5 * IQR;
        double upper_bound = Q3 + 1.5 * IQR;

        // Collect outliers
        vector<double> outliers;
        for (double value : column_data) { if (value < lower_bound || value > upper_bound) { outliers.push_back(value); } }
        if (!outliers.empty())  cout << "no outliers \n";
        else  return outliers;
    }


    static void remove_outliers(DataFrame& df, const string& column_name) {
        // Detect outliers
        vector<double> outliers = detect_outliers(df, column_name);

        // Find column index
        int col_index = -1;
        for (int i = 0; i < df[0].size(); i++) {
            if (df[0][i] == column_name) {
                col_index = i;
                break;
            }
        }

        if (col_index == -1) {
            cerr << "Column not found" << endl;
            return;
        }

        // Remove rows containing outliers
        df.erase(remove_if(df.begin() + 1, df.end(), [&](const vector<string>& row) {
            double value = stod(row[col_index]);
            return find(outliers.begin(), outliers.end(), value) != outliers.end();
            }), df.end());

        // Print number of rows removed
        cout << "Removed " << outliers.size() << " rows containing outliers in column \"" << column_name << "\"" << endl;
    }

    static void replace_outliers_with_mean(DataFrame& df, const string& column_name, const vector<double>& outliers) {
        // Find column index
        int col_index = -1;
        for (int i = 0; i < df[0].size(); i++)
        {
            if (df[0][i] == column_name)
            {
                col_index = i;
                break;
            }
        }

        if (col_index == -1) {
            cerr << "Column not found" << endl;
            return;
        }

        // Calculate mean of the column excluding outliers
        vector<double> column_data;
        for (int row = 1; row < df.size(); row++) {
            double value = stod(df[row][col_index]);
            if (find(outliers.begin(), outliers.end(), value) == outliers.end()) {
                column_data.push_back(value);
            }
        }
        double mean = accumulate(column_data.begin(), column_data.end(), 0.0) / column_data.size();

        // Replace outliers with mean
        int replaced_count = 0;
        for (int row = 1; row < df.size(); row++)
        {
            double value = stod(df[row][col_index]);
            if (find(outliers.begin(), outliers.end(), value) != outliers.end())
            {
                df[row][col_index] = to_string(mean);
                replaced_count++;
            }
        }

        // Print number of rows modified
        cout << "Replaced " << replaced_count << " outliers in column \"" << column_name << "\" with mean value " << mean << endl;
    }
};
