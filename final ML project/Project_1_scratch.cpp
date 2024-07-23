//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <string>
//#include <algorithm>
//#include <unordered_map>
//#include <cmath>
//#include <random>
//#include "CSVProcessor.h"
//#include "Scalers.h"
//#include "GnuplotVisualizer.h"
//#include "StatMeasures.h"
//#include "CorrFeatureSelection.h"
//#include "KFold.h"
//#include "LogisticRegression.h"
//#include "KNN.h"
//#include "RFE_FeatureSelection.h"
//#include "Encoding.h"
//#include "TrainTestSplit.h"
//#include "EvaluationMetrics.h"
//#include "PCA.h"
//#include "NaiveBayes.h"
//
//using namespace std;
//
//
//void calculateConfusionMatrix(const std::vector<double>& y_true, const std::vector<double>& y_pred, int& TP, int& TN, int& FP, int& FN) {
//	TP = TN = FP = FN = 0;
//
//	for (size_t i = 0; i < y_true.size(); ++i) {
//		if (y_true[i] == 1 && y_pred[i] == 1) {
//			TP++;
//		}
//		else if (y_true[i] == 0 && y_pred[i] == 0) {
//			TN++;
//		}
//		else if (y_true[i] == 0 && y_pred[i] == 1) {
//			FP++;
//		}
//		else if (y_true[i] == 1 && y_pred[i] == 0) {
//			FN++;
//		}
//	}
//}
//
//
//
//int main()
//{
//
//	CSVProcessor csvProcessor;
//	Scalers scaler;
//	StatMeasures stat;
//	Encoding encode;
//
//
//
//	string filename = "E:\\ODC\\project\\no_attack_csv.csv";
//	vector<vector<string>> noAttack = csvProcessor.readCSV(filename);
//
//	filename = "E:\\ODC\\project\\attack_csv.csv";
//	vector<vector<string>> attack = csvProcessor.readCSV(filename);
//
//	//// count null values in each feature
//	cout << "Null values for each column:" << endl;
//	csvProcessor.countNulls(noAttack);
//	cout << endl;
//
//	cout << "attack " << endl;
//
//	cout << "Null values for each column:" << endl;
//	csvProcessor.countNulls(attack);
//	cout << endl;
//
//
//	vector <int> dropcolumns = { 25,24,23,22,18,17,0 };
//
//	for (auto columnindex : dropcolumns)
//	{
//		csvProcessor.removeColumn(attack, columnindex);
//		csvProcessor.removeColumn(noAttack, columnindex);
//	}
//
//
//	cout << "Null values for each column:" << endl;
//	csvProcessor.countNulls(noAttack);
//	cout << endl;
//
//	cout << "attack " << endl;
//
//	cout << "Null values for each column:" << endl;
//	csvProcessor.countNulls(attack);
//	cout << endl;
//
//	//csvProcessor.displayCSV(attack);
//
//	csvProcessor.fillNulls(attack, 5, "47095");
//	csvProcessor.fillNulls(attack, 7, "5469");
//	csvProcessor.fillNulls(noAttack, 5, "45360.5");
//	csvProcessor.fillNulls(noAttack, 7, "53");
//
//
//	cout << "Null values for each column:" << endl;
//	csvProcessor.countNulls(noAttack);
//	cout << endl;
//
//	cout << "attack " << endl;
//
//	cout << "Null values for each column:" << endl;
//	csvProcessor.countNulls(attack);
//	cout << endl;
//
//
//	vector<vector<string>> concdata = csvProcessor.appendCSV(attack, noAttack);
//
//	cout << "new shape : " << concdata.size() << " " << concdata[0].size() << endl;
//
//	csvProcessor.displayCSV(concdata);
//
//	// 2 3 4 6  10
//
//	vector <string> labels = csvProcessor.extractColumn(concdata, 26);
//	vector <double> y = csvProcessor.convertColumnToDouble(labels);
//
//
//	csvProcessor.removeColumn(concdata, 26);
//
//	vector<vector<string>> tempdata = concdata;
//
//	vector < int> encodeInedx = { 10 ,6,4,3,2 };
//
//	for (auto column : encodeInedx)
//	{
//		csvProcessor.removeColumn(concdata, column);
//	}
//
//	vector<vector<double>> doubleData = csvProcessor.convertToDouble(concdata);
//
//	for (int i = doubleData[0].size(); i >= 0; i--)
//	{
//		vector<double> scaledVec = scaler.standardScaler(csvProcessor.extractColumn(doubleData, i));
//		csvProcessor.removeColumn(doubleData, i);
//		csvProcessor.addColumn(doubleData, scaledVec);
//	}
//
//
//	for (auto column : encodeInedx)
//	{
//		vector<vector<double>> temp = encode.oneHotEncode(csvProcessor.extractColumn(tempdata, column));
//		csvProcessor.addColumns(doubleData, temp);
//
//	}
//
//	csvProcessor.displayCSV(doubleData);
////----------------------------------------------------------------------------------------------------------------
//	vector<vector<double>> X_Copy = doubleData;
//
//	int k = 4; // Number of folds
//
//	// Logistic Regression model
//	LogisticRegression lrModel_1(X_Copy[0].size());
//	KFold<LogisticRegression> kfoldLR(k, lrModel_1, true); // true for classification
//
//	cout << "Performing k-fold cross-validation with Logistic Regression (test case 1)" << endl;
//	kfoldLR.kFoldCrossValidation(X_Copy, y);
////----------------------------------------------------------------------------------------------------
//	vector<vector<double>> X_Copy2 = doubleData;
//	TrainTestSplit tts;
//	double testRatio = 0.3;
//	csvProcessor.addColumn(X_Copy2, y); // y added in column number 282
//	vector<vector<double>> trainData, testData;
//	tie(trainData, testData) = tts.trainTestSplit(X_Copy2, testRatio);
//
//	vector<double> y_train1 = csvProcessor.extractColumn(trainData, 282);
//	vector<double> y_test1 = csvProcessor.extractColumn(testData, 282);
//
//	csvProcessor.removeColumn(trainData, 282);
//	csvProcessor.removeColumn(testData, 282);
//
//	// remove the selected features from a copy of the dataset
//	vector<vector<double>> trainCopy1 = trainData;
//	vector<vector<double>> testCopy1 = testData;
//	
//	LogisticRegression lrModel_2(trainCopy1[0].size());
//	// Train the classifier
//	lrModel_2.fit(trainCopy1, y_train1);
//	
//	// Calculate accuracy for test case 2 (test data)
//	vector<double> predictions1;
//	for (const auto& sample : testCopy1) {
//	    double prediction1 = lrModel_2.predict(sample);
//	    predictions1.push_back(prediction1);
//	}
//	double acc1_test = EvaluationMetrics().accuracy(y_test1, predictions1);
//	cout << "Accuracy for test case 2 : " << acc1_test << endl;
//
//
//	int TP, TN, FP, FN;
//	calculateConfusionMatrix(y_test1, predictions1, TP, TN, FP, FN);
//
//	std::cout << "Confusion Matrix:" << std::endl;
//	std::cout << "TP: " << TP << std::endl;
//	std::cout << "TN: " << TN << std::endl;
//	std::cout << "FP: " << FP << std::endl;
//	std::cout << "FN: " << FN << std::endl;
//
//	return 0;
//}