#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <armadillo>
#include <stdexcept> 
#include <limits>    
#include <cctype>    
#include <mlpack/core.hpp>
#include <mlpack/mlpack.hpp>
#include <mlpack/core/cv/k_fold_cv.hpp>
#include <mlpack/core/cv/metrics/accuracy.hpp>
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>
#include <mlpack/core/data/split_data.hpp>
#include <mlpack/methods/decision_tree/decision_tree.hpp>
#include <mlpack/methods/random_forest/random_forest.hpp>
#include <mlpack/methods/naive_bayes/naive_bayes_classifier.hpp>
#include <mlpack/methods/linear_svm/linear_svm.hpp>
#include <mlpack/methods/pca/pca.hpp>
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
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <random>
#include "CSVProcessor.h"
#include "Operations1.h"
#include "Encoding.h"
#include "StatMeasures.h"
#include "Scalers.h"

using namespace std;
using namespace arma;
using namespace mlpack::cv;
using namespace mlpack::regression;
using namespace mlpack::data;
using namespace mlpack::pca;
using namespace mlpack::tree;
using namespace mlpack::neighbor;
using namespace mlpack::svm;


double accuracy(const Row<size_t>& predictions, const Row<size_t>& labels)
{
	size_t correct = 0;
	for (size_t i = 0; i < predictions.n_elem; ++i)
	{
		if (predictions[i] == labels[i])
			correct++;
	}

	return static_cast<double>(correct) / predictions.n_elem;
}


void ConvertLabelsToSVM(arma::Row<size_t>& labels)
{
	labels.transform([](size_t val) { return val == 0 ? -1 : 1; });
}




int main()
{
	CSVProcessor csvProcessor;
	Scalers scaler;
	StatMeasures stat;
	Encoding encode;

	string filename = "C:\\Users\\Be_Bo\\Downloads\\no_attack_csv.csv";
	vector<vector<string>> noAttack = csvProcessor.readCSV(filename);

	filename = "C:\\Users\\Be_Bo\\Downloads\\attack_csv.csv";
	vector<vector<string>> attack = csvProcessor.readCSV(filename);

	//// count null values in each feature
	cout << "no attack \nNull values for each column:" << endl;
	csvProcessor.countNulls(noAttack);
	cout << endl;


	cout << "attack\nNull values for each column:" << endl;
	csvProcessor.countNulls(attack);
	cout << endl;


 // drop the empty columns
	vector <int> dropcolumns = { 25,24,23,22,18,17,0 };
	for (auto columnindex : dropcolumns)
	{
		csvProcessor.removeColumn(attack, columnindex);
		csvProcessor.removeColumn(noAttack, columnindex);
	}


	cout << "no attack\nNull values for each column after drop empty columns:" << endl;
	csvProcessor.countNulls(noAttack);
	cout << endl;


	cout << "attack\nNull values for each column after drop empty columns:" << endl;
	csvProcessor.countNulls(attack);
	cout << endl;

	// fill null values
	csvProcessor.fillNulls(attack, 5, "47095");
	csvProcessor.fillNulls(attack, 7, "5469");
	csvProcessor.fillNulls(noAttack, 5, "45360.5");
	csvProcessor.fillNulls(noAttack, 7, "53"); 
	cout << endl;


	cout << "no attack \nNull values for each column after replace nulls:" << endl;
	csvProcessor.countNulls(noAttack);
	cout << endl;


	cout << "attack\nNull values for each column after replace nulls:" << endl;
	csvProcessor.countNulls(attack);
	cout << endl;

	// concate the 2 datasets
	vector<vector<string>> concdata = csvProcessor.appendCSV(attack, noAttack);

	dataOperations::count_Duplicated(concdata);
	dataOperations::remove_duplicated(concdata);

	cout << "new shape "; dataOperations::shape(concdata); 
	cout << "\ndata after concate\n";

	csvProcessor.displayCSV(concdata);
	

	vector <string> labels = csvProcessor.extractColumn(concdata, 26);
	vector <double> y = csvProcessor.convertColumnToDouble(labels);

	csvProcessor.removeColumn(concdata, 26); // drop the target

	vector<vector<string>> tempdata = concdata;

	vector < int> encodeInedx = { 10 ,6,4,3,2 };   // 2 3 4 6  10  these columns have string values

	for (auto column : encodeInedx)
	{
		csvProcessor.removeColumn(concdata, column);
	}

	vector<vector<double>> doubleData = csvProcessor.convertToDouble(concdata);

	for (int i = doubleData[0].size(); i >= 0; i--)
	{
		vector<double> scaledVec = scaler.standardScaler(csvProcessor.extractColumn(doubleData, i));
		csvProcessor.removeColumn(doubleData, i);
		csvProcessor.addColumn(doubleData, scaledVec);
	}


	for (auto column : encodeInedx)
	{
		vector<vector<double>> temp = encode.oneHotEncode(csvProcessor.extractColumn(tempdata, column));
		csvProcessor.addColumns(doubleData, temp);

	}
	cout << "our final data \n";
	csvProcessor.displayCSV(doubleData);

	//dataOperations::export_to_csv(doubleData, "data.csv");
	//dataOperations::export_to_csv_vector(y, "labels.csv");



// //======================================= Models MLpack ==================================
//
//
		  arma::mat data;
		  mlpack::data::Load("C:\\Users\\Be_Bo\\source\\repos\\temp\\temp\\data.csv", data, true);

		  // Load the responses 
		  arma::Row<size_t> target;
		  mlpack::data::Load("C:\\Users\\Be_Bo\\source\\repos\\temp\\temp\\labels.csv", target, true);

		  // Split the data into training and test sets
		  arma::mat trainData, testData;
		  arma::Row<size_t> trainResponses, testResponses;
		  const double testRatio = 0.2; // 20% of the data for testing
		  data::Split(data, target, trainData, testData, trainResponses, testResponses, testRatio);

		  LogisticRegression<> lr(trainData, trainResponses, 0.001); // 0.001 is the regularization parameter

		  // Predict on training data
		  arma::Row<size_t> trainPredictions;
		  lr.Classify(trainData, trainPredictions);

		  // Predict on test data
		  arma::Row<size_t> testPredictions; 
		  lr.Classify(testData, testPredictions);

		  // Calculate accuracies
		  double trainAccuracy = accuracy(trainPredictions, conv_to<Row<size_t>>::from(trainResponses));
		  double testAccuracy = accuracy(testPredictions, conv_to<Row<size_t>>::from(testResponses));

		  cout << "Logistic-Training Accuracy: " << trainAccuracy * 100 << "%" << endl;
		  cout << "Logistic-Test Accuracy: " << testAccuracy * 100 << "%" << endl;

//
// //=========================================================================
//
		  DecisionTree<> dt;
		  int numclasses = 2;
		  dt.Train(trainData, trainResponses, numclasses, 5);

		  // Predict on training data
		  arma::Row<size_t> dtTrainPredictions;
		  dt.Classify(trainData, dtTrainPredictions);

		  // Predict on test data
		  arma::Row<size_t> dtTestPredictions;
		  dt.Classify(testData, dtTestPredictions);

		  // Calculate accuracies
		  double dtTrainAccuracy = accuracy(dtTrainPredictions, conv_to<Row<size_t>>::from(trainResponses));
		  double dtTestAccuracy = accuracy(dtTestPredictions, conv_to<Row<size_t>>::from(testResponses));

		  cout << "Decision Tree - Training Accuracy: " << dtTrainAccuracy * 100 << "%" << endl;
		  cout << "Decision Tree - Test Accuracy: " << dtTestAccuracy * 100 << "%" << endl;
//
//
//
// //=========================================================================
//
//
//
//
		  RandomForest<GiniGain, RandomDimensionSelect> rf;
		  size_t numTrees = 10; 
		  rf.Train(trainData, trainResponses, numclasses, numTrees);

		  // Predict on training data
		  arma::Row<size_t> rfTrainPredictions;
		  rf.Classify(trainData, rfTrainPredictions);

		  // Predict on test data
		  arma::Row<size_t> rfTestPredictions;
		  rf.Classify(testData, rfTestPredictions);

		  // Calculate accuracies
		  double rfTrainAccuracy = accuracy(rfTrainPredictions, conv_to<Row<size_t>>::from(trainResponses));
		  double rfTestAccuracy = accuracy(rfTestPredictions, conv_to<Row<size_t>>::from(testResponses));

		  cout << "Random Forest - Training Accuracy: " << rfTrainAccuracy * 100 << "%" << endl;
		  cout << "Random Forest - Test Accuracy: " << rfTestAccuracy * 100 << "%" << endl;
//
//
//
////=========================================================================
//
		  ConvertLabelsToSVM(trainResponses); 
		  ConvertLabelsToSVM(testResponses); 

		  LinearSVM<> svm;
		  svm.Train(trainData, trainResponses); // train  

		  // Predict on training data using SVM
		  arma::Row<size_t> svmTrainPredictions; 
		  svm.Classify(trainData, svmTrainPredictions); 

		  // Predict on test data using SVM 
		  arma::Row<size_t> svmTestPredictions; 
		  svm.Classify(testData, svmTestPredictions); 

		  // Calculate accuracies for SVM
		  double svmTrainAccuracy = accuracy(svmTrainPredictions, trainResponses);
		  double svmTestAccuracy = accuracy(svmTestPredictions, testResponses);

		  cout << "Linear SVM Training Accuracy: " << svmTrainAccuracy * 100 << "%" << endl;
		  cout << "Linear SVM Test Accuracy: " << svmTestAccuracy * 100 << "%" << endl;


}
