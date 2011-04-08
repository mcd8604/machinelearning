#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

using namespace std;

int main(int argc, char **argv)
{
	if(argc < 2)
		return 1;

	int i;
	double trainPercent, testPercent;
	string linestring, classLabel, inStr;
	map<string, vector<string> > classes;
	map<string, vector<string> >::iterator it;

	ifstream inFile(argv[1]);
	while (!inFile.eof())
	{
		getline(inFile, linestring);
		if(linestring.length() > 0)
		{
			i = linestring.find_last_of(',');
			classLabel = linestring.substr(i + 1);
			classes[classLabel].push_back(linestring);
		}
	}
	inFile.close();

	cout << "Enter percentage of records for training data: ";
	cin >> inStr;
	trainPercent = atof(inStr.c_str()) / 100.0;

	cout << "Enter percentage of remaining records for testing data: ";
	cin >> inStr;
	testPercent = atof(inStr.c_str()) / 100.0;

	char fileName[64];
     	strcpy(fileName, argv[1]);
	ofstream trainFile(strcat(fileName, "Train"));
     	strcpy(fileName, argv[1]);
	ofstream testFile(strcat(fileName, "Test"));
     	strcpy(fileName, argv[1]);
	ofstream validateFile(strcat(fileName, "Validate"));
	for(it = classes.begin(); it != classes.end(); it++)
	{
		vector<string> classRecords = it->second;
		int n = classRecords.size();	
		int trainNum = n * trainPercent;
		int testNum = (n - trainNum) * testPercent;
		for(i = 0; i < trainNum; ++i)
			trainFile << classRecords[i] + "\n";
		for(; i < trainNum + testNum; ++i)
			testFile << classRecords[i] + "\n";
		for(; i < n; ++i)
			validateFile << classRecords[i] + "\n";
	}	
	trainFile.close();
	testFile.close();
	validateFile.close();

	return 0;
}
