#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define RAND_MAX 10

//#define PI 3.14159265 
#define INV_SQRT_2PI 0.39894228
#define DEFAULT_BANDWIDTH 1.0

using namespace std;

// e^((x * x) / 2)
float standardNormalDensity(float x) {
	return exp((x * x) / -2.0) * INV_SQRT_2PI;
}

int main(int argc, char **argv) {
	fstream trainingFile(argv[1]);	
	fstream testFile(argv[2]);	
	string linestring;
	int i;
	float *trainingData;
	float *testData;
	float sum;
	float mean;
	float variance;
	float stddev;
	float bandwidth = DEFAULT_BANDWIDTH;

	// read in training data
	//trainingFile.clear();
	//trainingFile.seekg(0, ios::beg);
	//getline(trainingFile, linestring);
	
	int N = 6;
	trainingData = (float *)malloc(N * sizeof(float));
	trainingData[0] = -2.1;
	trainingData[1] = -1.3;
	trainingData[2] = 0.4; 
	trainingData[3] = 1.9;
	trainingData[4] = 5.1;
	trainingData[5] = 6.2;

	// randomly generate test data
/*	data = (float *)malloc(sizeof(float) * N);
	srand(time(NULL));
	for(i = 0; i < N; ++i)
		data[i] = rand() % RAND_MAX;*/

	// calc mean
	sum = 0;
	for(i = 0; i < N; ++i)
		sum += trainingData[i];
	mean = sum / N;

	cout << "Training Data Statistics:\n";
	cout << "Mean: " << mean << "\n";

	// calc variance, stddev
	float sum2 = 0;
	for(i = 0; i < N; ++i)
		sum2 += pow(mean - trainingData[i], 2);
	variance = sum2 / N;
	stddev = sqrt(variance);

	cout << "Variance: " << variance << "\n";
	cout << "Standard Deviation: " << stddev << "\n\n";

	// test data using kernel density estimation
	float test = 2.1;
	float density = 0;
	for(i = 0; i < N; ++i) {
		density += standardNormalDensity((test - trainingData[i]) / bandwidth);
	}
	density /= N * bandwidth;

	cout << "Test value: " << test << "\n";
	cout << "Estimated Density: " << density << "\n";

	return 0;
}
