#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

double lowest_ever = 1.0;

int maxpoints = 1;
const int maxfiles = 1500;
const int MAXFEATURES = 45;
const int MAXTYPES = 100;
const int MAXBINS = 1024;
const int MAXCLASSIFIERS = 25;
const int RATIO = 5; //  =  #cases/#bins
int BINS;
int MINCASES = 5 * RATIO;

const int INV_SQRT_2PI = 0.39894228;
const int DEFAULT_BANDWIDTH = 1.0;
const int KERNEL_STD_DEV = 2.25;

double weights[MAXFEATURES];
double CNweights[MAXCLASSIFIERS][MAXFEATURES];
double CNselected[MAXCLASSIFIERS][MAXFEATURES];

typedef struct testcase {
	string record_number;
	string filename;
	string class_label;
	string predicted_labels[MAXCLASSIFIERS];
	testcase * next;
};

typedef struct valuestruct {
	double value;
	valuestruct * next;
};

typedef struct histostruct {
	string type;
	//double storedvalues[MAXFEATURES];
	valuestruct * storedvalues[MAXFEATURES];
	double histogram[MAXFEATURES][MAXBINS];
	histostruct * next;
};

const double outliermult = 4.5;
double globalmean[MAXFEATURES];
double globalstddev[MAXFEATURES];

typedef struct typenode {
	string type;
	unsigned count;
	typenode * next;
};

typedef struct filestats {
	string filename;
	float totalsamples;
	float percentmisclass;
	typenode * typelist;
};

const double tailwidth = 0.005;
const int windowsize = 1;

void updatecmatrix(string predicted_label, string class_label,
		int numberoftypes, string typelist[MAXTYPES],
		int confmatrix[MAXTYPES][MAXTYPES]);

void updatehistos(double examplebuffer[MAXFEATURES], string class_label,
		int numberoffields, double ranges[MAXFEATURES],
		double min_values[MAXFEATURES], double max_values[MAXFEATURES],
		histostruct * &histolist);

void updateglobalhistos(double examplebuffer[], int numberoffields,
		double ranges[], double min_values[], double max_values[],
		double globalhistos[MAXFEATURES][MAXBINS]);

string prediction(double test_cases[windowsize][MAXFEATURES],
		histostruct * histolist, int numberoffields, double mins[MAXFEATURES],
		double maxes[MAXFEATURES], double ranges[MAXFEATURES], bool &found);

int main(int argc, char* argv[]) {
	ifstream examplefile(argv[1]);
	ifstream testfile(argv[2]);
	ofstream misclassrecords(argv[3]);
	string linestring, fieldstring;
	string class_labels[windowsize];
	//string test_case_label;
	string classlabel;
	int numberoffields;
	int fieldend;
	int examplecases = 0;
	int numberoftypes = 0;
	int numberofexampletypes = 0;
	double scores[MAXCLASSIFIERS];
	string typelist[MAXTYPES];
	string exampletypelist[MAXTYPES];
	int typecount[MAXTYPES];
	int exampletypecount[MAXTYPES];
	int confmatrix[MAXTYPES][MAXTYPES];
	histostruct * auxhisto;
	histostruct * histolist = 0;
	valuestruct examplevalues[MAXFEATURES];
	int testcases = 0;
	double correct = 0.0;
	double errors = 0.0;
	int i, j;
	for (i = 0; i < MAXTYPES; i++)
		for (j = 0; j < MAXTYPES; j++)
			confmatrix[i][j] = 0;
	bool selected[MAXFEATURES];
	double test_cases[windowsize][MAXFEATURES];
	double min_values[MAXFEATURES];
	double max_values[MAXFEATURES];
	double newmin[MAXFEATURES];
	double newmax[MAXFEATURES];
	double ranges[MAXFEATURES];
	double totalcumulative = 0;
	double examplebuffer[MAXFEATURES];
	string class_label;
	string predicted_label;
	char const *buffer;
	double field_value;
	int numberofclassifiers;
	cout << "\nEnter number of fields:\n";
	cin >> numberoffields;
	cout << "\nEnter number of classifiers:\n";
	cin >> numberofclassifiers;
	int CNnumberselected[MAXCLASSIFIERS];
	double CNweights[MAXCLASSIFIERS][MAXFEATURES];
	int CNselected[MAXCLASSIFIERS][MAXFEATURES];
	int CN;
	int numberselected;
	char response;
	for (CN = 0; CN < numberofclassifiers; CN++) {
		cout << "\nEnter score for classifier " << CN << ":\n";
		cin >> scores[CN];
	}
	for (CN = 0; CN < numberofclassifiers; CN++) {
		numberselected = 0;
		cout << "\nFor classifier " << CN
				<< " and for each field number listed, indicate with 'y' or 'n'\n"
				<< "whether or not to use the field.\n";
		for (i = 0; i < numberoffields; i++) {
			cout << "\nField " << i << " ?:\n";
			cin >> response;
			if (response == 'y') {
				numberselected++;
				CNselected[CN][i] = true;
			} else
				CNselected[CN][i] = false;
		}
		CNnumberselected[CN] = numberselected;
		cout << "\nFor classifier " << CN << ":\n";
		for (i = 0; i < numberselected; i++) {
			cout << "\nEnter weight for field " << i << endl;
			cin >> CNweights[CN][i];
		}
	}
	typedef struct globalh {
		double globalhistos[MAXFEATURES][MAXBINS];
	};

	globalh * gptr = 0;
	histostruct * aux;
	double freqsperfield[MAXFEATURES];

	testcase * testcaselist = 0;
	testcase * auxtestcase;
	bool labelfound;
	getline(testfile, linestring);
	while (!testfile.eof()) {
		auxtestcase = new testcase;
		auxtestcase->next = testcaselist;
		testcaselist = auxtestcase;
		fieldend = linestring.find_first_of(",", 0);
		fieldstring.assign(linestring, 0, fieldend);
		auxtestcase->record_number = fieldstring;
		//erase record number
		linestring.erase(0, fieldend + 1);
		for (i = 0; i < numberoffields; i++) {
			fieldend = linestring.find_first_of(",", 0);
			linestring.erase(0, fieldend + 1);
		}
		//get object name
		fieldend = linestring.find_first_of(",", 0);
		fieldstring.assign(linestring, 0, fieldend);
		class_label = fieldstring;
		auxtestcase->class_label = class_label;
		labelfound = false;
		for (int z = 0; z < numberoftypes; z++)
			if (typelist[z] == class_label) {
				labelfound = true;
				break;
			}
		if (labelfound == false) {
			typecount[numberoftypes] = 0;//actually counted in examplefile
			typelist[numberoftypes++] = class_label;
		}
		getline(testfile, linestring);
	}
	cout << "\nclass_labels in testcaselist:\n";
	auxtestcase = testcaselist;
	while (auxtestcase) {
		cout << auxtestcase->class_label << endl;
		auxtestcase = auxtestcase->next;
	}

	int callstopredict = 0;
	filestats filestatlist[maxfiles];
	int numfiles = 0;
	string filename;
	bool found;
	bool right;
	int findex;
	typenode * tptr;
	for (CN = 0; CN < numberofclassifiers; CN++) //outer CN loop
	{
		numberselected = CNnumberselected[CN];
		for (i = 0; i < numberoffields; i++)
			selected[i] = CNselected[CN][i];
		for (i = 0; i < numberselected; i++)
			weights[i] = CNweights[CN][i];

		gptr = new globalh;
		examplecases = 0;
		for (i = 0; i < numberoffields; i++)
			freqsperfield[i] = 0;

		examplefile.clear();
		examplefile.seekg(0, ios::beg);
		getline(examplefile, linestring);
		while (!examplefile.eof()) {
			j = -1;
			fieldend = linestring.find_first_of(",", 0);
			//erase record number
			linestring.erase(0, fieldend + 1);
			for (i = 0; i < numberoffields; i++) {
				fieldend = linestring.find_first_of(",", 0);
				fieldstring.assign(linestring, 0, fieldend);
				buffer = fieldstring.c_str();
				field_value = atof(buffer);
				if (selected[i])
					examplebuffer[++j] = field_value;
				linestring.erase(0, fieldend + 1);
			}
			for (i = 0; i < numberselected; i++) {
				if (freqsperfield[i] == 0) {
					freqsperfield[i] = 1;
					min_values[i] = examplebuffer[i];
					max_values[i] = examplebuffer[i];
				} else {
					if (examplebuffer[i] < min_values[i])
						min_values[i] = examplebuffer[i];
					if (examplebuffer[i] > max_values[i])
						max_values[i] = examplebuffer[i];
					freqsperfield[i]++;
				}
			}
			getline(examplefile, linestring);
		}

		for (i = 0; i < numberselected; i++)
			ranges[i] = max_values[i] - min_values[i];
		//display max and min and ranges and mean and std dev
		for (i = 0; i < numberselected; i++)
			cout << endl << i << ": " << max_values[i] << " " << min_values[i]
					<< " " << ranges[i];

		//create a global histogram
		BINS = MAXBINS;
		bool labelfound;
		for (j = 0; j < numberselected; j++)
			for (i = 0; i < BINS; i++)
				gptr->globalhistos[j][i] = 0;
		examplefile.clear();
		examplefile.seekg(0, ios::beg);
		getline(examplefile, linestring);
		while (!examplefile.eof()) {
			j = -1;
			fieldend = linestring.find_first_of(",", 0);
			//erase record number
			linestring.erase(0, fieldend + 1);
			for (i = 0; i < numberoffields; i++) {
				fieldend = linestring.find_first_of(",", 0);
				fieldstring.assign(linestring, 0, fieldend);
				buffer = fieldstring.c_str();
				field_value = atof(buffer);
				if (selected[i])
					examplebuffer[++j] = field_value;
				linestring.erase(0, fieldend + 1);
			}
			updateglobalhistos(examplebuffer, numberselected, ranges,
					min_values, max_values, gptr->globalhistos);
			getline(examplefile, linestring);
		}

		int kk;
		double total;
		//convert to relative frequencies
		for (j = 0; j < numberselected; j++) {
			total = 0.0;
			for (kk = 0; kk < BINS; kk++)
				total += gptr->globalhistos[j][kk];
			for (kk = 0; kk < BINS; kk++)
				gptr->globalhistos[j][kk] /= total;
		}
		//display globalhistos
		for (j = 0; j < numberselected; j++) {
			cout << endl;
			for (i = 0; i < BINS; i++)
				cout << gptr->globalhistos[j][i] << " ";
		}
		cout << endl;
		//for each feature, find tailwidth and 1-tailwidth values
		for (j = 0; j < numberselected; j++) {
			i = 0;
			totalcumulative = 0;
			totalcumulative = gptr->globalhistos[j][i];
			while (totalcumulative <= tailwidth)
				totalcumulative += gptr->globalhistos[j][++i];
			//calculate newmin[j] from i,ranges,min_values,max_values
			newmin[j] = (i / (double) BINS) * ranges[j] + min_values[j];
			cout << "\nnewmin[" << j << "] = " << newmin[j];
			while (totalcumulative <= 1 - tailwidth)
				totalcumulative += gptr->globalhistos[j][++i];
			newmax[j] = (i / (double) BINS) * ranges[j] + min_values[j];
			cout << "\nnewmax[" << j << "] = " << newmax[j];
			min_values[j] = newmin[j];
			max_values[j] = newmax[j];
			ranges[j] = max_values[j] - min_values[j];
		}
		delete gptr;

		//create histograms
		while (histolist) {
			auxhisto = histolist;
			histolist = histolist->next;
			delete auxhisto;
		}
		histolist = 0;
		numberofexampletypes = 0;
		examplefile.clear();
		examplefile.seekg(0, ios::beg);
		getline(examplefile, linestring);
		while (!examplefile.eof()) {
			j = -1;
			fieldend = linestring.find_first_of(",", 0);
			//erase record number
			linestring.erase(0, fieldend + 1);
			for (i = 0; i < numberoffields; i++) {
				fieldend = linestring.find_first_of(",", 0);
				fieldstring.assign(linestring, 0, fieldend);
				buffer = fieldstring.c_str();
				field_value = atof(buffer);
				if (selected[i])
					examplebuffer[++j] = field_value;
				linestring.erase(0, fieldend + 1);
			}
			//get object name
			fieldend = linestring.find_first_of(",", 0);
			fieldstring.assign(linestring, 0, fieldend);
			class_label = fieldstring;
			labelfound = false;
			for (int z = 0; z < numberoftypes; z++)
				if (typelist[z] == class_label) {
					labelfound = true;
					typecount[z]++;
					break;
				}
			if (labelfound == false) {
				typelist[numberoftypes] = class_label;
				typecount[numberoftypes++] = 1;
			}
			labelfound = false;
			for (int zz = 0; zz < numberofexampletypes; zz++)
				if (exampletypelist[zz] == class_label) {
					labelfound = true;
					exampletypecount[zz]++;
					break;
				}
			if (labelfound == false) {
				exampletypelist[numberofexampletypes] = class_label;
				exampletypecount[numberofexampletypes++] = 1;
			}
			updatehistos(examplebuffer, class_label, numberselected, ranges,
					min_values, max_values, histolist);
			getline(examplefile, linestring);
		}
		if (maxpoints > numberofexampletypes)
			maxpoints = numberofexampletypes;
		//display typelist, typecount
		cout << endl;
		int k;
		for (k = 0; k < numberoftypes; k++)
			cout << typelist[k] << ", " << typecount[k] << endl;
		//convert to relative frequencies
		int maxcount = typecount[0];
		for (j = 1; j < numberoftypes; j++)
			if (typecount[j] > maxcount)
				maxcount = typecount[j];
		aux = histolist;
		while (aux != 0) {
			for (j = 0; j < numberselected; j++)
				for (k = 0; k < BINS; k++)
					aux -> histogram[j][k] += 1.0 / (maxcount + 1);
			for (j = 0; j < numberselected; j++) {
				total = 0.0;
				for (k = 0; k < BINS; k++)
					total += aux -> histogram[j][k];
				for (k = 0; k < BINS; k++)
					aux -> histogram[j][k] /= total;
			}
			aux = aux -> next;
		}
		//display histograms
		aux = histolist;
		//while (aux != 0)
		//{
		//cout<<"\n\n\nType = "<<aux -> type<<endl;
		//for (j=0; j<numberselected; j++)
		//{
		//cout<<"\nField number "<<j<<" :\n";
		//for (k=0; k<BINS; k++)
		//cout<<aux -> histogram[j][k]<<"  ";
		//}
		//aux = aux -> next;
		//}
		//read test cases, predict
		testfile.clear();
		testfile.seekg(0, ios::beg);
		correct = 0.0;
		errors = 0.0;
		//get first windowsize samples
		for (int ii = 1; ii < windowsize; ii++) {
			getline(testfile, linestring);
			j = -1;
			fieldend = linestring.find_first_of(",", 0);
			//erase record number
			linestring.erase(0, fieldend + 1);
			for (i = 0; i < numberoffields; i++) {
				fieldend = linestring.find_first_of(",", 0);
				fieldstring.assign(linestring, 0, fieldend);
				buffer = fieldstring.c_str();
				field_value = atof(buffer);
				if (selected[i])
					test_cases[ii][++j] = field_value;
				linestring.erase(0, fieldend + 1);
			}
		}
		string savedrecord;
		string recordnumber;
		testcases = 0;
		getline(testfile, linestring);
		while (!testfile.eof()) {
			savedrecord = linestring;
			testcases++;
			//shift down cases
			for (i = 0; i < windowsize - 1; i++)
				for (j = 0; j < numberoffields; j++)
					test_cases[i][j] = test_cases[i + 1][j];
			//parse linestring, assign to test_cases[]
			j = -1;
			fieldend = linestring.find_first_of(",", 0);
			fieldstring.assign(linestring, 0, fieldend);
			recordnumber = fieldstring;
			//cout<<"\nIn testfile loop, recordnumber = "<<recordnumber;
			//erase record number
			linestring.erase(0, fieldend + 1);
			for (i = 0; i < numberoffields; i++) {
				fieldend = linestring.find_first_of(",", 0);
				fieldstring.assign(linestring, 0, fieldend);
				buffer = fieldstring.c_str();
				field_value = atof(buffer);
				if (selected[i])
					test_cases[windowsize - 1][++j] = field_value;
				linestring.erase(0, fieldend + 1);
			}
			//get object name
			fieldend = linestring.find_first_of(",", 0);
			fieldstring.assign(linestring, 0, fieldend);
			class_label = fieldstring;
			callstopredict++;
			//cout<<"\nCorrect type = "<<class_label;
			predicted_label = prediction(test_cases, histolist, numberselected,
					min_values, max_values, ranges, found);
			//if (predicted_label == class_label)
			//cout<<"Was correct\n";
			//else
			//cout<<"Was incorrect\n";
			//find savedrecord in testcaselist, add predicted_label
			//to predicted_labels[]
			//final classification, misclassrecords output after CN loop
			auxtestcase = testcaselist;
			while (auxtestcase->record_number != recordnumber)
				auxtestcase = auxtestcase->next;
			auxtestcase->predicted_labels[CN] = predicted_label;
			getline(testfile, linestring);
		}
	}//end CN loop
	//for each testcase in testcaselist, find predicted_label
	//by majority vote, call updatecmatrix, etc.
	typedef struct votestruct {
		string type;
		double points;
	};
	votestruct votes[MAXTYPES];
	int totaltypes;
	double maxpoints;
	string winningtype;
	string current_type;
	auxtestcase = testcaselist;
	while (auxtestcase) {
		filename = auxtestcase->filename;
		class_label = auxtestcase->class_label;
		totaltypes = 0;
		for (j = 0; j < numberofclassifiers; j++) {
			current_type = auxtestcase->predicted_labels[j];
			found = false;
			for (i = 0; i < totaltypes; i++)
				if (votes[i].type == current_type) {
					//votes[i].points++;
					votes[i].points += scores[j];
					found = true;
					break;
				}
			if (!found) {
				votes[totaltypes].type = current_type;
				//votes[totaltypes++].points = 1;
				votes[totaltypes++].points = scores[j];
			}
		}
		maxpoints = votes[0].points;
		winningtype = votes[0].type;
		for (i = 1; i < totaltypes; i++)
			if (votes[i].points > maxpoints) {
				maxpoints = votes[i].points;
				winningtype = votes[i].type;
			}
		predicted_label = winningtype;
		cout << "\npredicted_label = " << predicted_label << ", class_label = "
				<< class_label;
		updatecmatrix(predicted_label, class_label, numberoftypes, typelist,
				confmatrix);
		if (predicted_label == class_label) {
			correct++;
			cout << "\ncorrect";
			right = true;
		} else {
			errors++;
			right = false;
			misclassrecords << auxtestcase->record_number << "," << class_label
					<< "," << predicted_label << endl;
		}
		found = false;
		for (i = 0; i < numfiles; i++) {
			if (filestatlist[i].filename == filename) {
				found = true;
				findex = i;
				filestatlist[i].totalsamples++;
				if (!right)
					filestatlist[i].percentmisclass++;
				break;
			}
		}
		if (!found)//filename not found
		{
			filestatlist[numfiles].typelist = new typenode;
			filestatlist[numfiles].typelist->next = 0;
			filestatlist[numfiles].typelist->type = predicted_label;
			filestatlist[numfiles].typelist->count = 1;
			filestatlist[numfiles].filename = filename;
			filestatlist[numfiles].totalsamples = 1;
			;
			if (!right)
				filestatlist[numfiles].percentmisclass = 1;
			else
				filestatlist[numfiles].percentmisclass = 0;
			numfiles++;
		} else {
			tptr = filestatlist[findex].typelist;
			found = false;
			while (tptr) {
				if (tptr->type == predicted_label) {
					found = true;
					tptr->count++;
					break;
				}
				tptr = tptr->next;
			}
			if (!found)//this class_label not found for this filename
			{
				tptr = new typenode;
				tptr->type = predicted_label;
				tptr->count = 1;
				tptr->next = filestatlist[findex].typelist;
				filestatlist[findex].typelist = tptr;
			}
		}
		auxtestcase = auxtestcase->next;
	}

	//calculate percentmisclass, sort by percentmisclass
	for (i = 0; i < numfiles; i++) {
		filestatlist[i].percentmisclass = 100.0
				* filestatlist[i].percentmisclass
				/ filestatlist[i].totalsamples;
	}
	int maxindex;
	filestats tempstats;
	for (i = 0; i < numfiles - 1; i++) {
		maxindex = i;
		for (j = i + 1; j < numfiles; j++)
			if (filestatlist[j].percentmisclass
					> filestatlist[maxindex].percentmisclass)
				maxindex = j;
		tempstats = filestatlist[i];
		filestatlist[i] = filestatlist[maxindex];
		filestatlist[maxindex] = tempstats;
	}

	cout << "\nCalls to predict() = " << callstopredict;
	cout
			<< "\n===================================================================";
	cout << "\nTotal number of testcases = " << testcases;
	cout << "\nNumber of correct classifications = " << correct;
	cout << "\nNumber of incorrect classifications = " << errors;
	cout << "\nNumber unclassified = " << testcases - correct - errors;
	cout << "\nPercent correct = " << 100.0 * correct / testcases;
	//cout<<"\nSmallest normalized max. likelihood value = "
	//<<lowest_ever<<endl;
	cout << "\n\nClassification matrix:\n\n";
	for (i = 1; i <= numberoftypes; i++)
		cout << setiosflags(ios::right) << setw(3) << i << "  ";
	cout << "    <-classified as\n";
	for (i = 1; i <= numberoftypes; i++)
		cout << "---  ";
	cout << endl;
	for (i = 0; i < numberoftypes; i++) {
		for (j = 0; j < numberoftypes; j++)
			cout << setiosflags(ios::right) << setw(3) << confmatrix[j][i]
					<< "  ";
		cout << "    (" << i + 1 << "): class " << typelist[i] << endl;
	}

	cout << "\nClassifications by type:";
	for (i = 0; i < numberoftypes; i++) {
		cout << endl << endl << "-------" << typelist[i] << "-------";
		for (j = 0; j < numberoftypes; j++)
			if (confmatrix[j][i] >= 1) {
				cout << endl << typelist[j] << ": " << confmatrix[j][i];
				if (j == i)
					cout << " <---";
			}
	}

	//cout<<"\n\nTotal number of files = "<<numfiles;
	//cout<<"\nFile statistics:\n";
	//for (i=0; i<numfiles; i++)
	//{
	//cout<<"\n------------------------------------------------------------------\n";
	//cout<<filestatlist[i].filename<<", percent incorrect = "<<
	//filestatlist[i].percentmisclass<<", total samples = "<<
	//filestatlist[i].totalsamples<<endl;
	//cout<<"Types classified as:\n";
	//tptr=filestatlist[i].typelist;
	//while (tptr)
	//{
	//cout<<tptr->type<<", "<<tptr->count<<endl;
	//tptr=tptr->next;
	//}
	//}


}

float normalKernel(float x) {
	return (INV_SQRT_2PI / KERNEL_STD_DEV) * exp(-pow(x, 2) / (2.0 * pow(
			KERNEL_STD_DEV, 2)));
}

float getDensityEstimation(float x, float n, valuestruct *data) {
	int i;
	float density = 0;
	valuestruct * v = data;
	do {
		density += normalKernel((x - v -> value) / DEFAULT_BANDWIDTH);
	} while (data -> next);
	//for(i = 0; i < n; ++i)
	//	density += normalKernel((x - data[i]) / DEFAULT_BANDWIDTH);
	return density / (n * DEFAULT_BANDWIDTH);
}

string prediction(double test_cases[windowsize][MAXFEATURES],
		histostruct * histolist, int numberoffields, double mins[MAXFEATURES],
		double maxes[MAXFEATURES], double ranges[MAXFEATURES], bool &found) {
	typedef struct typeinfo {
		string name;
		double probability;
		typeinfo * next;
	};
	typeinfo * labels[windowsize];
	typeinfo * typeptr;
	typeinfo * before;
	typeinfo * after;
	double threshold = 0.02; //for no match
	double max_likelihood;
	double current_likelihood;
	double current_likelihood_KDE;
	double proportion;
	int i, j;
	int index;
	int maxvoteindex;
	found = true;
	string best_type;
	histostruct * aux;
	double likelihoodtotal;
	double likelihoodtotal_KDE;
	for (i = 0; i < windowsize; i++) {
		labels[i] = 0;
		aux = histolist;
		likelihoodtotal = 0.0;
		likelihoodtotal_KDE = 0.0;
		while (aux != 0) {
			current_likelihood = 1.0;
			current_likelihood_KDE = 1.0;
			for (j = 0; j < numberoffields; j++) {
				if (!weights[j])
					continue;
				if (test_cases[i][j] > maxes[j])
					proportion = 1.0;
				else if (test_cases[i][j] < mins[j])
					proportion = 0.0;
				else if (maxes[j] == mins[j])
					proportion = 1.0;
				else
					proportion = (test_cases[i][j] - mins[j]) / ranges[j];
				if (proportion >= (1.0 * BINS - 1.0) / BINS)
					index = BINS - 1;
				else
					index = ((int) (proportion * BINS)) % BINS;
				//current_likelihood *= aux -> histogram[j][index];
				current_likelihood *= pow(aux -> histogram[j][index],
						weights[j]);
				current_likelihood_KDE *= getDensityEstimation(
						test_cases[i][j], numberoffields,
						aux -> storedvalues[j]);
			}
			likelihoodtotal += current_likelihood;
			likelihoodtotal_KDE += current_likelihood_KDE;
			aux = aux -> next;
		}
		if (likelihoodtotal == 0.0) {
			if (likelihoodtotal_KDE != 0.0) {
				cout << "KDE found, histo didn't";
			} else {
				found = false;
				return best_type;
			}
		}
		aux = histolist;
		while (aux != 0) {
			current_likelihood = 1.0;
			current_likelihood_KDE = 1.0;
			for (j = 0; j < numberoffields; j++) {
				if (!weights[j])
					continue;
				if (test_cases[i][j] > maxes[j])
					proportion = 1.0;
				else if (test_cases[i][j] < mins[j])
					proportion = 0.0;
				else if (maxes[j] == mins[j])
					proportion = 1.0;
				else
					proportion = (test_cases[i][j] - mins[j]) / ranges[j];
				if (proportion >= (1.0 * BINS - 1.0) / BINS)
					index = BINS - 1;
				else
					index = ((int) (proportion * BINS)) % BINS;
				//current_likelihood *= aux -> histogram[j][index];
				current_likelihood *= pow(aux -> histogram[j][index],
						weights[j]);
				current_likelihood_KDE *= getDensityEstimation(
						test_cases[i][j], numberoffields,
						aux -> storedvalues[j]);
			}
			current_likelihood /= likelihoodtotal;
			current_likelihood_KDE /= likelihoodtotal_KDE;
			//insert new type in descending order by probability
			typeptr = new typeinfo;
			typeptr -> name = aux -> type;
			typeptr -> probability = current_likelihood;
			if (labels[i] == 0) {
				typeptr -> next = 0;
				labels[i] = typeptr;
			} else {
				before = labels[i];
				after = 0;
				while (before && (before->probability > current_likelihood)) {
					after = before;
					before = before->next;
				}
				if (!after) {
					typeptr->next = before;
					labels[i] = typeptr;
				} else {
					typeptr->next = before;
					after->next = typeptr;
				}
			}
			aux = aux -> next;
		}
	}
	//cout<<"\nWindow contents:\n";
	//for (i=0; i<windowsize; i++)
	//{
	//cout<<i<<":"<<endl;
	//typeptr = labels[i];
	//int countitems = 0;
	//while (typeptr && countitems < maxpoints)
	//while (typeptr && countitems < 1)
	//{
	//cout<<typeptr->name<<", "<<typeptr->probability<<endl;
	//typeptr=typeptr->next;
	//countitems++;
	//}
	//}

	typedef struct votestruct {
		string type;
		double points;
	};
	votestruct votes[MAXTYPES];
	bool typefound;
	int typesfound = 0;
	int points, pos;
	string currentlabel;
	for (i = 0; i < windowsize; i++) {
		typeptr = labels[i];
		for (pos = 0; pos < maxpoints; pos++) {
			currentlabel = typeptr->name;
			typefound = 0;
			for (j = 0; j < typesfound; j++)
				if (votes[j].type == currentlabel) {
					typefound = 1;
					votes[j].points += maxpoints - pos;
					break;
				}
			if (!typefound) {
				votes[typesfound].type = currentlabel;
				votes[typesfound].points = maxpoints - pos;
				typesfound++;
			}
			//for (int jj=0; jj<typesfound; jj++)
			//cout<<"\nvotes["<<jj<<"].type = "<<votes[jj].type
			//<<", points = "<<votes[jj].points;
			typeptr = typeptr -> next;
		}
	}

	maxvoteindex = 0;
	for (i = 1; i < typesfound; i++)
		if (votes[i].points > votes[maxvoteindex].points)
			maxvoteindex = i;
	best_type = votes[maxvoteindex].type;
	//cout<<"winner = "<<best_type<<endl;
	return best_type;
}

void updatehistos(double examplebuffer[], string class_label,
		int numberoffields, double ranges[], double min_values[],
		double max_values[], histostruct * &histolist) {
	//cout<<"\nIncoming class_label: "<<class_label;
	//cout<<"\nfield values: \n";
	//for (int z=0; z<numberoffields; z++)
	//cout<<examplebuffer[z]<<", ";
	int i, j, k;
	double proportion;
	int index;
	histostruct * aux;
	bool found = false;
	aux = histolist;
	while (aux != 0) {
		if (class_label == aux -> type) {
			found = true;
			break;
		}
		aux = aux -> next;
	}
	if (!found) {
		aux = new histostruct;
		aux -> type = class_label;
		aux -> next = histolist;
		histolist = aux;
		for (j = 0; j < numberoffields; j++)
			for (k = 0; k < BINS; k++)
				//aux -> histogram[j][k] = 0.125;
				//aux -> histogram[j][k] = 0.03125;
				aux -> histogram[j][k] = 0;
	}
	//else
	for (j = 0; j < numberoffields; j++) {
		if (ranges[j] == 0.0)
			proportion = 1.0;
		else
			proportion = (examplebuffer[j] - min_values[j]) / ranges[j];
		if (proportion >= (1.0 * BINS - 1.0) / BINS)
			index = BINS - 1;
		else if (proportion < 0)
			index = 0;
		else
			index = ((int) (proportion * BINS)) % BINS;
		aux -> histogram[j][index]++;
		valuestruct *v = aux -> storedvalues[k];
		if (!v)
			v = new valuestruct;
		else {
			while (v -> next)
				v = v -> next;
			v -> next = new valuestruct;
			v = v -> next;
		}
		v -> value = examplebuffer[j];
	}
}

void updatecmatrix(string predicted_label, string class_label,
		int numberoftypes, string typelist[MAXTYPES],
		int confmatrix[MAXTYPES][MAXTYPES]) {
	int row, col;
	row = 0;
	while (predicted_label != typelist[row])
		row++;
	col = 0;
	while (class_label != typelist[col])
		col++;
	confmatrix[row][col]++;
}

void updateglobalhistos(double examplebuffer[], int numberoffields,
		double ranges[], double min_values[], double max_values[],
		double globalhistos[MAXFEATURES][MAXBINS]) {
	int i, j, k;
	double proportion;
	int index;
	for (j = 0; j < numberoffields; j++) {
		if (ranges[j] == 0.0)
			proportion = 1.0;
		else
			proportion = (examplebuffer[j] - min_values[j]) / ranges[j];
		if (proportion >= (1.0 * BINS - 1.0) / BINS)
			index = BINS - 1;
		else
			index = ((int) (proportion * BINS)) % BINS;
		globalhistos[j][index]++;
	}
}

