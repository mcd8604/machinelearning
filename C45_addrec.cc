#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
	if(argc < 2)
		return 1;

	ifstream inFile(argv[1]);
	ofstream outFile(strcat(argv[1], "_rec"));
	string linestring;
	int i = 0;
	while (!inFile.eof())
	{
		getline(inFile, linestring);
		if(linestring.length() > 0)
			outFile << ++i << "," + linestring + "\n";
	}

	inFile.close();
	outFile.close();

	return 0;
}
