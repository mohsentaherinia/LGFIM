#include "stdafx.h"

void NScommunity::SaveCommunityInfoTofile(const TCnComV &CmtyV, const string& DataSetName, const string& CmtyAlg) {
	char  Filename[80] = "Z_CD_";
	strcat(Filename, DataSetName.c_str());
	strcat(Filename, "_");
	strcat(Filename, CmtyAlg.c_str());
	strcat(Filename, ".cmt");
	ofstream outputFile;
	outputFile.open(Filename, ios::out);
	int x = 12;
	if (outputFile.is_open())
	{
		outputFile << CmtyV.Len() << endl;
		for (int c = 0; c < CmtyV.Len(); c++)
		{
			for (int i = 0; i < CmtyV[c].Len(); i++)
				outputFile << CmtyV[c][i].Val << ",";
			outputFile << endl;
		}
		outputFile.close();
	}
	else
		cout << "Unable to open file";
}


vector<int>* NScommunity::readCommunityInfoFromFile(const char*DataSetName, const char* CmtyAlg, int & cmtyCount) {
	char  Filename[80] = "Z_CD_";
	strcat(Filename, DataSetName);
	strcat(Filename, "_");
	strcat(Filename, CmtyAlg);
	strcat(Filename, ".cmt");

	ifstream inputFile(Filename, ios::in);
	string line;
	int s = 0;
	int lineNumber = 1;
	vector<int> *cmtyVec = new vector<int>[2];
	if (inputFile.is_open())
	{
		while (getline(inputFile, line, '\n')) {
			if (line.empty())
				continue;
			//process line assuming it is read as a string
			for (string::iterator it = line.begin(); it != line.end(); ++it)
			{
				if (*it != ',') {
					int x = (int)*it - 48;
					s = s * 10 + x;
				}
				else {
					cmtyVec[lineNumber - 2].push_back(s);
					s = 0;
				}
			}
			if (lineNumber == 1) {
				//cout << s ;
				delete[]cmtyVec;
				cmtyVec = new vector<int>[s];
				cmtyCount = s;
				s = 0;
			}
			lineNumber++;
		}
		//PrintVector(ver);
		inputFile.close();
	}
	else
		cout << "Unable to open the my Community Vertex file";
	return cmtyVec;
}


void NScommunity::readLargeCommunityNodesFromFile(TIntV &NIdV) {
	ifstream inputFile("zLargeCmt.txt", ios::in);
	string line;
	int s = 0;
	if (inputFile.is_open())
	{
		while (getline(inputFile, line, '\n')) {
			if (line.empty())
				continue;
			//process line assuming it is read as a string
			for (string::iterator it = line.begin(); it != line.end(); ++it)
			{
				int x = (int)*it - 48;
				s = s * 10 + x;

			}
			NIdV.Add(s);
			s = 0;
		}
		inputFile.close();
	}
	else {
		cout << "\nUnable to open the my Community Vertex file\n";
		system("pause");
	}
}
