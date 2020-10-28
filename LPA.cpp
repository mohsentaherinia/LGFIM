#include "stdafx.h"


int NSlpa::findMaximumLabelorg(vector <int>arr){
	// Sort the array 
	int n = arr.size();
	if (n == 0)
		return -1;
	std::sort(arr.begin(), arr.end());

	// find the max frequency using linear traversal 
	int max_count = 1, res = arr[0], curr_count = 1;
	for (int i = 1; i < n; i++) {
		if (arr[i] == arr[i - 1])
			curr_count++;
		else {
			if (curr_count > max_count) {
				max_count = curr_count;
				res = arr[i - 1];
			}
			curr_count = 1;
		}
	}

	// If last element is most frequent 
	if (curr_count > max_count)
	{
		max_count = curr_count;
		res = arr[n - 1];
	}

	if (res==-1) {
		res = arr[n - 1];
	}

	return res;
}

int NSlpa::findMaximumLabel(vector <int>arr, vector <int>vertexCoverSet) {
	vector<int> maxSet;
	// Sort the array 
	int n = arr.size();
	if (n == 0)
		return -1;
	std::sort(arr.begin(), arr.end());

	// find the max frequency using linear traversal 
	int max_count = 1, curr_count = 1;
	//maxSet.push_back(arr[0]);
	for (int i = 1; i < n; i++) {
		if (arr[i] == arr[i - 1])
			curr_count++;
		else {
			if (curr_count > max_count) {
				max_count = curr_count;
				maxSet.clear();
				maxSet.push_back(arr[i - 1]);
			}
			else if (curr_count == max_count) {
				max_count = curr_count;
				maxSet.push_back(arr[i - 1]);
			}
			curr_count = 1;
		}
	}

	// If last element is most frequent 
	if (curr_count > max_count)
	{
		max_count = curr_count;
		maxSet.clear();
		maxSet.push_back(arr[n - 1]);
	}
	else if (curr_count == max_count) {
		max_count = curr_count;
		maxSet.push_back(arr[n - 1]);
	}

	if (maxSet.size() == 0) {
		maxSet.push_back(arr[n - 1]);
	}

	if (maxSet.size() == 1)
		return maxSet[0];
	else {
		for (vector<int>::iterator maxi = maxSet.begin(); maxi != maxSet.end(); ++maxi) {
			for (vector<int>::iterator u = vertexCoverSet.begin(); u != vertexCoverSet.end(); ++u) {
				if (*maxi == *u)
					return *maxi;
			}
		}
		srand(time(0));
		int rnd = rand() % maxSet.size();
		return maxSet[rnd];
	}
}

void NSlpa::updateLabelTable(TIntIntH& org,const TIntIntH&  temp) {
	for (size_t i = 0; i < temp.Len(); i++)
	{
		int _key = temp.GetKey(i);
		org.AddDat(_key) = temp.GetDat(_key);
	}
}

void NSlpa::printCommunityinfo(const TCnComV&CmtyV) {
	cout << "\n\n# comuunity = " << CmtyV.Len();
	for (int c = 0; c < CmtyV.Len(); c++) {
		/*if (CmtyV[c].Len() <= 50)
			continue;
		*/
		printf("\nCmty(%d) ,  #nodes= %d =>", c, CmtyV[c].Len());
		for (int i = 0; i < CmtyV[c].Len(); i++) {
			cout << CmtyV[c][i] << ", ";
		}
		cout << endl;
	}
}
