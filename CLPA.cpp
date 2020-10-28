#include "stdafx.h"


int NSclpa::CummunityCapacity(const int& t, const  int & N) {
	return (int((K_coeff *t) / ConvergeRep) + 1)*N / K_coeff;
}

int NSclpa::findMaximumLabelRandom(TIntIntH FvL, vector <int> A) {
	//Error With x64 bits Debuging
	//error accued randomly when line 4 is executing   -> FvL.SortByDat(false);
	//findMaximumLabelRandom is defined to solve this problem with x64bits Debugger
	int candidate;
	bool flag = true;
	vector <int> S;
	FvL.SortByDat(false);
	
	for (size_t i = 0; i < FvL.Len() && flag == true; i++) {
		for (vector<int>::iterator it = A.begin(); it != A.end(); ++it) {
			int label = FvL.GetKey(i);
			if (label == *it) {
				S.push_back(label);
				break;
			}
		}
		if (i < FvL.Len() - 1)
			if (FvL.GetDat(FvL.GetKey(i)) != FvL.GetDat(FvL.GetKey(i + 1)) && S.size() > 0)
				flag = false;
	}
	
	if (S.size() == 0)
		return -1;
	//srand(time(0));
	int rnd = rand() % S.size();
	return S[rnd];
}

int NSclpa::findMaximumLabelRandom2(TIntIntH &FvL, vector <int> &A) {
	int maxCounter=-1;
	vector <int> CandiateList;
	for (vector<int>::iterator it = A.begin(); it != A.end(); ++it) {
		if (FvL.GetDat(*it) > maxCounter) {
			CandiateList.clear();
			CandiateList.push_back(*it);
		}
		if (FvL.GetDat(*it) == maxCounter) {
			CandiateList.push_back(*it);
		}
	}
	if (CandiateList.size() == 0)
		return -1;
	srand(time(0));
	int rnd = rand() % CandiateList.size();
	return CandiateList[rnd];
}