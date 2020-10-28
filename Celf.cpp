#include "stdafx.h"

vector<int>  NScelf::callCELFGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedCELF;
	NScelf::CELFGlobal(Graph, GC, Model, seedCELF, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("CELFGlobal", GC, seedCELF);
	return seedCELF;
}

vector<int>  NScelf::callCELFGlobalOfLocal(const vector<int>& candidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedCELF;
	NScelf::CELFGlobalOfLocal(candidateSet, Graph, GC, Model, seedCELF, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("CELFGL", GC, seedCELF);
	return seedCELF;
}

vector<int>  NScelf::callCELFLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC) {
	vector<int> seedCELF;
	NScelf::CELFLocal(Graph, GC, Model, seedCELF, SeedSize, ICProb, MC);
	NStools::SaveSeedSetToFileLocal("CELFLocal", GC, seedCELF);
	return seedCELF;
}

vector<int>  NScelf::callCELFPPPGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedCELFPPP;
	NScelf::CELFPPPGlobal(Graph, GC, Model, seedCELFPPP, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("CELFPPPGlobal", GC, seedCELFPPP);
	return seedCELFPPP;
}

vector<int>  NScelf::callCELFPPPGlobalOfLocal(const vector<int>& candidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedCELFPPP;
	NScelf::CELFPPPGlobalOfLocal(candidateSet, Graph, GC, Model, seedCELFPPP, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("CELFPPPGL", GC, seedCELFPPP);
	return seedCELFPPP;
}

vector<int>  NScelf::callCELFPPPLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC) {
	vector<int> seedCELFPPP;
	NScelf::CELFPPPLocal(Graph, GC, Model, seedCELFPPP, SeedSize, ICProb, MC);
	NStools::SaveSeedSetToFileLocal("CELFPPPLocal", GC, seedCELFPPP);
	return seedCELFPPP;
}

void NScelf::QueueSortFirstK(queue<TIntIntH>& q, const int & N)
{
	TIntIntH tempH;
	tempH.Gen(N);
	TInt _key,_data;
	int i = 0;
	while (!q.empty()) {
		TIntIntH x;
		x = q.front();
		x.GetKeyDat(0, _key, _data);
		tempH.AddDat(_key)= _data;
		q.pop();
		i++;
	}
	////////////////////////////
	tempH.SortByDat(false);
	///////////////////////////
	for (size_t i = 0; i < N; i++)
	{
		tempH.GetKeyDat(i, _key, _data);
		TIntIntH x;
		x.AddDat(_key) = _data;
		q.push(x);
	}
}

void NScelf::QueueSortOthersK(queue<TIntIntH>& q, const int & N)
{
	//printqueue(q);
	vector <TIntIntH> temp;
	temp.push_back(q.back());
	while (!q.empty()) {
		temp.push_back(q.front());
		q.pop();
		}
	temp.pop_back();
	//cout << "\nvector\n";
	/*for (size_t it = 0; it < N; it++) {
		TIntIntH x = temp[it];
		int key = x.GetKey(0);
		cout << key << " " << x.GetDat(key) << endl;
	}
*/
	
	////////////////////////////
	//if ((temp.size() == 0) || (temp.size() == 1))
	if ((temp.size() == 0))
		return;
	for (size_t it = 0; it < N; it++)
	{
		if ((it + 1) == N)
			break;
		/*TIntIntH x = temp[it];
		int keyx = x.GetKey(0);
		TIntIntH y = temp[it+1];
		int keyy = y.GetKey(0);
		if (x.GetDat(keyx) >= y.GetDat(keyy))
		*/
		if (temp[it].GetDat(temp[it].GetKey(0)) >= temp[it+1].GetDat(temp[it+1].GetKey(0)))
			break;
		else {
			TIntIntH t = temp[it];
			temp[it] = temp[it + 1];
			temp[it + 1] = t;
		}
	}
		
	///////////////////////////
	for (size_t i = 0; i < N; i++)
		q.push(temp[i]);
	//printqueue(q);
}

void NScelf::printqueue(queue<TIntIntH> q) {
	cout << "\n\n\n";
	while (!q.empty()) {
		TIntIntH item2 = q.front();
		cout << "\nnode " << item2.GetKey(0) << "\t\tinfectedNodes = " << item2.GetDat(item2.GetKey(0));
		q.pop();
	}
}