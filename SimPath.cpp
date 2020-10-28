#include "stdafx.h"
#include <algorithm>

vector<int> NSsimPath::callSimPath1(const PNGraph & Graph, const GlobalConst & GC ,const char* Model , const int& SeedSize, const double& ICProbb, const int& MCS, const double& R, const double& L,vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedSP;// (SeedSize, 0);
	NSsimPath::GetSimPath1(Graph,GC, seedSP, SeedSize, R, L, curveInfoISV, curveInfoRTV, indexPlot,Model,ICProbb,MCS);
	cout << "\nseedSet SimPath1\n";
	NStools::PrintVector(seedSP);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("SimPath1", GC, seedSP);
	return seedSP;
}

vector<int> NSsimPath::callSimPath1_v2(const vector<int> vertexSet, const PNGraph & Graph, const GlobalConst & GC ,const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& R, const double& L) {
	vector<int> seedSP;// (SeedSize, 0);
	NSsimPath::GetSimPath1_v2(vertexSet,Graph, GC, seedSP, SeedSize, R, L, Model, ICProbb, MCS);
	return seedSP;
}

vector<int> NSsimPath::callSimPath2(const PNGraph & Graph, const GlobalConst & GC ,const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& R, const double& L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedSP;// (SeedSize, 0);
	NSsimPath::GetSimPath2(Graph, GC, seedSP, SeedSize, R, L, curveInfoISV, curveInfoRTV, indexPlot, Model, ICProbb, MCS);
	cout << "\nseedSet SimPath2\n";
	NStools::PrintVector(seedSP);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("SimPath2", GC, seedSP);
	return seedSP;
}

vector<int> NSsimPath::callSimPath2_v2(const vector<int> vertexSet, const PNGraph & Graph, const GlobalConst & GC ,const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& R, const double& L) {
	vector<int> seedSP;// (SeedSize, 0);
	NSsimPath::GetSimPath2_v2(vertexSet, Graph, GC, seedSP, SeedSize, R, L, Model, ICProbb, MCS);
	return seedSP;
}

vector<int> NSsimPath::callSimPath3(const PNGraph & Graph, const GlobalConst & GC ,const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& R, const double& L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedSP;// (SeedSize, 0);
	NSsimPath::GetSimPath3(Graph, GC, seedSP, SeedSize, R, L, curveInfoISV, curveInfoRTV, indexPlot, Model, ICProbb, MCS);
	cout << "\nseedSet SimPath3\n";
	NStools::PrintVector(seedSP);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("SimPath3", GC, seedSP);
	return seedSP;
}

vector<int> NSsimPath::callSimPath3_v2(const vector<int> vertexSet, const PNGraph & Graph, const GlobalConst & GC ,const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& R, const double& L) {
	vector<int> seedSP;// (SeedSize, 0);
	NSsimPath::GetSimPath3_v2(vertexSet, Graph, GC, seedSP, SeedSize, R, L, Model, ICProbb, MCS);
	return seedSP;
}

void NSsimPath::DifferenceVector(vector<int>& S1, vector<int>& S2, vector<int>& S3) {
	S3.resize(S1.size());
	vector<int>::iterator it;
	std::sort(S1.begin(), S1.end());
	std::sort(S2.begin(), S2.end());

	it = std::set_difference(S1.begin(), S1.end(), S2.begin(), S2.end(), S3.begin());
	S3.resize(it - S3.begin()); 
}

void NSsimPath::SetIntersection(vector<int>& S1, vector<int>& S2, vector<int>& S3) {
	S3.resize(S1.size());
	vector<int>::iterator it;
	std::sort(S1.begin(), S1.end());
	std::sort(S2.begin(), S2.end());

	it = std::set_intersection(S1.begin(), S1.end(), S2.begin(), S2.end(), S3.begin());
	S3.resize(it - S3.begin()); 
}

void NSsimPath::QueueSort(queue<TIntFltH>& q)
{
	TIntFltH tempH;
	const int  N = q.size();
	tempH.Gen(N);
	TInt _key;
	TFlt _data;
	int i = 0;
	while (!q.empty()) {
		TIntFltH x;
		x = q.front();
		x.GetKeyDat(0, _key, _data);
		tempH.AddDat(_key) = _data;
		q.pop();
		i++;
	}
	////////////////////////////
	tempH.SortByDat(false);
	///////////////////////////
	for (size_t i = 0; i < N; i++)
	{
		tempH.GetKeyDat(i, _key, _data);
		TIntFltH x;
		x.AddDat(_key) = _data;
		q.push(x);
	}
}

void NSsimPath::printqueue(queue<TIntFltH> q) {
	cout << "\n\n\n";
	int i = 1;
	while (!q.empty()) {
		TIntFltH item2 = q.front();
		cout << "\nnode " << item2.GetKey(0) << "\t\tinfectedNodes = " << item2.GetDat(item2.GetKey(0))<<"\t"<<i++;
		q.pop();
	}
}

bool NSsimPath::SatisfyCondition(const int &x, const int &y, const stack<int>&Q, list<int> * D, const vector<int>& Wset) {
	//if (y NOT exist in W) return false;
	bool IsExist1 = false;
	for (size_t i = 0; i < Wset.size(); i++)
		if (Wset[i] == y)
			IsExist1 = true;
	if (IsExist1 == false)
		return false;
	//if (y exist in Q) return false;
	stack<int> tempQ = Q;
	while (!tempQ.empty()) {
		if (y == tempQ.top())
			return false;
		tempQ.pop();
	}
	//if (y exist in D[x]) return false;
	list <int> ::iterator it;
	for (it = D[x].begin(); it != D[x].end(); ++it)
		if (*it == y)
			return false;
	return true;
}

void NSsimPath::PeekUPTopL(vector<int>& Uset, queue<TIntFltH> qCELF, const int &L) {
	int i = 0;
	while (i<L && !qCELF.empty()) {
		TIntFltH x = qCELF.front();
		Uset.push_back(x.GetKey(0));
		i++;
		qCELF.pop();
	}
}

void NSsimPath::RemoveXfromQCELF(queue<TIntFltH> &qCELF, const int &item) {
	queue<TIntFltH> tempQ;
	while (!qCELF.empty()) {
		TIntFltH x = qCELF.front();
		qCELF.pop();
		if (x.GetKey(0) == item)
			continue;
		tempQ.push(x);
		}
	//qCELF = tempQ;
	while (!tempQ.empty()) {
		TIntFltH x = tempQ.front();
		qCELF.push(x);
		tempQ.pop();
	}
}

void NSsimPath::qCELFUpdate(queue<TIntFltH> &qCELF, const int &item , const float & value) {
	queue<TIntFltH> tempQ;
	while (!qCELF.empty()) {
		TIntFltH x = qCELF.front();
		qCELF.pop();
		if (x.GetKey(0) == item)
			x.AddDat(item) = value;
		tempQ.push(x);
	}
	//qCELF = tempQ;
	while (!tempQ.empty()) {
		TIntFltH x = tempQ.front();
		qCELF.push(x);
		tempQ.pop();
	}
}