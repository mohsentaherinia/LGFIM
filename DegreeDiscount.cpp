#include "stdafx.h"
vector<int>  NSDegreeDiscount::callDegDisGlobal(const PNGraph &Graph, const GlobalConst &GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC,const double& ratio, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedDegDis;
	NSDegreeDiscount::DegDisGlobal(Graph,GC, Model, seedDegDis, SeedSize, ICProb, MC,ratio, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("DegDisGlobal",GC, seedDegDis);
	return seedDegDis;
}

vector<int>  NSDegreeDiscount::callDegDisGlobalOfLocal(const vector<int>& candidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, const double& ratio, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedDegDis;
	NSDegreeDiscount::DegDisGlobalOfLocal(candidateSet, Graph,GC, Model, seedDegDis, SeedSize, ICProb, MC, ratio, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("DegDisGL", GC, seedDegDis);
	return seedDegDis;
}

vector<int>  NSDegreeDiscount::callDegDisLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, const double& ratio) {
	vector<int> seedDegDis;
	NSDegreeDiscount::DegDisLocal(Graph,GC, Model, seedDegDis, SeedSize, ICProb, MC,ratio);
	NStools::SaveSeedSetToFileLocal("DegDisLocal", GC, seedDegDis);
	return seedDegDis;
}
