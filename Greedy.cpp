#include "stdafx.h"

vector<int>  NSgreedy::callGreedyGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedGreedy;
	NSgreedy::GreedyGlobal(Graph, GC, Model ,seedGreedy, SeedSize, ICProb, MC,curveInfoISV,curveInfoRTV,indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("GreedyGlobal", GC, seedGreedy);
	return seedGreedy;
}

vector<int>  NSgreedy::callGreedyGlobalOfLocal(const vector<int>& candidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedGreedy;
	NSgreedy::GreedyGlobalOfLocal(candidateSet,Graph, GC, Model, seedGreedy, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("GreedyGL", GC, seedGreedy);
	return seedGreedy;
}

vector<int>  NSgreedy::callGreedyLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC) {
	vector<int> seedGreedy;
	NSgreedy::GreedyLocal(Graph, GC, Model, seedGreedy, SeedSize, ICProb, MC);
	NStools::SaveSeedSetToFileLocal("GreedyLocal", GC, seedGreedy);
	return seedGreedy;
}


