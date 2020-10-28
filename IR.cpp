#include "stdafx.h"


vector<int> NSIR::callIRGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntFltH IRRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSIR::SortedIR(Graph, IRRankH, alpha);
	IRRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(IRRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "IRGlobal   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("IRGlobal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "IRGlobal", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "IRGlobal", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("IRGlobal", GC, seedSet);
	return seedSet;
}

vector<int> NSIR::callIRLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha) {
	TIntFltH IRRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSIR::SortedIR(Graph, IRRankH, alpha);
	IRRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(IRRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		NStools::SaveToFile("IRLocal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		
	}
	
	NStools::SaveSeedSetToFileLocal("IRLocal", GC, seedSet);
	return seedSet;
}

vector<int> NSIR::callIRGlobalofLocal(const vector<int>& CandidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntFltH IRRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSIR::SortedIR(Graph, IRRankH, alpha);
	TIntFltH NewRankH;
	for (int i = 0; i < CandidateSet.size(); i++)
	{
		NewRankH.AddDat(CandidateSet[i], IRRankH.GetDat((CandidateSet[i])));
	}
	NewRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(NewRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "IRGL   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("IRGL", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "IRGL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "IRGL", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("IRGL", GC, seedSet);
	return seedSet;
}
