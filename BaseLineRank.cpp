#include "stdafx.h"
vector<int> NSbaseLineRank::callPageRank1Global(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntFltH PRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MyPageRank1(Graph, PRankH, 0.85, 0.0001, 100);
	PRankH.SortByDat(false);
	vector<int> seedPR;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedPR.push_back(PRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedPR, Model, ICProbb, MCS);
		cout << endl << "PageRank1Global   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent"<< exactRunTime;
		NStools::SaveToFile("PageRank1Global", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "PageRank1Global", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "PageRank1Global", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("PageRank1Global", GC, seedPR);
	return seedPR;
}

vector<int> NSbaseLineRank::callPageRank1Local(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS) {
	TIntFltH PRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MyPageRank1(Graph, PRankH, 0.85, 0.0001, 100);
	PRankH.SortByDat(false);
	vector<int> seedPR;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedPR.push_back(PRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seedPR, Model, ICProbb, MCS);
		NStools::SaveToFile("callPageRank1Local", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NStools::SaveSeedSetToFileLocal("PageRank1Local", GC, seedPR);
	return seedPR;
}

vector<int> NSbaseLineRank::callPageRank1GlobalofLocal(const vector<int>& CandidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntFltH PRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MyPageRank1(Graph, PRankH, 0.85, 0.0001, 100);
	TIntFltH NewRankH;
	for (int i = 0; i < CandidateSet.size(); i++)
	{
		NewRankH.AddDat(CandidateSet[i], PRankH.GetDat((CandidateSet[i])));
	}
	NewRankH.SortByDat(false);
	vector<int> seedPR;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedPR.push_back(NewRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedPR, Model, ICProbb, MCS);
		cout << endl << "PageRank1GL   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("PageRank1GL", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "PageRank1GL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "PageRank1GL", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("PageRank1GL", GC, seedPR);
	return seedPR;
}

vector<int> NSbaseLineRank::callPageRank2Global(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntFltH PRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MyPageRank2(Graph, PRankH, 0.85, 0.0001, 100);
	PRankH.SortByDat(false);
	vector<int> seedPR;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedPR.push_back(PRankH.GetKey(k-1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime =  EndTime - ElapsedSimTimes;
		int infectedNodes = 999;
		if (GC._simMode == 0)
			infectedNodes =  NSIS::callInfluenceSpreadModel(Graph, seedPR, Model, ICProbb, MCS);
		cout << endl << "PageRank2Global   K=" << k  << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent " << exactRunTime;
		NStools::SaveToFile("PageRank2Global", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "PageRank", "Influence Spraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "PageRank", "Running Time (Sec)");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("PageRank2Global", GC, seedPR);
	return seedPR;
}

vector<int> NSbaseLineRank::callPageRank2Local(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS) {
	TIntFltH PRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MyPageRank2(Graph, PRankH, 0.85, 0.0001, 100);
	PRankH.SortByDat(false);
	vector<int> seedPR;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedPR.push_back(PRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seedPR, Model, ICProbb, MCS);
		NStools::SaveToFile("callPageRank2Local", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NStools::SaveSeedSetToFileLocal("PageRank2Local", GC, seedPR);
	return seedPR;
}

vector<int> NSbaseLineRank::callPageRank2GlobalofLocal(const vector<int>& CandidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntFltH PRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MyPageRank2(Graph, PRankH, 0.85, 0.0001, 100);
	TIntFltH NewRankH;
	for (int i = 0; i < CandidateSet.size(); i++)
	{
		NewRankH.AddDat(CandidateSet[i], PRankH.GetDat((CandidateSet[i])));
	}
	NewRankH.SortByDat(false);

	vector<int> seedPR;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedPR.push_back(NewRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedPR, Model, ICProbb, MCS);
		cout << endl << "PageRank2GL   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent"<< exactRunTime;
		NStools::SaveToFile("PageRank2GL", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "PageRank2GL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "PageRank2GL", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("PageRank2GL", GC, seedPR);
	return seedPR;
}

vector<int> NSbaseLineRank::callMaxDegreeRankGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MaximimDegree(Graph, MRankH);
	MRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(MRankH.GetKey(k-1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;
		if (GC._simMode == 0)
			 infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "MaxDegreeGlobal   K=" << k  << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("MaxDegreeGlobal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "Max-Degree", "Influence Spraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "Max-Degree", "Running Time (Sec)");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("MaxDegreeGlobal", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callMaxDegreeRankLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MaximimDegree(Graph, MRankH);
	MRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(MRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		NStools::SaveToFile("MaxDegreeLocal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NStools::SaveSeedSetToFileLocal("MaxDegreeLocal", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callMaxDegreeRankGlobalofLocal(const vector<int>& CandidateSet,const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime; 
	NSbaseLineRank::MaximimDegree(Graph, MRankH);
	TIntIntH NewRankH;
	for (int i = 0; i < CandidateSet.size(); i++)
	{
		NewRankH.AddDat(CandidateSet[i], MRankH.GetDat((CandidateSet[i])));
	}
	NewRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(NewRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "MaxDegreeGL   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("MaxDegreeGL", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "MaxDegreeGL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "MaxDegreeGL", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("MaxDegreeGL", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callRandomCentrality(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		while (true) {
			int RandomNode = Graph->GetRndNId();
			if (!NStools::findNodeInVector(seedSet, RandomNode)) {// ignarance Duplicated Node
				seedSet.push_back(RandomNode);
				break;
			}
		}
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "Random   K=" << k  << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("Random", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "Random", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "Random", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("Random", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callRandomCentrality_v2(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS) {
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		while (true) {
			int RandomNode = Graph->GetRndNId();
			if (!NStools::findNodeInVector(seedSet, RandomNode)) {// ignarance Duplicated Node
				seedSet.push_back(RandomNode);
				break;
			}
		}
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		NStools::SaveToFile("Random_v2", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	return seedSet;
}

vector<int> NSbaseLineRank::callMinDegreeRankGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MinimumDegree(Graph, MRankH);
	MRankH.SortByDat(true);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(MRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "MinDegreeGlobal   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("MinDegreeGlobal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "MinDegreeGlobal", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "MinDegreeGlobal", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("MinDegreeGlobal", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callMinDegreeRankLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MinimumDegree(Graph, MRankH);
	MRankH.SortByDat(true);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(MRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		NStools::SaveToFile("MinDegreeLocal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NStools::SaveSeedSetToFileLocal("MinDegreeLocal", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callMinDegreeRankGlobalofLocal(const vector<int>& CandidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::MinimumDegree(Graph, MRankH);
	TIntIntH NewRankH;
	for (int i = 0; i < CandidateSet.size(); i++)
	{
		NewRankH.AddDat(CandidateSet[i], MRankH.GetDat((CandidateSet[i])));
	}
	NewRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(NewRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "MinDegreeGL   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("MinDegreeGL", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "MinDegreeGL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "MinDegreeGL", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("MinDegreeGL", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callKunduRankGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::KunduRank(Graph, MRankH);
	MRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(MRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "KunduGlobal   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("KunduGlobal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "KunduGlobal", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "KunduGlobal", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("KunduGlobal", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callKunduRankLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::KunduRank(Graph, MRankH);
	MRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(MRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		NStools::SaveToFile("KunduLocal", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NStools::SaveSeedSetToFileLocal("KunduLocal", GC, seedSet);
	return seedSet;
}

vector<int> NSbaseLineRank::callKunduRankGlobalofLocal(const vector<int>& CandidateSet, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int &indexPlot) {
	TIntIntH MRankH;
	TExeTm ExeTmR;
	double ElapsedSimTimes = 0, EndTime;
	NSbaseLineRank::KunduRank(Graph, MRankH);
	TIntIntH NewRankH;
	for (int i = 0; i < CandidateSet.size(); i++)
	{
		NewRankH.AddDat(CandidateSet[i], MRankH.GetDat((CandidateSet[i])));
	}
	NewRankH.SortByDat(false);
	vector<int> seedSet;
	for (size_t k = 1; k <= SeedSize; k++)
	{
		seedSet.push_back(NewRankH.GetKey(k - 1));
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedSet, Model, ICProbb, MCS);
		cout << endl << "KunduGL   K=" << k << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
		NStools::SaveToFile("KunduGL", k, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k, infectedNodes, Model, "KunduGL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k, exactRunTime, Model, "KunduGL", "RunningTime");
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	}
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("KunduGL", GC, seedSet);
	return seedSet;
}
