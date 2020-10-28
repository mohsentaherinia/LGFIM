#pragma once
#include "stdafx.h"


namespace NSgreedy {
	vector<int> callGreedyGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callGreedyGlobalOfLocal(const vector<int>& candidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callGreedyLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC);

	
	template<class PGraph>
	void GreedyGlobal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		int BestSpread, infectedNodes;
		TExeTm ExeTm;
		for (size_t k = 0; k < SeedSize; k++)
		{
			int HighVer, curVer;
			BestSpread = -1, infectedNodes = 0;
			//////////////////////////////////////////////////////////////////////////////
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				curVer = NI.GetId();
				if (NStools::findNodeInVector(seeds, curVer))//If current Vertex exist in (set of discovered influential Nodes in prev Step)=> goto next vertex
					continue;
				seeds.push_back(curVer);
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				if (infectedNodes > BestSpread) {
					BestSpread = infectedNodes;
					HighVer = curVer;
				}
			}
			seeds.push_back(HighVer);
			cout << endl << "GreedyGlobal   K=" << k + 1 << " " << BestSpread << " is " << (double)BestSpread / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm.GetSecs()<<" Sec";
			NStools::SaveToFile("GreedyGlobal", k + 1, BestSpread, Graph->GetNodes(), MC, Model, ExeTm.GetSecs(), GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, BestSpread, Model, "GreedyGlobal", "InfluenceSpraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm.GetSecs(), Model, "GreedyGlobal", "RunningTime");
		}
	}

	template<class PGraph>
	void GreedyGlobalOfLocal(const vector<int> &candidateSet, const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTm;
		int BestSpread, infectedNodes;
		for (size_t k = 0; k < SeedSize; k++)
		{
			int HighVer, curVer;
			BestSpread = -1, infectedNodes = 0;
			//////////////////////////////////////////////////////////////////////////////
			//for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			//curVer = NI.GetId();
			for (size_t i = 0; i < candidateSet.size(); i++) {
				curVer = candidateSet[i];
				if (NStools::findNodeInVector(seeds, curVer))//If current Vertex exist in (set of discovered influential Nodes in prev Step)=> goto next vertex
					continue;
				seeds.push_back(curVer);
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				if (infectedNodes > BestSpread) {
					BestSpread = infectedNodes;
					HighVer = curVer;
				}
			}
			seeds.push_back(HighVer);
			cout << endl << "GreedyGL   K=" << k + 1 << " " << BestSpread << " is " << (double)BestSpread / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm.GetSecs()<<" Sec";
			NStools::SaveToFile("GreedyGL", k + 1, BestSpread, Graph->GetNodes(), MC, Model, ExeTm.GetSecs(), GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, BestSpread, Model, "GreedyGL", "InfluenceSpraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm.GetSecs(), Model, "GreedyGL", "RunningTime");
		}
	}

	template<class PGraph>
	void GreedyLocal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC) {
		TExeTm ExeTm; 
		int BestSpread, infectedNodes;
		for (size_t k = 0; k < SeedSize; k++)
		{
			int HighVer, curVer;
			BestSpread = -1, infectedNodes = 0;
			//////////////////////////////////////////////////////////////////////////////
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				curVer = NI.GetId();
				if (NStools::findNodeInVector(seeds, curVer))//If current Vertex exist in (set of discovered influential Nodes in prev Step)=> goto next vertex
					continue;
				seeds.push_back(curVer);
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				if (infectedNodes > BestSpread) {
					BestSpread = infectedNodes;
					HighVer = curVer;
				}
			}
			seeds.push_back(HighVer);
			NStools::SaveToFile("GreedyLocal", k + 1, BestSpread, Graph->GetNodes(), MC, Model, ExeTm.GetSecs(), GC);
		}
	}
	
}