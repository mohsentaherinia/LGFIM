#pragma once
#include "stdafx.h"


namespace NSCOFIM {
	vector<int> callCofimGreedyGlobal(const PNGraph *SubG, const int& CmtyNum,const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callCofimCELFGlobal(const PNGraph *SubG, const int& CmtyNum, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	long CountOf_N_S(const PNGraph & Graph, const vector<int>& seed);
	long CountOf_NC_S(const PNGraph & Graph, const vector<int>& seed, const TIntIntH& NodeCmtyH, set<int>& uniqeCmty);
	inline long CountOf_NC_CurNode(const PNGraph & Graph, const int & curNode, const TIntIntH & NodeCmtyH, const set<int>& uniqeCmty);
	inline long CountOf_N_CurNode(const PNGraph & Graph, const int & vertex);


	
	template<class PGraph>
	void CofimGreedyGlobal(const TIntIntH& NodeCmtyH,const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		set<int> uniqeCmty;
		int Lambda = 3;
		long currentSeedOutNeigh, currentSeedWithcurVerOutNeigh;
		long currentSeedOutCmtyNeigh, currentSeedWithcurVerOutCmtyNeigh;
		long BestSeedOutDegree; 
		int infectedNodess=999;
		TExeTm ExeTmR;
		BestSeedOutDegree = -1;
		double ElapsedSimTimes = 0, EndTime;
		for (size_t k = 0; k < SeedSize; k++)
		{
			int HighVer, curVer;
			currentSeedOutNeigh = NSCOFIM::CountOf_N_S(Graph, seeds);
			currentSeedOutCmtyNeigh = NSCOFIM::CountOf_NC_S(Graph, seeds, NodeCmtyH, uniqeCmty);
			//////////Create A set of discovered influential Nodes in prev Step////////////
			//////////////////////////////////////////////////////////////////////////////
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				curVer = NI.GetId();
				if (NStools::findNodeInVector(seeds, curVer))//If current Vertex exist in (set of discovered influential Nodes in prev Step)=> goto next vertex
					continue;
				currentSeedWithcurVerOutNeigh = currentSeedOutNeigh + NSCOFIM::CountOf_N_CurNode(Graph, curVer);
				currentSeedWithcurVerOutCmtyNeigh = currentSeedOutCmtyNeigh + NSCOFIM::CountOf_NC_CurNode(Graph, curVer, NodeCmtyH, uniqeCmty);
				int xxx = currentSeedWithcurVerOutNeigh + Lambda * currentSeedWithcurVerOutCmtyNeigh;
				if ( xxx > BestSeedOutDegree) {
					BestSeedOutDegree = xxx;
					HighVer = curVer;
				}
			}
			seeds.push_back(HighVer);
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			if (GC._simMode == 0)
				infectedNodess = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
			cout << endl << "CofimGreedyGlobal   K=" << k + 1 << " " << infectedNodess << " is " << (double)infectedNodess / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
			NStools::SaveToFile("CofimGreedyGlobal", k + 1, infectedNodess, Graph->GetNodes(), MC, Model, exactRunTime, GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodess, Model, "COFIM", "Influence Spraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "COFIM", "Running Time (Sec)");
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		}
	}

	template<class PGraph>
	void CofimCELFGlobal(const TIntIntH& NodeCmtyH, const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;

		set<int> uniqeCmty;
		int Lambda = 5;
		long currentSeedOutNeigh, currentSeedWithcurVerOutNeigh;
		long currentSeedOutCmtyNeigh, currentSeedWithcurVerOutCmtyNeigh;
		currentSeedOutNeigh = 0;// NSCOFIM::CountOf_N_S(Graph, seeds);
		currentSeedOutCmtyNeigh = 0;// NSCOFIM::CountOf_NC_S(Graph, seeds, NodeCmtyH, uniqeCmty);
		int xxx;

		queue<TIntIntH> q;
		
		int infectedNodes=999, curVer;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			TIntIntH item1;
			curVer = NI.GetId();
		
			currentSeedWithcurVerOutNeigh = currentSeedOutNeigh + NSCOFIM::CountOf_N_CurNode(Graph, curVer);
			currentSeedWithcurVerOutCmtyNeigh = currentSeedOutCmtyNeigh + NSCOFIM::CountOf_NC_CurNode(Graph, curVer, NodeCmtyH, uniqeCmty);
			xxx = currentSeedWithcurVerOutNeigh + Lambda * currentSeedWithcurVerOutCmtyNeigh;

			item1.AddDat(curVer) = xxx;
			q.push(item1);
		}
		NScelf::QueueSortFirstK(q, q.size());
		//# Select the first node and remove from candidate list
		int k = 0;
		TIntIntH item2 = q.front();
		seeds.push_back(item2.GetKey(0));
		xxx = item2.GetDat(seeds[0]);
		EndTime = ExeTmR.GetSecs();
		double exactRunTime = EndTime - ElapsedSimTimes;
		if (GC._simMode == 0)
			infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
		q.pop();
		cout << endl << "CofimCELFGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
		NStools::SaveToFile("CofimCELFGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, exactRunTime, GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "COFIM-CELF", "Influence Spraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "COFIM-CELF", "Running Time (Sec)");
		//	# Find the next k - 1 nodes using the list - sorting procedure
		ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int Newxxx;
			currentSeedOutNeigh = NSCOFIM::CountOf_N_S(Graph, seeds);
			currentSeedOutCmtyNeigh = NSCOFIM::CountOf_NC_S(Graph, seeds, NodeCmtyH, uniqeCmty);
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = q.front();
				int current = item3.GetKey(0);
				
				q.pop();

				
				seeds.push_back(current);

				currentSeedWithcurVerOutNeigh = currentSeedOutNeigh + NSCOFIM::CountOf_N_CurNode(Graph, current);
				currentSeedWithcurVerOutCmtyNeigh = currentSeedOutCmtyNeigh + NSCOFIM::CountOf_NC_CurNode(Graph, current, NodeCmtyH, uniqeCmty);
				Newxxx = currentSeedWithcurVerOutNeigh + Lambda * currentSeedWithcurVerOutCmtyNeigh;

				seeds.pop_back();

				TIntIntH item4;
				item4.AddDat(current) = Newxxx - xxx;
				q.push(item4);
				NScelf::QueueSortOthersK(q, q.size());
				TIntIntH item5 = q.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = q.front();
			seeds.push_back(item6.GetKey(0));
			xxx = item6.GetDat(seeds[k]);
			//seeds[k] = item6.GetKey(0);
			EndTime = ExeTmR.GetSecs();
		    exactRunTime = EndTime - ElapsedSimTimes;
			if (GC._simMode == 0)
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
			q.pop();
			cout << endl << "CofimCELFGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
			NStools::SaveToFile("CofimCELFGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, exactRunTime, GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "COFIM-CELF", "Influence Spraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "COFIM-CELF", "Running Time (Sec)");
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		}//end for
	
	}

	



}