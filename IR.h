#pragma once
#include "stdafx.h"
#define NUM_LOOP_IR  20



namespace NSIR {
	vector<int> callIRGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callIRLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha);
	vector<int> callIRGlobalofLocal(const vector<int>& CandidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	
	template<class PGraph>
	void SortedIR(const PGraph& Graph, TIntFltH& IRRankH, const double& alpha) {
		NSIR::GetIR1(Graph, IRRankH, alpha);
		//IRRankH.SortByDat(false); //Decressing data value  //Asc
								  //MohsenTNT::PrintHaShTable(IRRankH, Graph->GetNodes());
	}

	template<class PGraph>
	void GetIR1(const PGraph& Graph, TIntFltH& IRRankH, const double& alpha)
	{
		int NNodes = Graph->GetNodes();
		IRRankH.Gen(NNodes);
		bool changed = true;
		bool saturated = false;
		int i = 0;
		double theta = 0.0001;
		double edgeProb;
		long count = 0;

		TIntFltH dpV, newdpV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			dpV.AddDat(NI.GetId(), 1.0);
			newdpV.AddDat(NI.GetId(), 0.0);
			IRRankH.AddDat(NI.GetId(), 0.0);
		}


		while (count++ < NUM_LOOP_IR && changed && !saturated)
		{
			changed = false;

			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				newdpV.AddDat(NI.GetId()) = 0.0;
			}


			float value;
			int curNodeiuID, neighbourjvID;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				curNodeiuID = NI.GetId();
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					neighbourjvID = NI.GetOutNId(e);
					PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
					edgeProb = (float)1.0 / NI2.GetInDeg();
					if (edgeProb)
					{
						//newdp[i] += alpha * edgeProb * dp[e.v];
						value = newdpV.GetDat(curNodeiuID);
						value += alpha * edgeProb * dpV.GetDat(neighbourjvID);
						newdpV.AddDat(curNodeiuID) = value;
					}

				}
				//	newdp[i] += 1;
				newdpV.AddDat(curNodeiuID) = newdpV.GetDat(curNodeiuID) + 1;
			}
			/////////////////////////////////////////////////////////////////////////////////////
			int curNodei;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				curNodei = NI.GetId();
				//if((newdp[i] < dp[i] - theta) || (newdp[i] > dp[i] + theta))
				if ((newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) < -theta) ||
					(newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) > theta))
					changed = true;
				//if(newdp[i] > NNodes)
				if (newdpV.GetDat(curNodei) > NNodes)
					saturated = true;
				//dp[i] = newdp[i];
				dpV.AddDat(curNodei) = newdpV.GetDat(curNodei);
			}
		}
		IRRankH = dpV;
	}
	
}
