#include "stdafx.h"

namespace NSLTModel {
	

	template<class PGraph>
	int LTinfluenceSpread1(const PGraph& Graph, const vector<int> &seeds) {
		int affected = 0;
		double sumThrsActNeibs;
		queue<int> q;
		TIntBoolH ActivedV;
		TIntFltH tetaV;
		TIntIntH  ActNeighbourCntV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			ActivedV.AddDat(NI.GetId(), false);
			tetaV.AddDat(NI.GetId(), ((double)(1 + rand() % LIMIT)) / LIMIT);
			ActNeighbourCntV.AddDat(NI.GetId(), 0);
		}
		for (int i = 0; i < seeds.size(); i++)
		{
			q.push(seeds[i]);
			ActivedV.AddDat(seeds[i]) = true;
		}
		affected += seeds.size();
		////////////////////////////////////////////////////////////////////////////////////////
		while (!q.empty()) {
			int activatedNode = q.front();
			q.pop();
			PGraph::TObj::TNodeI NI = Graph->GetNI(activatedNode);
			//////////////////baraye tamame hamsayehaye node fa'al ke az safe fa'alha dar amade
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int curNode = NI.GetOutNId(e); 
				if (ActivedV.GetDat(curNode))
					continue;
				
				int ActNeibCnt = ActNeighbourCntV.GetDat(curNode);
				ActNeighbourCntV.AddDat(curNode)= ActNeibCnt + 1;
				ActNeibCnt++;


				PGraph::TObj::TNodeI NI2 = Graph->GetNI(curNode);
				sumThrsActNeibs = (float)ActNeibCnt / NI2.GetInDeg();
				
				float tettaValue = tetaV.GetDat(curNode);
				
				if (sumThrsActNeibs >= tettaValue)
				{
					ActivedV.AddDat(curNode) = true;
					affected++;
					q.push(curNode);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int LTinfluenceSpread(const PGraph& Graph, const vector<int> &seeds, const int &MC) {
		int SumAffected = 0, counter = 0;;
		for (size_t i = 0; i < MC; i++) {
			SumAffected += LTinfluenceSpread1(Graph, seeds);
			counter++;
		}
		return SumAffected / counter;
	}

}; // namespace MohsenTNT


