#include "stdafx.h"


namespace NSLTModelP {

	template<class PGraph>
	int LTinfluenceSpreadPar1(const PGraph& Graph, const vector<int> &seeds) {
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
				ActNeighbourCntV.AddDat(curNode) = ActNeibCnt + 1;
				ActNeibCnt++;

				PGraph::TObj::TNodeI NI2 = Graph->GetNI(curNode);
				sumThrsActNeibs = (float)ActNeibCnt / NI2.GetInDeg();

				float tettaValue = tetaV.GetDat(curNode);

				if (sumThrsActNeibs >= tettaValue)
				{
					ActivedV.AddDat(curNode, true);
					affected++;
					q.push(curNode);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int LTinfluenceSpreadPar2(const PGraph& Graph, const vector<int> &seeds) {
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
				ActNeighbourCntV.AddDat(curNode, ActNeibCnt + 1);
				ActNeibCnt++;

				PGraph::TObj::TNodeI NI2 = Graph->GetNI(curNode);
				sumThrsActNeibs = (float)ActNeibCnt / NI2.GetInDeg();

				float tettaValue = tetaV.GetDat(curNode);

				if (sumThrsActNeibs >= tettaValue)
				{
					ActivedV.AddDat(curNode, true);
					affected++;
					q.push(curNode);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int LTinfluenceSpreadPar3(const PGraph& Graph, const vector<int> &seeds) {
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
				ActNeighbourCntV.AddDat(curNode, ActNeibCnt + 1);
				ActNeibCnt++;

				PGraph::TObj::TNodeI NI2 = Graph->GetNI(curNode);
				sumThrsActNeibs = (float)ActNeibCnt / NI2.GetInDeg();

				float tettaValue = tetaV.GetDat(curNode);

				if (sumThrsActNeibs >= tettaValue)
				{
					ActivedV.AddDat(curNode, true);
					affected++;
					q.push(curNode);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int LTinfluenceSpreadPar4(const PGraph& Graph, const vector<int> &seeds) {
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
				ActNeighbourCntV.AddDat(curNode, ActNeibCnt + 1);
				ActNeibCnt++;

				PGraph::TObj::TNodeI NI2 = Graph->GetNI(curNode);
				sumThrsActNeibs = (float)ActNeibCnt / NI2.GetInDeg();

				float tettaValue = tetaV.GetDat(curNode);

				if (sumThrsActNeibs >= tettaValue)
				{
					ActivedV.AddDat(curNode, true);
					affected++;
					q.push(curNode);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int LTinfluenceSpreadParallel4(const PGraph& Graph, const vector<int> &seeds, const int &MC) {
		int SumAffected[4] = { 0, 0, 0, 0 };
		int counter = 0;
		///////////
#pragma omp parallel sections
		{
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[0] += LTinfluenceSpreadPar1(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[1] += LTinfluenceSpreadPar2(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[2] += LTinfluenceSpreadPar3(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[3] += LTinfluenceSpreadPar4(Graph, seeds);
				counter++;
			}
		}


		}
		///////////
		return (SumAffected[0] + SumAffected[1] + SumAffected[2] + SumAffected[3]) / counter;
	}


}; // namespace MohsenTNT

