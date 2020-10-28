#include "stdafx.h"

namespace NSLTModelPF {

	template<class PGraph>
	int LTinfluenceSpreadParfull1(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull2(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull3(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull4(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull5(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull6(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull7(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull8(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull9(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull10(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull11(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull12(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull13(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull14(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull15(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull16(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull17(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull18(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull19(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull20(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull21(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull22(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull23(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParfull24(const PGraph& Graph, const vector<int> &seeds) {
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
	int LTinfluenceSpreadParallelFull(const PGraph& Graph, const vector<int> &seeds, const int &MC) {
		int SumAffected[24] = { 0 };
		int counter = 0;
		///////////
#pragma omp parallel sections
		{
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[0] += LTinfluenceSpreadParfull1(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[1] += LTinfluenceSpreadParfull2(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[2] += LTinfluenceSpreadParfull3(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[3] += LTinfluenceSpreadParfull4(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[4] += LTinfluenceSpreadParfull5(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[5] += LTinfluenceSpreadParfull6(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[6] += LTinfluenceSpreadParfull7(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[7] += LTinfluenceSpreadParfull8(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[8] += LTinfluenceSpreadParfull9(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[9] += LTinfluenceSpreadParfull10(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[10] += LTinfluenceSpreadParfull11(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[11] += LTinfluenceSpreadParfull12(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[12] += LTinfluenceSpreadParfull13(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[13] += LTinfluenceSpreadParfull14(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[14] += LTinfluenceSpreadParfull15(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[15] += LTinfluenceSpreadParfull16(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[16] += LTinfluenceSpreadParfull17(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[17] += LTinfluenceSpreadParfull18(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[18] += LTinfluenceSpreadParfull19(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[19] += LTinfluenceSpreadParfull20(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[20] += LTinfluenceSpreadParfull21(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[21] += LTinfluenceSpreadParfull22(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[22] += LTinfluenceSpreadParfull23(Graph, seeds);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[23] += LTinfluenceSpreadParfull24(Graph, seeds);
				counter++;
			}
		}


		}
		///////////
		int sum = 0;
		for (size_t k = 0; k < 24; k++)
		{
			sum += SumAffected[k];
		}
		return sum / counter;
	}


}; // namespace MohsenTNT

