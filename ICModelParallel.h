#include "stdafx.h"


namespace NSICModelP {

	template<class PGraph>
	int ICinfluenceSpreadPar1(const PGraph& Graph, const vector<int> &seeds, const double &prob){
		int affected = 0;
		queue<int> q;
		TIntBoolH Completed;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Completed.AddDat(NI.GetId(), false);
		for (int i = 0; i < seeds.size(); i++)
		{
			q.push(seeds[i]);
			Completed.AddDat(seeds[i]) = true;
		}
		affected += seeds.size();
		while (!q.empty())
		{
			int curnode = q.front();
			q.pop();
			PGraph::TObj::TNodeI NI = Graph->GetNI(curnode);
			int neighbour;
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbour = NI.GetOutNId(e);

				if (Completed.GetDat(neighbour))
					continue;
				double gen = ((double)(1 + rand() % LIMIT)) / LIMIT;
				if (gen <= prob)
				{
					Completed.AddDat(neighbour, true);
					affected++;
					q.push(neighbour);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int ICinfluenceSpreadPar2(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
		int affected = 0;
		queue<int> q;
		TIntBoolH Completed;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Completed.AddDat(NI.GetId(), false);
		for (int i = 0; i < seeds.size(); i++)
		{
			q.push(seeds[i]);
			Completed.AddDat(seeds[i]) = true;
		}
		affected += seeds.size();
		while (!q.empty())
		{
			int curnode = q.front();
			q.pop();
			PGraph::TObj::TNodeI NI = Graph->GetNI(curnode);
			int neighbour;
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbour = NI.GetOutNId(e);

				if (Completed.GetDat(neighbour))
					continue;
				double gen = ((double)(1 + rand() % LIMIT)) / LIMIT;
				if (gen <= prob)
				{
					Completed.AddDat(neighbour, true);
					affected++;
					q.push(neighbour);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int ICinfluenceSpreadPar3(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
		int affected = 0;
		queue<int> q;
		TIntBoolH Completed;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Completed.AddDat(NI.GetId(), false);
		for (int i = 0; i < seeds.size(); i++)
		{
			q.push(seeds[i]);
			Completed.AddDat(seeds[i]) = true;
		}
		affected += seeds.size();
		while (!q.empty())
		{
			int curnode = q.front();
			q.pop();
			PGraph::TObj::TNodeI NI = Graph->GetNI(curnode);
			int neighbour;
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbour = NI.GetOutNId(e);

				if (Completed.GetDat(neighbour))
					continue;
				double gen = ((double)(1 + rand() % LIMIT)) / LIMIT;
				if (gen <= prob)
				{
					Completed.AddDat(neighbour, true);
					affected++;
					q.push(neighbour);
				}
			}
		}
		return affected;
	}

	template<class PGraph>
	int ICinfluenceSpreadPar4(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
		int affected = 0;
		queue<int> q;
		TIntBoolH Completed;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Completed.AddDat(NI.GetId(), false);
		for (int i = 0; i < seeds.size(); i++)
		{
			q.push(seeds[i]);
			Completed.AddDat(seeds[i]) = true;
		}
		affected += seeds.size();
		while (!q.empty())
		{
			int curnode = q.front();
			q.pop();
			PGraph::TObj::TNodeI NI = Graph->GetNI(curnode);
			int neighbour;
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbour = NI.GetOutNId(e);

				if (Completed.GetDat(neighbour))
					continue;
				double gen = ((double)(1 + rand() % LIMIT)) / LIMIT;
				if (gen <= prob)
				{
					Completed.AddDat(neighbour, true);
					affected++;
					q.push(neighbour);
				}
			}
		}
		return affected;
	}


	template<class PGraph>
	int ICinfluenceSpreadParallel4(const PGraph& Graph, const vector<int> &seeds, const double &prob, const int &MC){
		int SumAffected[4] = { 0, 0, 0, 0 };
		int counter = 0;
		///////////
#pragma omp parallel sections
		{
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[0] += ICinfluenceSpreadPar1(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[1] += ICinfluenceSpreadPar2(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[2] += ICinfluenceSpreadPar3(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 4; i++) {
				SumAffected[3] += ICinfluenceSpreadPar4(Graph, seeds, prob);
				counter++;
			}
		}
		}
		///////////
		return (SumAffected[0] + SumAffected[1] + SumAffected[2] + SumAffected[3]) / counter;
	}


}; // namespace MohsenTNT


