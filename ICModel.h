#include "stdafx.h"
#define LIMIT 10000

namespace NSICModel {

	template<class PGraph>
	int ICinfluenceSpread(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
		int affected = 0;
		queue<int> q;
		TIntBoolH Completed;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Completed.AddDat(NI.GetId(), false);
		for (int i = 0; i < seeds.size(); i++)
		{
			q.push(seeds[i]);
			Completed.AddDat(seeds[i])=true;
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
	int ICinfluenceSpread(const PGraph& Graph, const vector<int> &seeds, const double &prob, const int &MC) {
		int SumAffected = 0;
		int counter = 0;
		for (size_t i = 0; i < MC; i++) {
			SumAffected += ICinfluenceSpread(Graph, seeds, prob);
			counter++;
		}
		return SumAffected / counter;
	}
};

