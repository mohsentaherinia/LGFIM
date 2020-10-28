#include "stdafx.h"


namespace NSICModelPF {

	template<class PGraph>
	int ICinfluenceSpreadParfull1(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull2(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull3(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull4(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull5(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull6(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull7(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull8(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull9(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull10(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull11(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull12(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull13(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull14(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull15(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull16(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull17(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull18(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull19(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull20(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull21(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull22(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull23(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParfull24(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	int ICinfluenceSpreadParallelFull(const PGraph& Graph, const vector<int> &seeds, const double &prob, const int &MC) {
		int SumAffected[24] = { 0 };
		int counter = 0;
		///////////
#pragma omp parallel sections
		{
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[0] += ICinfluenceSpreadParfull1(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[1] += ICinfluenceSpreadParfull2(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[2] += ICinfluenceSpreadParfull3(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[3] += ICinfluenceSpreadParfull4(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[4] += ICinfluenceSpreadParfull5(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[5] += ICinfluenceSpreadParfull6(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[6] += ICinfluenceSpreadParfull7(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[7] += ICinfluenceSpreadParfull8(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[8] += ICinfluenceSpreadParfull9(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[9] += ICinfluenceSpreadParfull10(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[10] += ICinfluenceSpreadParfull11(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[11] += ICinfluenceSpreadParfull12(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[12] += ICinfluenceSpreadParfull13(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[13] += ICinfluenceSpreadParfull14(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[14] += ICinfluenceSpreadParfull15(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[15] += ICinfluenceSpreadParfull16(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[16] += ICinfluenceSpreadParfull17(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[17] += ICinfluenceSpreadParfull18(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[18] += ICinfluenceSpreadParfull19(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[19] += ICinfluenceSpreadParfull20(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[20] += ICinfluenceSpreadParfull21(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[21] += ICinfluenceSpreadParfull22(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[22] += ICinfluenceSpreadParfull23(Graph, seeds, prob);
				counter++;
			}
		}
#pragma omp section
		{
			for (int i = 0; i < MC / 24; i++) {
				SumAffected[23] += ICinfluenceSpreadParfull24(Graph, seeds, prob);
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


