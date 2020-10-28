#pragma once
#include "stdafx.h"
#define IcProbCGA .5
#define MC_CGA 6
#define theta 0.1


namespace NSCGA {
	void integrate(const int & cm, const int & cl, TCnComV & CmtyV);
	bool IsexistVertexInL(const int & u, const set<int> L[], const int & cm);
	void DifferenceVector(vector<int>& S1, vector<int>& S2, vector<int>& S3);
	int getCmtyLAbelOfvertex(const TCnComV & CmtyV, const int & vertex);

	template<class PGraph>
	int isInfluenceICinfluenceSpread(const PGraph& Graph, const vector<int> &seeds, const double &prob, const int & u) {
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
					if (neighbour == u)
						return 1;
				}
			}
		}
		return 0;
	}

	template<class PGraph>
	bool isInfluence(const int&v, const int&u, const PGraph& Graph) {
		int SumAffected = 0;
		vector<int> seeds;
		seeds.push_back(v);
		for (size_t i = 0; i < MC_CGA; i++)
			SumAffected += isInfluenceICinfluenceSpread(Graph, seeds, IcProbCGA, u);
		if (SumAffected >(MC_CGA / 2))
			return true;
		return false;
	}

	template<class PGraph>
	double callCGA(const PGraph& Graph) {
		TIntIntVH H;
		TIntIntH labelV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int v = NI.GetId();
			labelV.AddDat(v, v);
			TIntV NodeVector;
			H.AddDat(v) = NodeVector;
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int u = NI.GetOutNId(e);
				if (isInfluence(v, u, Graph))
				{
					NodeVector = H.GetDat(v);
					NodeVector.Add(u);
					H.AddDat(v) = NodeVector;
				}
			}
		}
		////////////////////////////////////////////////////////////
		int counter = 0;
		while (counter++ < ConvergeRep) {
			TIntIntH templabelV;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				int CurNode = NI.GetId();
				vector <int> NLV;
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					int neighbourID = NI.GetOutNId(e);
					if (H.GetDat(CurNode).SearchBin(neighbourID) != -1)
						NLV.push_back(labelV.GetDat(neighbourID));
				}
				int newlabel = NSlpa::findMaximumLabelorg(NLV);
				templabelV.AddDat(CurNode) = newlabel;
			}
			NSlpa::updateLabelTable(labelV, templabelV);
		}
		//NStools::PrintHaShTable(labelV, labelV.Len());
		TCnComV CmtyV;
		NSlpa::ExtractCommunity2(Graph, labelV, CmtyV);
		NSlpa::printCommunityinfo(CmtyV);
		combinationStep(Graph, CmtyV, H);
		NSlpa::printCommunityinfo(CmtyV);
		double Q = 0;
		Q = NSlpa::GetModularity2(Graph, CmtyV);
		return Q;
	}

	template<class PGraph>
	void combinationStep(const PGraph& Graph, TCnComV& CmtyV, const TIntIntVH& H) {
		//for calculate coEntropy function
		vector<int> graphVertexSet;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			graphVertexSet.push_back(NI.GetId());
		}
		/////////
		const int cmtyCount = CmtyV.Len();
		set<int>*L = new set<int>[cmtyCount]; //List of adj vertex of other community neigbour
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int vi = NI.GetId();
			int cmtyvi = getCmtyLAbelOfvertex(CmtyV, vi);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int uj = NI.GetOutNId(e);
				int cmtyuj = getCmtyLAbelOfvertex(CmtyV, uj);
				if (cmtyvi == -1 || cmtyuj == -1)
					break;
				if (cmtyvi != cmtyuj) //v and j  bieng to different community
					if (H.GetDat(vi).SearchBin(uj) != -1)
						L[cmtyvi].insert(uj);
			}
		}
		//////////////
		//////////////
		//////////////
		bool IsComb = true;
		while (IsComb) {
			for (int cm = 0; cm < CmtyV.Len(); cm++) {
				set<int> s;
				for (set<int>::iterator u = L[cm].begin(); u != L[cm].end(); ++u) {
					int cmtylabel = getCmtyLAbelOfvertex(CmtyV, *u);
					s.insert(cmtylabel);
				}
				IsComb = false;
				for (set<int>::iterator cl = s.begin(); cl != s.end(); ++cl) {
					double coEntropyLM = coEntropyLMF(Graph, CmtyV, L, cm, *cl, graphVertexSet);
					if (coEntropyLM > theta) {
						if (cm == *cl)
							continue;
						integrate(cm, *cl, CmtyV);
						NSlpa::printCommunityinfo(CmtyV);
						IsComb = true;
					}
				}
			}
		}
		//	NStools::PrintVectorOfVector(L, cmtyCount);
		delete[]L;
	}

	template<class PGraph>
	double coEntropyLMF(const PGraph& Graph, const TCnComV& CmtyV, const set<int>L[], const int& cm, const int& cl, vector<int>& graphVertexSet) {
		double maxEntropy = 0.0;
		vector<int> cmtyCMvertexSet;
		for (size_t ccc = 0; ccc < CmtyV[cm].Len(); ccc++)
			cmtyCMvertexSet.push_back(CmtyV[cm][ccc]);
		for (size_t v = 0; v < CmtyV[cm].Len(); v++) {
			vector<int> seedv;
			seedv.push_back(CmtyV[cm][v]);
			double Rmv = CustomICinfluenceSpread(Graph, seedv, cmtyCMvertexSet);
			Rmv /= (double)cmtyCMvertexSet.size();

			for (size_t u = 0; u < CmtyV[cl].Len(); u++) {
				if (IsexistVertexInL(CmtyV[cl][u], L, cm)) {
					vector<int> seedu, OutCMVertexSet;
					seedu.push_back(CmtyV[cl][u]);
					//calc vertexset Outside cm
					DifferenceVector(graphVertexSet, cmtyCMvertexSet, OutCMVertexSet);
					double Rmu = CustomICinfluenceSpread(Graph, seedu, OutCMVertexSet);
					Rmu /= (double)OutCMVertexSet.size();
					double div = Rmu / Rmv;
					if (div > maxEntropy)
						maxEntropy = div;
				}//end of if
			}//End Of Inner for
		}//End Of Outer for
		return maxEntropy;
	}

	template<class PGraph>
	int CustomICinfluenceSpreadin(const PGraph& Graph, const vector<int> &seeds, const double &prob) {
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
	double CustomICinfluenceSpread(const PGraph& Graph, const vector<int> &seeds, vector<int> &vertexSet) {
		TIntV NIdV;
		for (vector<int>::iterator it = vertexSet.begin(); it != vertexSet.end(); ++it)
			NIdV.Add(*it);
		if (NIdV.SearchBin(seeds[0]) == -1)
			NIdV.Add(seeds[0]);
		PNGraph SubG = TSnap::GetSubGraph(Graph, NIdV);
		
		double SumAffected = 0;
		for (size_t i = 0; i < MC_CGA; i++)
			SumAffected += CustomICinfluenceSpreadin(SubG, seeds, IcProbCGA);
		return SumAffected / (double)MC_CGA;
	}


}