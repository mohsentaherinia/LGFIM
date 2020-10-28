#pragma once
#include "stdafx.h"
#define ConvergeRep 10
namespace NSlpa {
	int findMaximumLabelorg(vector<int> arr);
	int findMaximumLabel(vector<int> arr, vector<int> vertexCoverSet);
	void updateLabelTable(TIntIntH & org, const TIntIntH & temp);
	void printCommunityinfo(const TCnComV & CmtyV);
	
	
	template<class PGraph>
	double callIMLPAdirected(const PGraph& Graph) {
		vector<int> vertexSet, coreSet;
		TIntBoolH visited;
		TIntIntH MaximumDegree;
		TIntIntH labelV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			vertexSet.push_back(NI.GetId());
			MaximumDegree.AddDat(NI.GetId(), NI.GetOutDeg());
			visited.AddDat(NI.GetId(), false);
			labelV.AddDat(NI.GetId(), -1);
		}
		MaximumDegree.SortByDat(false);  //true == increasing  data value  Afzayeshii
		//NStools::PrintHaShTable(MaximumDegree, MaximumDegree.Len());
		////////////////////////////////////////////////////
		int i = 0;
		while (i < Graph->GetNodes()) {
			int currentCoreNode = MaximumDegree.GetKey(i++);
			if (!visited.GetDat(currentCoreNode)) {
				coreSet.push_back(currentCoreNode);
				visited.AddDat(currentCoreNode, true);
				PGraph::TObj::TNodeI NI = Graph->GetNI(currentCoreNode);
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					int neighbourID = NI.GetOutNId(e);
					visited.AddDat(neighbourID, true);
				}
			}
		}
		//NStools::PrintVector(coreSet);
		//////////////////////////////////////////////////////////
		for (vector<int>::iterator it = coreSet.begin(); it != coreSet.end(); ++it)
			labelV.AddDat(*it) = *it;
		//NStools::PrintHaShTable(labelV, labelV.Len());
		////////////////////////////////////////////////////////////
		TIntIntH templabelV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int CurNode = NI.GetId();
			if (labelV.GetDat(CurNode) != -1)
				continue;
			vector <int> NLV;
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbourID = NI.GetOutNId(e);
				NLV.push_back(labelV.GetDat(neighbourID));
			}
			int newlabel = findMaximumLabel(NLV, coreSet);
			if (newlabel != -1)
				templabelV.AddDat(CurNode) = newlabel;
		}
		updateLabelTable(labelV, templabelV);
		//NStools::PrintHaShTable(labelV, labelV.Len());
		////////////////////////////////////////////////////////////
		int counter = 0;
		while (counter++ < ConvergeRep) {
			TIntIntH templabelV;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				int CurNode = NI.GetId();
				vector <int> NLV;
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					int neighbourID = NI.GetOutNId(e);
					NLV.push_back(labelV.GetDat(neighbourID));
				}
				int newlabel = findMaximumLabel(NLV, coreSet);
				if (newlabel != -1)
					templabelV.AddDat(CurNode) = newlabel;
			}
			updateLabelTable(labelV, templabelV);
		}
		TCnComV CmtyV;
		ExtractCommunity2(Graph, labelV, CmtyV);
		printCommunityinfo(CmtyV);
		double Q = 0;
		Q = GetModularity2(Graph, CmtyV);
		NScommunity::SaveCommunityInfoTofileDetails(Graph, CmtyV, "dataset0", "LPADirected", Q);
		NScommunity::SaveCommunityInfoTofile(CmtyV, "dataset0", "LPADirected");
		return Q;

	}
			
	template<class PGraph>
	double callLPA2(const PGraph& Graph) {
		TIntIntH labelV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			labelV.AddDat(NI.GetId(), NI.GetId());
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

					NLV.push_back(labelV.GetDat(neighbourID));
				}
				int newlabel = findMaximumLabelorg(NLV);
				if (newlabel != -1)
					templabelV.AddDat(CurNode) = newlabel;
			}
			updateLabelTable(labelV, templabelV);
		}
		//NStools::PrintHaShTable(labelV, labelV.Len());
		TCnComV CmtyV;
		ExtractCommunity2(Graph, labelV, CmtyV);
		printCommunityinfo(CmtyV);
		double Q = 0;
		Q = GetModularity2(Graph, CmtyV);
		return Q;
	}

	template<class PGraph>
	void ExtractCommunity2(const PGraph& Graph,const TIntIntH &labelV, TCnComV& CmtyV) {
		//labelV.SortByDat();
		TIntIntVH s;
		for (int i = 0; i < labelV.Len(); i++) {
			TIntV t;
			s.AddDat(labelV.GetDat(labelV.GetKey(i))) = t;
		}
		
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int l = labelV.GetDat(NI.GetId());
			TIntV NodeVector;
			NodeVector = s.GetDat(l);
			NodeVector.Add(NI.GetId());
			s.AddDat(l) = NodeVector;
		}

		for (size_t i = 0; i < s.Len(); i++){
			int _key = s.GetKey(i);
			TCnCom t;
			t = s.GetDat(_key);
			CmtyV.Add(t);
		}
	}
	
	template<class PGraph>
	double GetModularityLescovich(const PGraph& Graph,TCnComV& CnComV) {
		TIntH OutDegH;
		const int OrigEdges = Graph->GetEdges();
		for (PGraph::TObj::TNodeI  NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			OutDegH.AddDat(NI.GetId(), NI.GetOutDeg());
		}
		///
		double Mod = 0;
		for (int c = 0; c < CnComV.Len(); c++) {
			const TIntV& NIdV = CnComV[c]();
			double EIn = 0, EEIn = 0;
			for (int i = 0; i < NIdV.Len(); i++) {
				PGraph::TObj::TNodeI  NI = Graph->GetNI(NIdV[i]);
				EIn += NI.GetOutDeg();
				EEIn += OutDegH.GetDat(NIdV[i]);
			}
			Mod += (EIn - EEIn*EEIn / (2.0*OrigEdges));
		}
		if (Mod == 0) { return 0; }
		else { return Mod / (2.0*OrigEdges); }
	}

	template<class PGraph>
	bool isExistEdge(const PGraph&Graph, int SrcNode, int desNode) {
		PGraph::TObj::TNodeI  NIi = Graph->GetNI(SrcNode);
		for (size_t e = 0; e < NIi.GetOutDeg(); e++)
			if (NIi.GetOutNId(e) == desNode)
				return true;
		return false;
	}

	template<class PGraph>
	double GetModularity2(const PGraph& Graph, TCnComV& CnComV) {
		TIntH OutDegH;
		const int OrigEdges = Graph->GetEdges();
		double Mod = 0;
		for (int c = 0; c < CnComV.Len(); c++) {
			const TIntV& NIdV = CnComV[c]();
			double AUV = 0;
			for (int i = 0; i < NIdV.Len() - 1; i++) {
				PGraph::TObj::TNodeI  NIi = Graph->GetNI(NIdV[i]);
				for (int j = i + 1; j < NIdV.Len(); j++) {
					PGraph::TObj::TNodeI  NIj = Graph->GetNI(NIdV[j]);
					bool IsexistEdgeB = isExistEdge(Graph,NIdV[i], NIdV[j]);
					if (IsexistEdgeB)
						AUV = 1;
					else
						AUV = 0;

					Mod += AUV - (NIi.GetOutDeg()*NIj.GetOutDeg() / (2.0*OrigEdges));
				}
			}
		}

		if (Mod == 0) { return 0; }
		else { return Mod / (2.0*OrigEdges); }
	}

}