#pragma once
#include "stdafx.h"

namespace NSbaseLineRank {
	vector<int> callPageRank1Global(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callPageRank1Local(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS);
	vector<int> callPageRank1GlobalofLocal(const vector<int>& CandidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callPageRank2Global(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callPageRank2Local(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS);
	vector<int> callPageRank2GlobalofLocal(const vector<int>& CandidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callMaxDegreeRankGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callMaxDegreeRankGlobalofLocal(const vector<int>& vertexSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callMaxDegreeRankLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS);
	vector<int> callRandomCentrality(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callRandomCentrality_v2(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS);
	vector<int> callMinDegreeRankGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callMinDegreeRankGlobalofLocal(const vector<int>& vertexSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callMinDegreeRankLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS);
	vector<int> callKunduRankGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	vector<int> callKunduRankLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS);
	vector<int> callKunduRankGlobalofLocal(const vector<int>& CandidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int & indexPlot);
	
	
	template<class PGraph>
	void MinimumDegree(const PGraph& Graph, TIntIntH& MinimumRank) {
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			MinimumRank.AddDat(NI.GetId(), NI.GetOutDeg());
		}
		//MinimumRank.SortByDat(true);  //true == increasing  data value  Afzayeshii
	}

	template<class PGraph>
	void MaximimDegree(const PGraph& Graph, TIntIntH& MaximumRank) {
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			MaximumRank.AddDat(NI.GetId(), NI.GetOutDeg());
		}
		//MaximumRank.SortByDat(false);  //true == increasing  data value  Afzayeshii
		//MohsenTNT::PrintHaShTable(MaximumRank, Graph->GetNodes());
	}

	template<class PGraph>
	void KunduRank(const PGraph& Graph, TIntIntH& KunduRank) {
		//its ranking equal to 
		//Summmution of  node v out-degree + out-neighbor of node v out-degree 

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			KunduRank.AddDat(NI.GetId(), NI.GetOutDeg());
		}

		TIntIntH TempKunduRank;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int curID = NI.GetId();
			int sumOutDeg = KunduRank.GetDat(curID);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbourID = NI.GetOutNId(e);
				sumOutDeg += KunduRank.GetDat(neighbourID);
			}
			TempKunduRank.AddDat(curID, sumOutDeg);
		}
		KunduRank = TempKunduRank;
		//KunduRank.SortByDat(false);  //true == increasing  data value  Afzayeshii
									 //NStools::PrintHaShTable(KunduRank, Graph->GetNodes());
	}

	template<class PGraph>
	void MyPageRank1(const PGraph& Graph, TIntFltH& PRankH, const double& C, const double& Eps, const int& MaxIter) {
		const int NNodes = Graph->GetNodes();
		TVec<typename PGraph::TObj::TNodeI> NV;
		PRankH.Gen(NNodes);
		int MxId = -1;
		for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			NV.Add(NI);
			PRankH.AddDat(NI.GetId(), 1.0 / NNodes);
			int Id = NI.GetId();
			if (Id > MxId) {
				MxId = Id;
			}
		}

		TFltV PRankV(MxId + 1);
		TIntV OutDegV(MxId + 1);

		for (int j = 0; j < NNodes; j++) {
			typename PGraph::TObj::TNodeI NI = NV[j];
			int Id = NI.GetId();
			PRankV[Id] = 1.0 / NNodes;
			OutDegV[Id] = NI.GetOutDeg();
		}

		TFltV TmpV(NNodes);

		for (int iter = 0; iter < MaxIter; iter++) {
			for (int j = 0; j < NNodes; j++) {
				typename PGraph::TObj::TNodeI NI = NV[j];
				TFlt Tmp = 0;
				for (int e = 0; e < NI.GetInDeg(); e++) {
				const int InNId = NI.GetInNId(e);
				//for (int e = 0; e < NI.GetOutDeg(); e++) {
					//const int InNId = NI.GetOutNId(e);
					const int OutDeg = OutDegV[InNId];
					if (OutDeg > 0) {
						Tmp += PRankV[InNId] / OutDeg;
					}
				}
				TmpV[j] = C*Tmp; // Berkhin (the correct way of doing it)
			}
			double sum = 0;
			for (int i = 0; i < TmpV.Len(); i++) {
				sum += TmpV[i];
			}
			const double Leaked = (1.0 - sum) / double(NNodes);

			double diff = 0;
			for (int i = 0; i < NNodes; i++) {
				typename PGraph::TObj::TNodeI NI = NV[i];
				double NewVal = TmpV[i] + Leaked; // Berkhin
				int Id = NI.GetId();
				diff += fabs(NewVal - PRankV[Id]);
				PRankV[Id] = NewVal;
			}
			if (diff < Eps) {
				break;
			}
		}
		for (int i = 0; i < NNodes; i++) {
			typename PGraph::TObj::TNodeI NI = NV[i];
			PRankH[i] = PRankV[NI.GetId()];
		}
		//PRankH.SortByDat(false); //Decressing data value  //Asc
	 //MohsenTNT::PrintHaShTable(PRankH, Graph->GetNodes());
	}

	template<class PGraph>
	void MyPageRank2(const PGraph& Graph, TIntFltH& PRankH, const double& C , const double& Eps, const int& MaxIter) {
		const int NNodes = Graph->GetNodes();
		TVec<typename PGraph::TObj::TNodeI> NV;
		PRankH.Gen(NNodes);
		int MxId = -1;
		for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			NV.Add(NI);
			PRankH.AddDat(NI.GetId(), 1.0 / NNodes);
			int Id = NI.GetId();
			if (Id > MxId) {
				MxId = Id;
			}
		}

		TFltV PRankV(MxId + 1);
		TIntV OutDegV(MxId + 1);

		for (int j = 0; j < NNodes; j++) {
			typename PGraph::TObj::TNodeI NI = NV[j];
			int Id = NI.GetId();
			PRankV[Id] = 1.0 / NNodes;
			OutDegV[Id] = NI.GetOutDeg();
		}

		TFltV TmpV(NNodes);

		for (int iter = 0; iter < MaxIter; iter++) {
			for (int j = 0; j < NNodes; j++) {
				typename PGraph::TObj::TNodeI NI = NV[j];
				TFlt Tmp = 0;
				//for (int e = 0; e < NI.GetInDeg(); e++) {
				//const int InNId = NI.GetInNId(e);
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					const int InNId = NI.GetOutNId(e);
					const int OutDeg = OutDegV[InNId];
					if (OutDeg > 0) {
						Tmp += PRankV[InNId] / OutDeg;
					}
				}
				TmpV[j] = C*Tmp; // Berkhin (the correct way of doing it)
			}
			double sum = 0;
			for (int i = 0; i < TmpV.Len(); i++) {
				sum += TmpV[i];
			}
			const double Leaked = (1.0 - sum) / double(NNodes);

			double diff = 0;
			for (int i = 0; i < NNodes; i++) {
				typename PGraph::TObj::TNodeI NI = NV[i];
				double NewVal = TmpV[i] + Leaked; // Berkhin
				int Id = NI.GetId();
				diff += fabs(NewVal - PRankV[Id]);
				PRankV[Id] = NewVal;
			}
			if (diff < Eps) {
				break;
			}
		}
		for (int i = 0; i < NNodes; i++) {
			typename PGraph::TObj::TNodeI NI = NV[i];
			PRankH[i] = PRankV[NI.GetId()];
		}
		//PRankH.SortByDat(false); //Decressing data value  //Asc
								 //MohsenTNT::PrintHaShTable(PRankH, Graph->GetNodes());
	}
}
