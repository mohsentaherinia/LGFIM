#include "stdafx.h"
vector<int>  NSCOFIM::callCofimGreedyGlobal(const PNGraph *SubG,const int& CmtyNum, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	TIntIntH NodeCmtyH;
	for (size_t i = 0; i < CmtyNum; i++)
		for (PNGraph::TObj::TNodeI NI = SubG[i]->BegNI(); NI < SubG[i]->EndNI(); NI++)
			NodeCmtyH.AddDat(NI.GetId(), i);
	
	vector<int> seedCOFIM;
	NSCOFIM::CofimGreedyGlobal(NodeCmtyH,Graph, GC, Model, seedCOFIM, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("CofimGreedyGlobal", GC, seedCOFIM);
	return seedCOFIM;
}

vector<int>  NSCOFIM::callCofimCELFGlobal(const PNGraph *SubG, const int& CmtyNum, const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	TIntIntH NodeCmtyH;
	for (size_t i = 0; i < CmtyNum; i++)
		for (PNGraph::TObj::TNodeI NI = SubG[i]->BegNI(); NI < SubG[i]->EndNI(); NI++)
			NodeCmtyH.AddDat(NI.GetId(), i);

	vector<int> seedCOFIM;
	NSCOFIM::CofimCELFGlobal(NodeCmtyH, Graph, GC, Model, seedCOFIM, SeedSize, ICProb, MC, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("CofimCELFGlobal", GC, seedCOFIM);
	return seedCOFIM;
}

long NSCOFIM::CountOf_N_S(const PNGraph &Graph, const vector<int>& seed) {
	long seedOutDegree = 0;
	for (size_t i = 0; i < seed.size(); i++) {
		PNGraph::TObj::TNodeI NI = Graph->GetNI(seed[i]);
		seedOutDegree += NI.GetOutDeg();
	}
	return seedOutDegree;
}

inline long NSCOFIM::CountOf_N_CurNode(const PNGraph &Graph, const int& vertex) {
	long seedOutDegree = 0;
	PNGraph::TObj::TNodeI NI = Graph->GetNI(vertex);
	seedOutDegree += NI.GetOutDeg();
	return seedOutDegree;
}

long NSCOFIM::CountOf_NC_S(const PNGraph &Graph, const vector<int>& seed, const TIntIntH& NodeCmtyH, set<int>& nbrsCmtyList) {
	//pass bt reference
	for (size_t i = 0; i < seed.size(); i++) {
		PNGraph::TObj::TNodeI NI = Graph->GetNI(seed[i]);
		for (int e = 0; e < NI.GetOutDeg(); e++) {
			 nbrsCmtyList.insert(NodeCmtyH.GetDat(NI.GetOutNId(e)));
		}
	}
	return nbrsCmtyList.size();
}

inline long NSCOFIM::CountOf_NC_CurNode(const PNGraph &Graph, const int& curNode, const TIntIntH& NodeCmtyH,const set<int>& nbrsCmtyList) {
	//Pass By Value
	set<int> CurNodeNbrsCmtyList;
	PNGraph::TObj::TNodeI NI = Graph->GetNI(curNode);
	for (int e = 0; e < NI.GetOutDeg(); e++) {
		if (nbrsCmtyList.find(NodeCmtyH.GetDat(NI.GetOutNId(e))) == nbrsCmtyList.end())
			CurNodeNbrsCmtyList.insert(NodeCmtyH.GetDat(NI.GetOutNId(e)));
	}
	return CurNodeNbrsCmtyList.size();
}



