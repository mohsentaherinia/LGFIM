#include "stdafx.h"

int NSIRIE::findinArray(const int vertex[], int neighbourID, int n)
{
	for (size_t i = 0; i < n; i++)
	{
		if (vertex[i] == neighbourID)
			return i;
	}
	return -1;
}

//vector<int> NSIRIE::callIRIEGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha , vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
//	vector<int> seedIRIE;
//	NSIRIE::IRIEGlobal(Graph, GC, seedIRIE, alpha, SeedSize,Model, ICProbb, MCS,curveInfoISV,curveInfoRTV,indexPlot);
//	NSplot::plotAllInfluenceSpread(curveInfoISV);
//	NSplot::plotAllInfluenceSpread(curveInfoRTV);
//	NStools::SaveSeedSetToFile("IRIEGlobal", GC, seedIRIE);
//	return seedIRIE;
//}
//
//vector<int> NSIRIE::callIRIEGlobalOfLocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
//	vector<int> seedIRIE;
//	NSIRIE::IRIEGlobalofLocal(Graph, GC, seedIRIE, alpha, SeedSize, Model, ICProbb, MCS, curveInfoISV, curveInfoRTV, indexPlot);
//	NSplot::plotAllInfluenceSpread(curveInfoISV);
//	NSplot::plotAllInfluenceSpread(curveInfoRTV);
//	NStools::SaveSeedSetToFile("IRIEGL", GC, seedIRIE);
//	return seedIRIE;
//}
//
//vector<int> NSIRIE::callIRIELocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha) {
//	vector<int> seedIRIE;
//	NSIRIE::IRIELocal(Graph, GC, seedIRIE, alpha, SeedSize, Model, ICProbb, MCS);
//	NStools::SaveSeedSetToFileLocal("IRIELocal", GC, seedIRIE);
//	return seedIRIE;
//}

vector<int> NSIRIE::callIRIEPMIAGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedIRIEPMIA;
	NSIRIE::IRIEPMIAGlobal(Graph, GC, seedIRIEPMIA, alpha, SeedSize, Model, ICProbb, MCS, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("IRIEPMIAGlobal", GC, seedIRIEPMIA);
	return seedIRIEPMIA;
}

vector<int> NSIRIE::callIRIEPMIAGlobalOfLocal(const vector<int>& candidateSet,const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedIRIEPMIA;
	NSIRIE::IRIEPMIAGlobalOfLocal(candidateSet,Graph, GC, seedIRIEPMIA, alpha, SeedSize, Model, ICProbb, MCS, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("IRIEPMIAGL", GC, seedIRIEPMIA);
	return seedIRIEPMIA;
}

vector<int> NSIRIE::callIRIEPMIALocal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha) {
	vector<int> seedIRIEPMIA;
	NSIRIE::IRIEPMIALocal(Graph, GC, seedIRIEPMIA, alpha, SeedSize, Model, ICProbb, MCS);
	NStools::SaveSeedSetToFileLocal("IRIEPMIALocal", GC, seedIRIEPMIA);
	return seedIRIEPMIA;
}



