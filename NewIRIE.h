#pragma once
#include "stdafx.h"
//#define NUM_LOOP_IR  20
#define NUM_SUBS_LOOP	3
#define EPSs 1e-10

class newIRIEclass{
	private:
		static TIntFltH dpV, AP;
		static TIntBoolH Visited;
		static int k;
		static bool saturated;

		static void NewIRIEGlobal(const PNGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);

		static int initialInfluence(const PNGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);

		static int residualInfluence(const PNGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);

		static int IRIEGetMax(const PNGraph& Graph, vector<int>& seed);

		static void computeAP(const PNGraph& Graph, const vector<int> &seeds, const int &src, const int& thetaa);

		static int findinArray(const int vertex[], int neighbourID, int n);


	public:
		static vector<int> callNewIRIEGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	};

