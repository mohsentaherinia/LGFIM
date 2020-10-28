#pragma once
#include "stdafx.h"

namespace NSIS {
	int callInfluenceSpreadModel(const PNGraph & Graph, const vector<int>& seed, const string & Model, const double & ICProbb, const int & MCS);
}

namespace NStools {
	bool findNodeInVector(vector<int>& seed, const int& Node);
	void SaveToFile(char * algName, const int & k, const int & infectedNodes, const int & GraphSize, const int & mc, const char * Model, const double& Runtime, const GlobalConst & GC);
	//void SaveToFile(char * algName, const int & k, const int & infectedNodes, const int & GraphSize, const int & mc, const char * Model, const char * Runtime, const GlobalConst & GC);
	void SaveSeedSetToFile(char * algName, const GlobalConst & GC, const vector<int>& SeedSet);
	void SaveSeedSetToFileLocal(char * algName, const GlobalConst & GC, const vector<int>& SeedSet);
	
	template<class PGraph>
	void SparceTrippleEdgeData(const PGraph& Graph, TIntPrFltH &NodNodDat) {
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int curNodeiuID = NI.GetId();
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbourjvID = NI.GetOutNId(e);
				PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
				float EdgeProbe = 1.0 / NI2.GetInDeg();
				TIntPr x;
				x.Val1 = curNodeiuID;
				x.Val2 = neighbourjvID;
				NodNodDat.AddDat(x) = EdgeProbe;
			}
		}
		for (size_t i = 0; i < NodNodDat.Len(); i++) {
			//TIntPr x = NodNodDat.GetKey(i);
			//cout << "\n" << x.Val1() << "\t" << x.Val2() << "\t" << NodNodDat[i].Val;
			cout << "\n" << NodNodDat.GetKey(i).Val1() << "\t" << NodNodDat.GetKey(i).Val2() << "\t" << NodNodDat[i].Val;
		
		}
	}

	template <class MYDT>
	void PrintHaShTable(const MYDT& x, const int &size) {
		for (size_t i = 0; i < size; i++)
		{
			int _key = x.GetKey(i);
			float _data = x.GetDat(_key);
			cout << "\n" << i << " " << _key << " " << _data;
		}
		cout << "\n";
	}

	template <class MYDT>
	void PrintVector(vector<MYDT>& x) {
		int i = 0;
		for (vector<MYDT>::iterator it = x.begin(); it != x.end(); ++it)
			cout<<i++ <<" "<< *it<<endl;
		cout << endl;
	}

	template <class MYDT>
	void PrintVectorOfVector(vector<MYDT> x[], const int size) {
		for (size_t i = 0; i < size; i++){
			cout << i << "=";
			for (vector<MYDT>::iterator it = x[i].begin(); it != x[i].end(); ++it)
				cout << " " << *it;
			cout << endl;
		}
	}

	template <class MYDT>
	void PrintVectorOfVector(set<MYDT> x[], const int size) {
		for (size_t i = 0; i < size; i++) {
			cout << i << "=";
			for (set<MYDT>::iterator it = x[i].begin(); it != x[i].end(); ++it)
				cout << " " << *it;
			cout << endl;
		}
	}

}//end of namespace


namespace NSplot {
	struct curveInfo {
		TVec<TPair<TFlt, TFlt >> XY;
		char  algName[20];   //Greedy or PageRank or ...
		char  modelName[20]; //LT or IC or ...
		char typeCurve[20]; //IS or RT
		char  dataSet[40];  //dataset name
 	};

	void plotAllInfluenceSpread(vector<curveInfo> &curveInfo);

	void setPlottingValue(vector<curveInfo> &curveInfoX, const int & index, const int & xAxis, const float &yAxis, const char* algName, const char* modelName, const char* typeCurve);

	//template <class MYDT>
	//void plotInfluenceSpread(const MYDT& XY, const char  *algName, const char  *modelName) {
	//	char  Filename[80] = "_IS";
	//	strcat(Filename, algName);
	//	strcat(Filename, modelName);
	//	TGnuPlot Gp(Filename, modelName);
	//	Gp.AddPlot(XY, gpwLinesPoints, algName);
	//	//Gp.AddPlot(XY2, gpwPoints, "curve2");
	//	Gp.SetXYLabel("K", "InfluenceSpread");
	//	Gp.SavePng(); //or Gp.SaveEps();
	//				  ///open Picture in Windows////////////////////
	//	char fname[80];
	//	strcpy(fname, Filename);
	//	strcat(fname, ".png");
	//	system(fname);
	//}

	//template <class MYDT>
	//void plotRunningTime(const MYDT& XY, const char  *algName, const char  *modelName) {
	//	char  Filename[80] = "_RT";
	//	strcat(Filename, algName);
	//	strcat(Filename, modelName);
	//	TGnuPlot Gp(Filename, modelName);
	//	Gp.AddPlot(XY, gpwLinesPoints, algName);
	//	//Gp.AddPlot(XY2, gpwPoints, "curve2");
	//	Gp.SetXYLabel("K", "Running Time");
	//	Gp.SavePng(); //or Gp.SaveEps();
	//				  ///open Picture in Windows////////////////////
	//	char fname[80];
	//	strcpy(fname, Filename);
	//	strcat(fname, ".png");
	//	system(fname);
	//}
}