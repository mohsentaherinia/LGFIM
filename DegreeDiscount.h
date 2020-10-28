#pragma once
#include "stdafx.h"
namespace NSDegreeDiscount {
	vector<int> callDegDisGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, const double & ratio, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callDegDisGlobalOfLocal(const vector<int> &candidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, const double & ratio, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callDegDisLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, const double & ratio);

	template<class PGraph>
	void DegDisGlobal(const PGraph& Graph, const GlobalConst &GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC,const double& ratio, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		long n = Graph->GetNodes();
		TIntBoolH used;
		TIntIntH countV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			countV.AddDat(NI.GetId(), 0.0);
			used.AddDat(NI.GetId(), false);
		}
	
		for (int k = 0; k < SeedSize; k++)
		{
			double max = -1000000.0;
			int CanNode = -1;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				int curVer = NI.GetId();
				if (!used.GetDat(curVer))
				{
					double tmp = NI.GetOutDeg() - 2 * countV.GetDat(curVer) - ratio*countV.GetDat(curVer) * (NI.GetOutDeg() - countV.GetDat(curVer));
						if (tmp > max)
						{
							max = tmp;
							CanNode = curVer;
						}
				}
			}
			used.AddDat(CanNode) = true;
			seeds.push_back(CanNode);
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			int infectedNodes = 999;
			if (GC._simMode == 0)
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
			cout << endl << "DegDisGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
			NStools::SaveToFile("DegDisGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, exactRunTime,GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "Deg-Dis", "Influence Spraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "Deg-Dis", "Running Time (Sec)");
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
			PGraph::TObj::TNodeI NI2 = Graph->GetNI(CanNode);
			for (int e = 0; e < NI2.GetOutDeg(); e++) {
				int neighbourCanNodeID = NI2.GetOutNId(e);
				countV.AddDat(neighbourCanNodeID) = countV.GetDat(neighbourCanNodeID) +  1;
			}
		}
	}

	template<class PGraph>
	void DegDisGlobalOfLocal(const vector<int> &candidateSet, const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, const double& ratio, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime; 
		long n = Graph->GetNodes();
		TIntBoolH used;
		TIntIntH countV;
		for (size_t i = 0; i < candidateSet.size(); i++) {
			countV.AddDat(candidateSet[i], 0.0);
			used.AddDat(candidateSet[i], false);
		}
		
		for (int k = 0; k < SeedSize; k++)
		{
			double max = -1000000.0;
			int CanNode = -1;
			for (size_t i = 0; i < candidateSet.size(); i++){
				int curVer = candidateSet[i];
				PGraph::TObj::TNodeI NI = Graph->GetNI(curVer);
				if (!used.GetDat(curVer))
				{
					double tmp = NI.GetOutDeg() - 2 * countV.GetDat(curVer) - ratio*countV.GetDat(curVer) * (NI.GetOutDeg() - countV.GetDat(curVer));
					if (tmp > max)
					{
						max = tmp;
						CanNode = curVer;
					}
				}
			}
			used.AddDat(CanNode) = true;
			seeds.push_back(CanNode);
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
			cout << endl << "DegDisGL   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
			NStools::SaveToFile("DegDisGL", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, exactRunTime,GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "DegDisGL", "InfluenceSpraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "DegDisGL", "RunningTime");
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime; 
			PGraph::TObj::TNodeI NI2 = Graph->GetNI(CanNode);
			for (int e = 0; e < NI2.GetOutDeg(); e++) {
				int neighbourCanNodeID = NI2.GetOutNId(e);
				if (countV.IsKey(neighbourCanNodeID))
					countV.AddDat(neighbourCanNodeID) = countV.GetDat(neighbourCanNodeID) + 1;
			}
		}
	}

	template<class PGraph>
	void DegDisLocal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, const double& ratio) {
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime; 
		long n = Graph->GetNodes();
		TIntBoolH used;
		TIntIntH countV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			countV.AddDat(NI.GetId(), 0.0);
			used.AddDat(NI.GetId(), false);
		}
		
		for (int k = 0; k < SeedSize; k++)
		{
			double max = -1000000.0;
			int CanNode = -1;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				int curVer = NI.GetId();
				if (!used.GetDat(curVer))
					{
					double tmp = NI.GetOutDeg() - 2 * countV.GetDat(curVer) - ratio*countV.GetDat(curVer) * (NI.GetOutDeg() - countV.GetDat(curVer));
					if (tmp > max)
					{
						max = tmp;
						CanNode = curVer;
					}
				}
			}
			used.AddDat(CanNode) = true;
			seeds.push_back(CanNode);
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
			NStools::SaveToFile("DegDisLocal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, exactRunTime, GC);
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
			PGraph::TObj::TNodeI NI2 = Graph->GetNI(CanNode);
			for (int e = 0; e < NI2.GetOutDeg(); e++) {
				int neighbourCanNodeID = NI2.GetOutNId(e);
				countV.AddDat(neighbourCanNodeID) = countV.GetDat(neighbourCanNodeID) + 1;
			}
		}
	}


}