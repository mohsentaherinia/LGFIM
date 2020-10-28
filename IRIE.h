#pragma once
#include "stdafx.h"
#define NUM_LOOP_IR  20
#define EPss 1e-10


namespace NSIRIE {
	

	int findinArray(const int vertex[], int neighbourID, int n);

	/*vector<int> callIRIEGlobal(const PNGraph &Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callIRIEGlobalOfLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callIRIELocal(const PNGraph &Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha);
	*/
	vector<int> callIRIEPMIAGlobal(const PNGraph &Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callIRIEPMIAGlobalOfLocal(const vector<int> &candidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callIRIEPMIALocal(const PNGraph &Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & alpha);

	

	template<class PGraph>
	int IRIEGetMax(const PGraph& Graph, vector<int>& seed, const TIntFltH & dp)
	{
		int MaxIndexKey = -9999999;
		float MaxValue = -99999999;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int curNode = NI.GetId();
			if (dp.GetDat(curNode) > MaxValue)
				if (!NStools::findNodeInVector(seed, curNode)) {
					MaxIndexKey = curNode;
					MaxValue = dp.GetDat(curNode);
				}
		}
		return MaxIndexKey;
	}

//	template<class PGraph>
//	void IRIEGlobal(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize,const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot)
//	{
//		TExeTm ExeTmR;
//		double ElapsedSimTimes = 0, EndTime;
//		int NNodes = Graph->GetNodes();
//		bool changed = true;
//		bool saturated = false;
//		int i = 0;
//		double theta = 0.0001;
//		double edgeProb;
//		long count = 0;
//
//		TIntFltH dpV, newdpV, APs;
//		TIntBoolH Visited;
//
//		
//		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//			dpV.AddDat(NI.GetId(), 1.0);		//r(v)=1
//			newdpV.AddDat(NI.GetId(), 0.0);		//r(u)=1
//			APs.AddDat(NI.GetId(), 0.0);		//APs(u)=0
//			Visited.AddDat(NI.GetId(), false);	//V = V - u
//		}
//		//////////////////////////////////////////////////////////////////////////////////////////
//		for (size_t k = 0; k < SeedSize; k++) {
//			////////// For u isExist in S      APs(u)=1
//			for (size_t index2 = 0; index2 < seed.size() ; index2++)
//				APs.AddDat(seed[index2]) = 0;
//			//////////
//			////////// For u isExist in V - S   APs(u)=SetValue
//			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//				int curNodet = NI.GetId();
//				if (!Visited.GetDat(curNodet))
//					APs.AddDat(curNodet) = 0;//???????????????????? 1 equal to IR   //  0  BETTER THAN ir
//			}
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//				int curNodet = NI.GetId();
//				dpV.AddDat(curNodet) = (1 - APs.GetDat(curNodet))*dpV.GetDat(curNodet);
//			}
//			////////////////////////////////////////////////////////////////////////////////////////////////////////
//			while (count++ < NUM_LOOP_IR && changed && !saturated)
//			{
//				changed = false;
//
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					newdpV.AddDat(NI.GetId(), 0.0);
//				}
//				//////////////////////////////////////////
//				float value;
//				int curNodeiuID, neighbourjvID;
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					curNodeiuID = NI.GetId();
//					for (int e = 0; e < NI.GetOutDeg(); e++) {
//						neighbourjvID = NI.GetOutNId(e);
//						PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
//						edgeProb = (float)1.0 / NI2.GetInDeg();
//						if (edgeProb)
//						{
//							//newdp[i] += alpha * edgeProb * dp[e.v];
//							value = newdpV.GetDat(curNodeiuID);
//							value += (1 - APs.GetDat(curNodeiuID)) *alpha * edgeProb * dpV.GetDat(neighbourjvID);
//							newdpV.AddDat(curNodeiuID) = value;
//						}
//
//					}
//					//newdp[i] += 1;
//					//newnewdp[i] += 1 - ap[i];
//					newdpV.AddDat(curNodeiuID) = newdpV.GetDat(curNodeiuID) + 1 - APs.GetDat(curNodeiuID);
//				}
//				/////////////////////////////////////////////////////////////////////////////////////
//				int curNodei;
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					curNodei = NI.GetId();
//					////if((newdp[i] < dp[i] - theta) || (newdp[i] > dp[i] + theta))
//					if ((newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) < -theta) ||
//						(newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) > theta))
//						changed = true;
//					////if(newdp[i] > NNodes)
//					if (newdpV.GetDat(curNodei) > NNodes)
//						saturated = true;
//					////dp[i] = newdp[i];
//					dpV.AddDat(curNodei) = newdpV.GetDat(curNodei);
//				}
//			}//end inner Loop
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			int u = NSIRIE::IRIEGetMax(Graph, seed, dpV);//U = Argmax dv(u)
//			seed.push_back(u);//S = S + u
//			Visited.AddDat(u) =true ;//V = V - u
//			//////////
//			EndTime = ExeTmR.GetSecs();
//			double exactRunTime = EndTime - ElapsedSimTimes;
//			int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
//			cout << endl << "IRIEGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
//			NStools::SaveToFile("IRIEGlobal", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
//			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "IRIEGlobal", "InfluenceSpraed");
//			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "IRIEGlobal", "RunningTime");
//			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
//		/////////
//		}//end outer Loop
//}
//
//	template<class PGraph>
//	void IRIEGlobalofLocal(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
//		TExeTm ExeTmR;
//		double ElapsedSimTimes = 0, EndTime;
//		int NNodes = Graph->GetNodes();
//		bool changed = true;
//		bool saturated = false;
//		int i = 0;
//		double theta = 0.0001;
//		double edgeProb;
//		long count = 0;
//
//		TIntFltH dpV, newdpV, APs;
//		TIntBoolH Visited;
//
//
//		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//			dpV.AddDat(NI.GetId(), 1.0);		//r(v)=1
//			newdpV.AddDat(NI.GetId(), 0.0);		//r(u)=1
//			APs.AddDat(NI.GetId(), 0.0);		//APs(u)=0
//			Visited.AddDat(NI.GetId(), false);	//V = V - u
//		}
//		//////////////////////////////////////////////////////////////////////////////////////////
//		for (size_t k = 0; k < 1; k++) {
//			////////// For u isExist in S      APs(u)=1
//			for (size_t index2 = 0; index2 < seed.size(); index2++)
//				APs.AddDat(seed[index2]) = 0;
//			//////////
//			////////// For u isExist in V - S   APs(u)=SetValue
//			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//				int curNodet = NI.GetId();
//				if (!Visited.GetDat(curNodet))
//					APs.AddDat(curNodet) = 0;//???????????????????? 1 equal to IR   //  0  BETTER THAN ir
//			}
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//				int curNodet = NI.GetId();
//				dpV.AddDat(curNodet) = (1 - APs.GetDat(curNodet))*dpV.GetDat(curNodet);
//			}
//			////////////////////////////////////////////////////////////////////////////////////////////////////////
//			while (count++ < NUM_LOOP_IR && changed && !saturated)
//			{
//				changed = false;
//
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					newdpV.AddDat(NI.GetId(), 0.0);
//				}
//				//////////////////////////////////////////
//				float value;
//				int curNodeiuID, neighbourjvID;
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					curNodeiuID = NI.GetId();
//					for (int e = 0; e < NI.GetOutDeg(); e++) {
//						neighbourjvID = NI.GetOutNId(e);
//						PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
//						edgeProb = (float)1.0 / NI2.GetInDeg();
//						if (edgeProb)
//						{
//							//newdp[i] += alpha * edgeProb * dp[e.v];
//							value = newdpV.GetDat(curNodeiuID);
//							value += (1 - APs.GetDat(curNodeiuID)) *alpha * edgeProb * dpV.GetDat(neighbourjvID);
//							newdpV.AddDat(curNodeiuID) = value;
//						}
//
//					}
//					//newdp[i] += 1;
//					//newnewdp[i] += 1 - ap[i];
//					newdpV.AddDat(curNodeiuID) = newdpV.GetDat(curNodeiuID) + 1 - APs.GetDat(curNodeiuID);
//				}
//				/////////////////////////////////////////////////////////////////////////////////////
//				int curNodei;
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					curNodei = NI.GetId();
//					////if((newdp[i] < dp[i] - theta) || (newdp[i] > dp[i] + theta))
//					if ((newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) < -theta) ||
//						(newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) > theta))
//						changed = true;
//					////if(newdp[i] > NNodes)
//					if (newdpV.GetDat(curNodei) > NNodes)
//						saturated = true;
//					////dp[i] = newdp[i];
//					dpV.AddDat(curNodei) = newdpV.GetDat(curNodei);
//				}
//			}//end inner Loop
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//		}//end outer Loop
//		dpV.SortByDat(false);
//		for (size_t k = 0; k < SeedSize;k++) {
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			seed.push_back(dpV.GetKey(k));
//			EndTime = ExeTmR.GetSecs();
//			double exactRunTime = EndTime - ElapsedSimTimes;
//			int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
//			cout << endl << "IRIEGL   K=" << k+1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
//			NStools::SaveToFile("IRIEGL", k+1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime, GC);
//			NSplot::setPlottingValue(curveInfoISV, indexPlot, k+1, infectedNodes, Model, "IRIEGL", "InfluenceSpraed");
//			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k+1, exactRunTime, Model, "IRIEGL", "RunningTime");
//			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
//			/////////
//		}//end outer Loop
//	}
//
//	template<class PGraph>
//	void IRIELocal(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS)
//	{
//		TExeTm ExeTmR;
//		double ElapsedSimTimes = 0, EndTime;
//
//		int NNodes = Graph->GetNodes();
//		bool changed = true;
//		bool saturated = false;
//		int i = 0;
//		double theta = 0.0001;
//		double edgeProb;
//		long count = 0;
//
//		TIntFltH dpV, newdpV, APs;
//		TIntBoolH Visited;
//
//
//		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//			dpV.AddDat(NI.GetId(), 1.0);		//r(v)=1
//			newdpV.AddDat(NI.GetId(), 0.0);		//r(u)=1
//			APs.AddDat(NI.GetId(), 0.0);		//APs(u)=0
//			Visited.AddDat(NI.GetId(), false);	//V = V - u
//		}
//		//////////////////////////////////////////////////////////////////////////////////////////
//		for (size_t k = 0; k < SeedSize; k++) {
//			////////// For u isExist in S      APs(u)=1
//			for (size_t index2 = 0; index2 < seed.size(); index2++)
//				APs.AddDat(seed[index2]) = 0;
//			//////////
//			////////// For u isExist in V - S   APs(u)=SetValue
//			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//				int curNodet = NI.GetId();
//				if (!Visited.GetDat(curNodet))
//					APs.AddDat(curNodet) = 0;//???????????????????? 1 equal to IR   //  0  BETTER THAN ir
//			}
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			/////////////////////////////////////////////////////////////////////////////////////////////////////////
//			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//				int curNodet = NI.GetId();
//				dpV.AddDat(curNodet) = (1 - APs.GetDat(curNodet))*dpV.GetDat(curNodet);
//			}
//			////////////////////////////////////////////////////////////////////////////////////////////////////////
//			while (count++ < NUM_LOOP_IR && changed && !saturated)
//			{
//				changed = false;
//
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					newdpV.AddDat(NI.GetId(), 0.0);
//				}
//				//////////////////////////////////////////
//				float value;
//				int curNodeiuID, neighbourjvID;
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					curNodeiuID = NI.GetId();
//					for (int e = 0; e < NI.GetOutDeg(); e++) {
//						neighbourjvID = NI.GetOutNId(e);
//						PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
//						edgeProb = (float)1.0 / NI2.GetInDeg();
//						if (edgeProb)
//						{
//							//newdp[i] += alpha * edgeProb * dp[e.v];
//							value = newdpV.GetDat(curNodeiuID);
//							value += (1 - APs.GetDat(curNodeiuID)) *alpha * edgeProb * dpV.GetDat(neighbourjvID);
//							newdpV.AddDat(curNodeiuID) = value;
//						}
//
//					}
//					//newdp[i] += 1;
//					//newnewdp[i] += 1 - ap[i];
//					newdpV.AddDat(curNodeiuID) = newdpV.GetDat(curNodeiuID) + 1 - APs.GetDat(curNodeiuID);
//				}
//				/////////////////////////////////////////////////////////////////////////////////////
//				int curNodei;
//				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//					curNodei = NI.GetId();
//					////if((newdp[i] < dp[i] - theta) || (newdp[i] > dp[i] + theta))
//					if ((newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) < -theta) ||
//						(newdpV.GetDat(curNodei) - dpV.GetDat(curNodei) > theta))
//						changed = true;
//					////if(newdp[i] > NNodes)
//					if (newdpV.GetDat(curNodei) > NNodes)
//						saturated = true;
//					////dp[i] = newdp[i];
//					dpV.AddDat(curNodei) = newdpV.GetDat(curNodei);
//				}
//			}//end inner Loop
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			 /////////////////////////////////////////////////////////////////////////////////////////////////////////
//			int u = NSIRIE::IRIEGetMax(Graph, seed, dpV);//U = Argmax dv(u)
//			seed.push_back(u);//S = S + u
//			Visited.AddDat(u) = true;//V = V - u
//			EndTime = ExeTmR.GetSecs();
//			double exactRunTime = EndTime - ElapsedSimTimes;
//			int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
//			NStools::SaveToFile("IRIELocal", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
//			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
//			////////
//		}//end outer Loop
//	}
//
	
	template<class PGraph>
	void IRIEPMIAGlobal(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		//int NNodes = Graph->GetNodes();
		double edgeProb;
		long count = 0;

		TIntFltH dpV, APs;

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			dpV.AddDat(NI.GetId(), 1.0);		//r(v)=1
			APs.AddDat(NI.GetId(), 0.0);		//APs(u)=0
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		for (size_t k = 0; k < SeedSize; k++) {
			//count = 0;
			////////////////////////////////////////////////////////////////////////////////////////////////////////
			while (count++ < NUM_LOOP_IR)
			{
				
				int curNodeiuID, neighbourjvID;
				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
					float value = 0;
					curNodeiuID = NI.GetId();
					for (int e = 0; e < NI.GetOutDeg(); e++) {
						neighbourjvID = NI.GetOutNId(e);
						PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
						edgeProb = (float)1.0 / NI2.GetInDeg();
						if (edgeProb)
							value += edgeProb * dpV.GetDat(neighbourjvID);
					}
					dpV.AddDat(curNodeiuID) =  (1 - APs.GetDat(curNodeiuID)) *(1 + alpha * value);
				}
			}//end inner Loop
			int u = NSIRIE::IRIEGetMax(Graph, seed, dpV);//U = Argmax dv(u)
			seed.push_back(u);//S = S + u
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			int infectedNodes = 0;
			if (GC._simMode == 0)
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
			cout << endl << "IRIEPMIAGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
			NStools::SaveToFile("IRIEPMIAGlobal", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "IRIE", "Influence Spraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "IRIE", "Running Time (Sec)");
			////////
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
			/*computeAP(Graph, seed, u, 320, APs);
			for (size_t i = 0; i < seed.size(); i++)
				APs.AddDat(seed[i]) = 1;*/
		}
	}

	template<class PGraph>
	void IRIEPMIAGlobalOfLocal(const vector<int>& CandidateSet,const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		//int NNodes = Graph->GetNodes();
		double edgeProb;
		long count = 0;

		TIntFltH dpV, APs;

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			dpV.AddDat(NI.GetId(), 1.0);		//r(v)=1
			APs.AddDat(NI.GetId(), 0.0);		//APs(u)=0
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		while (count++ < 3)
		{
			int curNodeiuID, neighbourjvID;
			for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				float value = 0;
				curNodeiuID = NI.GetId();
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					neighbourjvID = NI.GetOutNId(e);
					PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
					edgeProb = (float)1.0 / NI2.GetInDeg();
					if (edgeProb)
						value += edgeProb * dpV.GetDat(neighbourjvID);
				}
				dpV.AddDat(curNodeiuID) = (1 - APs.GetDat(curNodeiuID)) *(1 + alpha * value);
			}
		}//end inner Loop
		while (count++ < NUM_LOOP_IR)
		{
			int curNodeiuID, neighbourjvID;
			//for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
				//curNodeiuID = NI.GetId();
			for (size_t i = 0; i < CandidateSet.size(); i++) {
				curNodeiuID = CandidateSet[i];
				PGraph::TObj::TNodeI NI = Graph->GetNI(curNodeiuID);
				float value = 0;
				for (int e = 0; e < NI.GetOutDeg(); e++) {
					neighbourjvID = NI.GetOutNId(e);
					PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
					edgeProb = (float)1.0 / NI2.GetInDeg();
					if (edgeProb)
						value += edgeProb * dpV.GetDat(neighbourjvID);
				}
				dpV.AddDat(curNodeiuID) = (1 - APs.GetDat(curNodeiuID)) *(1 + alpha * value);
			}
		}//end inner Loop

		dpV.SortByDat(false);
		/*TIntFltH dpVCandidateSet;
		for (size_t i = 0; i < CandidateSet.size(); i++) {
			dpVCandidateSet.AddDat(CandidateSet[i]) = dpV.GetDat(CandidateSet[i]);
		}
		dpVCandidateSet.SortByDat(false);
*/
		for (size_t k = 0; k < SeedSize; k++) {
			////////////////////////////////////////////////////////////////////////////////////////////////////////
			//seed.push_back(dpVCandidateSet.GetKey(k));
			seed.push_back(dpV.GetKey(k));
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			int infectedNodes = 0;
			if (GC._simMode == 0)
				infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
			cout << endl << "IRIEPMIAGL   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
			NStools::SaveToFile("IRIEPMIAGL", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime, GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "LGFIM", "Influence Spraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "LGFIM", "Running Time (Sec)");
			////////
			FILE *FRT = fopen("__Compare_RT_IRIEGLwithCELF.not", "at+");
			fprintf(FRT, "%Lf,\t", ExeTmR.GetSecs());
			fclose(FRT);
			////////
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
		}//end outer Loop
	}

	template<class PGraph>
	void IRIEPMIALocal(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		//int NNodes = Graph->GetNodes();
		double edgeProb;
		long count = 0;

		TIntFltH dpV, APs;

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			dpV.AddDat(NI.GetId(), 1.0);		//r(v)=1
			APs.AddDat(NI.GetId(), 0.0);		//APs(u)=0
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		for (size_t k = 0; k < SeedSize; k++) {
			//count = 0;
			////////////////////////////////////////////////////////////////////////////////////////////////////////
			while (count++ < NUM_LOOP_IR)
			{

				int curNodeiuID, neighbourjvID;
				for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
					float value = 0;
					curNodeiuID = NI.GetId();
					for (int e = 0; e < NI.GetOutDeg(); e++) {
						neighbourjvID = NI.GetOutNId(e);
						PGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
						edgeProb = (float)1.0 / NI2.GetInDeg();
						if (edgeProb)
							value += edgeProb * dpV.GetDat(neighbourjvID);
					}
					dpV.AddDat(curNodeiuID) = (1 - APs.GetDat(curNodeiuID)) *(1 + alpha * value);
				}
			}//end inner Loop
			int u = NSIRIE::IRIEGetMax(Graph, seed, dpV);//U = Argmax dv(u)
			seed.push_back(u);//S = S + u
			EndTime = ExeTmR.GetSecs();
			double exactRunTime = EndTime - ElapsedSimTimes;
			int infectedNodes = 999;// NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
			//cout << endl << "IRIEPMIAGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
			NStools::SaveToFile("IRIEPMIALocal", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime, GC);
			//NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "IRIEPMIAGlobal", "InfluenceSpraed");
			//NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "IRIEPMIAGlobal", "RunningTime");
			////////
			ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
			/*computeAP(Graph, seed, u, 320, APs);
			for (size_t i = 0; i < seed.size(); i++)
			APs.AddDat(seed[i]) = 1;*/
		}
	}

	template<class PGraph>
	void computeAP(const PGraph& Graph, const vector<int> &seeds, const int &src, const int& theta, TIntFltH& APs)
	{
		const int n = Graph->GetNodes();
		int* vertex = new int[n];
		bool* Visited = new bool[n];
		Heap h;
		HeapNode hNode;
		HeapNode curNode;
		int minNode;
		double minDist;
		double threshold = log((double)theta);

		initHeap(&h, n + 1);

		int z = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			Visited[z] = false;
			vertex[z] = NI.GetId();
			z++;
		}
		for (size_t index2 = 0; index2 < seeds.size() - 1 && seeds[index2] != NULL; index2++)
			Visited[findinArray(vertex, seeds[index2], n)] = true;
		for (int i = 0; i < n; i++) {
			if (!Visited[i] && vertex[i] != src) {
				hNode.key = i;
				hNode.value = DBL_MAX;
				insertHeap(&h, hNode);
			}
		}

		hNode.key = findinArray(vertex, src, n);
		hNode.value = 0;
		insertHeap(&h, hNode);

		while (h.count)
		{
			curNode = h.elements[1];
			minNode = curNode.key;
			minDist = curNode.value;

			if (minDist >= threshold)
				break;
			PNGraph::TObj::TNodeI NI2 = Graph->GetNI(vertex[minNode]);
			for (int e = 0; e < NI2.GetOutDeg(); e++) {
				int neighbourID = NI2.GetOutNId(e);
				int L = findinArray(vertex, neighbourID, n);
				if (!Visited[L])
				{
					PGraph::TObj::TNodeI NI3 = Graph->GetNI(neighbourID);
					double edgeProb = (float)1.0 / NI3.GetInDeg();

					if (h.elements[h.index[L]].value > minDist + edgeProb + EPss)
						decreaseKeyHeap(&h, h.index[L], minDist + edgeProb);
				}
			}

			APs.AddDat(vertex[minNode]) = APs.GetDat(vertex[minNode]) + exp(-minDist);
			if (APs.GetDat(vertex[minNode]) > 1 - EPss)
				APs.AddDat(vertex[minNode]) = 1;
			removeMinHeap(&h);
		}
		freeHeap(&h);
	}


}

