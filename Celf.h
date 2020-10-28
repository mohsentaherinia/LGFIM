#pragma once
#include "stdafx.h"


namespace NScelf {
	vector<int> callCELFGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callCELFGlobalOfLocal(const vector<int>& candidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callCELFLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC);
	vector<int> callCELFPPPGlobal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callCELFPPPGlobalOfLocal(const vector<int>& candidateSet, const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callCELFPPPLocal(const PNGraph & Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProb, const int & MC);

	void QueueSortFirstK(queue<TIntIntH>& q, const int & N);
	void QueueSortOthersK(queue<TIntIntH>& q, const int & N);
	void printqueue(queue<TIntIntH> q);

	template<class PGraph>
	void CELFGlobal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTm1;
		queue<TIntIntH> q;
		vector <int> seedT1(1, -1111);
		int infectedNodes, curVer;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			TIntIntH item1;
			curVer = NI.GetId();
			seedT1[0] = curVer;
			infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
			item1.AddDat(curVer) = infectedNodes;
			q.push(item1);
		}
		NScelf::QueueSortFirstK(q, q.size());
		//# Select the first node and remove from candidate list
		TIntIntH item2 = q.front();
		int k = 0;
		seeds.push_back(item2.GetKey(0));
		infectedNodes = item2.GetDat(seeds[k]);
		q.pop();
		cout << endl << "CELFGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
		NStools::SaveToFile("CELFGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELF", "Influence Spraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELF", "Runnin gTime (SEC)");
		//	# Find the next k - 1 nodes using the list - sorting procedure
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int infectedNodesLocal;
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = q.front();
				int current = item3.GetKey(0);
				//Q[0] = (current, IC(g, S + [current], p, mc) - spread)
				//remove front of q (current) => then Evaluate the IC(s+current) => push cur and its influence to q
				q.pop();
				seeds.push_back(current);
				infectedNodesLocal = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				TIntIntH item4;
				item4.AddDat(current) = infectedNodesLocal - infectedNodes;
				q.push(item4);
				NScelf::QueueSortOthersK(q, q.size());
				TIntIntH item5 = q.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = q.front();
			seeds.push_back(item6.GetKey(0));
			infectedNodes += item6.GetDat(seeds[k]);
			q.pop();
			cout << endl << "CELFGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
			NStools::SaveToFile("CELFGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELF", "Influence Spraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELF", "Runnin gTime (SEC)");
		}//end for
	}

	template<class PGraph>
	void CELFLocal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC) {
		TExeTm ExeTm1;
		queue<TIntIntH> q;
		vector <int> seedT1(1, -1111);
		int infectedNodes, curVer;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			TIntIntH item1;
			curVer = NI.GetId();
			seedT1[0] = curVer;
			infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
			item1.AddDat(curVer) = infectedNodes;
			q.push(item1);
		}
		NScelf::QueueSortFirstK(q, q.size());
		//# Select the first node and remove from candidate list
		TIntIntH item2 = q.front();
		int k = 0;
		seeds.push_back(item2.GetKey(0));
		infectedNodes = item2.GetDat(seeds[k]);
		q.pop();
		NStools::SaveToFile("CELFLocal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		//	# Find the next k - 1 nodes using the list - sorting procedure
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int infectedNodesLocal;
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = q.front();
				int current = item3.GetKey(0);
				//Q[0] = (current, IC(g, S + [current], p, mc) - spread)
				//remove front of q (current) => then Evaluate the IC(s+current) => push cur and its influence to q
				q.pop();
				seeds.push_back(current);
				infectedNodesLocal = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back(); 
				TIntIntH item4;
				item4.AddDat(current) = infectedNodesLocal - infectedNodes;
				q.push(item4);
				NScelf::QueueSortOthersK(q, q.size());
				TIntIntH item5 = q.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = q.front();
			seeds.push_back(item6.GetKey(0));
			infectedNodes += item6.GetDat(seeds[k]);
			q.pop();
			NStools::SaveToFile("CELFLocal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
			}//end for
	}

	template<class PGraph>
	void CELFGlobalOfLocal(const vector<int>& candidateSet, const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTm1;
		queue<TIntIntH> q;
		vector <int> seedT1(1, -1111);
		int infectedNodes, curVer;
		for (size_t i = 0; i < candidateSet.size(); i++)
		{
			TIntIntH item1;
			curVer = candidateSet[i];
			seedT1[0] = curVer;
			infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
			item1.AddDat(curVer) = infectedNodes;
			q.push(item1);
		}
		/*for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		TIntIntH item1;
		curVer = NI.GetId();
		seedT1[0] = curVer;
		infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
		item1.AddDat(curVer) = infectedNodes;
		q.push(item1);
		}*/
		NScelf::QueueSortFirstK(q, q.size());
		//# Select the first node and remove from candidate list
		TIntIntH item2 = q.front();
		int k = 0;
		seeds.push_back(item2.GetKey(0));
		infectedNodes = item2.GetDat(seeds[k]);
		q.pop();
		cout << endl << "CELFGL   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
		NStools::SaveToFile("CELFGL", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELFGL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELFGL", "RunningTime");
		//	# Find the next k - 1 nodes using the list - sorting procedure
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int infectedNodesLocal;
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = q.front();
				int current = item3.GetKey(0);
				//Q[0] = (current, IC(g, S + [current], p, mc) - spread)
				//remove front of q (current) => then Evaluate the IC(s+current) => push cur and its influence to q
				q.pop();
				seeds.push_back(current);
				infectedNodesLocal = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back(); 
				TIntIntH item4;
				item4.AddDat(current) = infectedNodesLocal - infectedNodes;
				q.push(item4);
				NScelf::QueueSortOthersK(q, q.size());
				TIntIntH item5 = q.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = q.front();
			seeds.push_back(item6.GetKey(0));
			infectedNodes += item6.GetDat(seeds[k]);
			q.pop();
			cout << endl << "CELFGL   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
			NStools::SaveToFile("CELFGL", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELFGL", "InfluenceSpraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELFGL", "RunningTime");

		}//end for
	}

	template<class PGraph>
	void CELFPPPGlobal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTm1;
		queue<TIntIntH> qCELF;
		queue<int> qParent;
		queue<int> qParentTemp;
		TIntBoolH Visited;
		TIntBoolH VisitedTemp;
		TIntFltH spread;
		int infectedNodes;

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			//spread.AddDat(NI.GetId()) = 0.0;
			spread.AddDat(NI.GetId()) = 1 + NI.GetOutDeg();
			Visited.AddDat(NI.GetId()) = false;
		}
		//NStools::PrintHaShTable(spread, spread.Len());
		vector<int> vertexCoverSet;
		NSsimPath::get_vertex_cover2(Graph, vertexCoverSet);

		//Calculate Influence Spread for coverVertex set
		vector <int> seedT1(1, -1111);
		for (size_t j = 0; j < vertexCoverSet.size(); j++) {
			int curNode = vertexCoverSet[j];
			seedT1[0] = curNode;
			spread.AddDat(curNode) = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
			Visited.AddDat(curNode) = true;
		}
		//NStools::PrintHaShTable(spread, spread.Len());
		for (size_t j = 0; j < vertexCoverSet.size(); j++) {
			int curNode = vertexCoverSet[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(curNode);
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int parentID = NI.GetInNId(e);
				if (!Visited.GetDat(parentID)) {
					qParent.push(parentID);
					Visited.AddDat(parentID) = true;
				}
			}
		}
		qParentTemp = qParent;
		VisitedTemp = Visited;
		for (size_t i = 0; i < 1; i++)
		{
			qParent = qParentTemp;
			Visited = VisitedTemp;
			////Calculate Influence Spread for  coverVertex set PARENTS
			while (!qParent.empty()) {
				int curNodeID = qParent.front();
				qParent.pop();
				PGraph::TObj::TNodeI CurNodeNI = Graph->GetNI(curNodeID);
				float x = 1;
				for (int e = 0; e < CurNodeNI.GetOutDeg(); e++) {
					int childID = CurNodeNI.GetOutNId(e);
					PGraph::TObj::TNodeI ChildNI = Graph->GetNI(childID);
					x += spread.GetDat(childID) * (1.0 / ChildNI.GetInDeg());
				}
				//spread.AddDat(curNodeID) = spread.GetDat(curNodeID) + x;
				spread.AddDat(curNodeID) = x;

				for (int e = 0; e < CurNodeNI.GetInDeg(); e++) {
					int parentID = CurNodeNI.GetInNId(e);
					if (!Visited.GetDat(parentID)) {
						qParent.push(parentID);
						Visited.AddDat(parentID) = true;
					}
				}
			}
		}


		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			TIntIntH item1;
			int curVer = NI.GetId();
			item1.AddDat(curVer) = spread.GetDat(curVer);
			qCELF.push(item1);
		}
		NScelf::QueueSortFirstK(qCELF, qCELF.size());
		//NScelf::printqueue(qCELF);
		//# Select the first node and remove from candidate list
		TIntIntH item2 = qCELF.front();
		int k = 0;
		seeds.push_back(item2.GetKey(0));
		infectedNodes = item2.GetDat(seeds[k]);
		qCELF.pop();
		cout << endl << "CELFPPPGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
		NStools::SaveToFile("CELFPPPGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELFPPPGlobal", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELFPPPGlobal", "RunningTime");
		//	# Find the next k - 1 nodes using the list - sorting procedure
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int infectedNodesLocal;
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = qCELF.front();
				int current = item3.GetKey(0);
				//Q[0] = (current, IC(g, S + [current], p, mc) - spread)
				//remove front of q (current) => then Evaluate the IC(s+current) => push cur and its influence to q
				qCELF.pop();
				seeds.push_back(current);
				infectedNodesLocal = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				TIntIntH item4;
				item4.AddDat(current) = infectedNodesLocal - infectedNodes;
				qCELF.push(item4);
				NScelf::QueueSortOthersK(qCELF, qCELF.size());
				TIntIntH item5 = qCELF.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = qCELF.front();
			seeds.push_back(item6.GetKey(0));
			infectedNodes += item6.GetDat(seeds[k]);
			qCELF.pop();
			cout << endl << "CELFPPPGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
			NStools::SaveToFile("CELFPPPGlobal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELFPPPGlobal", "InfluenceSpraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELFPPPGlobal", "RunningTime");
		}//end for
	}

	template<class PGraph>
	void CELFPPPLocal(const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC) {
		TExeTm ExeTm1;
		queue<TIntIntH> qCELF;
		queue<int> qParent;
		queue<int> qParentTemp;
		TIntBoolH Visited;
		TIntBoolH VisitedTemp;
		TIntFltH spread;
		int infectedNodes;

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			//spread.AddDat(NI.GetId()) = 0.0;
			spread.AddDat(NI.GetId()) = 1 + NI.GetOutDeg();
			Visited.AddDat(NI.GetId()) = false;
		}
		//NStools::PrintHaShTable(spread, spread.Len());
		vector<int> vertexCoverSet;
		NSsimPath::get_vertex_cover2(Graph, vertexCoverSet);

		//Calculate Influence Spread for coverVertex set
		vector <int> seedT1(1, -1111);
		for (size_t j = 0; j < vertexCoverSet.size(); j++) {
			int curNode = vertexCoverSet[j];
			seedT1[0] = curNode;
			spread.AddDat(curNode) = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
			Visited.AddDat(curNode) = true;
		}
		//NStools::PrintHaShTable(spread, spread.Len());
		for (size_t j = 0; j < vertexCoverSet.size(); j++) {
			int curNode = vertexCoverSet[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(curNode);
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int parentID = NI.GetInNId(e);
				if (!Visited.GetDat(parentID)) {
					qParent.push(parentID);
					Visited.AddDat(parentID) = true;
				}
			}
		}
		qParentTemp = qParent;
		VisitedTemp = Visited;
		for (size_t i = 0; i < 1; i++)
		{
			qParent = qParentTemp;
			Visited = VisitedTemp;
			////Calculate Influence Spread for  coverVertex set PARENTS
			while (!qParent.empty()) {
				int curNodeID = qParent.front();
				qParent.pop();
				PGraph::TObj::TNodeI CurNodeNI = Graph->GetNI(curNodeID);
				float x = 1;
				for (int e = 0; e < CurNodeNI.GetOutDeg(); e++) {
					int childID = CurNodeNI.GetOutNId(e);
					PGraph::TObj::TNodeI ChildNI = Graph->GetNI(childID);
					x += spread.GetDat(childID) * (1.0 / ChildNI.GetInDeg());
				}
				//spread.AddDat(curNodeID) = spread.GetDat(curNodeID) + x;
				spread.AddDat(curNodeID) = x;

				for (int e = 0; e < CurNodeNI.GetInDeg(); e++) {
					int parentID = CurNodeNI.GetInNId(e);
					if (!Visited.GetDat(parentID)) {
						qParent.push(parentID);
						Visited.AddDat(parentID) = true;
					}
				}
			}
		}


		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			TIntIntH item1;
			int curVer = NI.GetId();
			item1.AddDat(curVer) = spread.GetDat(curVer);
			qCELF.push(item1);
		}
		NScelf::QueueSortFirstK(qCELF, qCELF.size());
		//NScelf::printqueue(qCELF);
		//# Select the first node and remove from candidate list
		TIntIntH item2 = qCELF.front();
		int k = 0;
		seeds.push_back(item2.GetKey(0));
		infectedNodes = item2.GetDat(seeds[k]);
		qCELF.pop();
		NStools::SaveToFile("CELFPPPLocal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		//	# Find the next k - 1 nodes using the list - sorting procedure
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int infectedNodesLocal;
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = qCELF.front();
				int current = item3.GetKey(0);
				//Q[0] = (current, IC(g, S + [current], p, mc) - spread)
				//remove front of q (current) => then Evaluate the IC(s+current) => push cur and its influence to q
				qCELF.pop();
				seeds.push_back(current);
				infectedNodesLocal = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				TIntIntH item4;
				item4.AddDat(current) = infectedNodesLocal - infectedNodes;
				qCELF.push(item4);
				NScelf::QueueSortOthersK(qCELF, qCELF.size());
				TIntIntH item5 = qCELF.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = qCELF.front();
			seeds.push_back(item6.GetKey(0));
			infectedNodes += item6.GetDat(seeds[k]);
			qCELF.pop();
			NStools::SaveToFile("CELFPPPLocal", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		}//end for
	}

	template<class PGraph>
	void CELFPPPGlobalOfLocal(const vector<int>& candidateSet, const PGraph& Graph, const GlobalConst & GC, const char* Model, vector <int>& seeds, const int& SeedSize, const double &ICProb, const int& MC, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
		TExeTm ExeTm1;
		queue<TIntIntH> qCELF;
		queue<int> qParent;
		queue<int> qParentTemp;
		TIntBoolH Visited;
		TIntBoolH VisitedTemp;
		TIntFltH spread;
		int infectedNodes;

		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			//spread.AddDat(NI.GetId()) = 0.0;
			spread.AddDat(NI.GetId()) = 1 + NI.GetOutDeg();
			Visited.AddDat(NI.GetId()) = false;
		}
		//NStools::PrintHaShTable(spread, spread.Len());
		vector<int> vertexCoverSet;
		NSsimPath::get_vertex_cover2(Graph, vertexCoverSet);

		//Calculate Influence Spread for coverVertex set
		vector <int> seedT1(1, -1111);
		for (size_t j = 0; j < vertexCoverSet.size(); j++) {
			int curNode = vertexCoverSet[j];
			seedT1[0] = curNode;
			spread.AddDat(curNode) = NSIS::callInfluenceSpreadModel(Graph, seedT1, Model, ICProb, MC);
			Visited.AddDat(curNode) = true;
		}
		//NStools::PrintHaShTable(spread, spread.Len());
		for (size_t j = 0; j < vertexCoverSet.size(); j++) {
			int curNode = vertexCoverSet[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(curNode);
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int parentID = NI.GetInNId(e);
				if (!Visited.GetDat(parentID)) {
					qParent.push(parentID);
					Visited.AddDat(parentID) = true;
				}
			}
		}
		qParentTemp = qParent;
		VisitedTemp = Visited;
		for (size_t i = 0; i < 1; i++)
		{
			qParent = qParentTemp;
			Visited = VisitedTemp;
			////Calculate Influence Spread for  coverVertex set PARENTS
			while (!qParent.empty()) {
				int curNodeID = qParent.front();
				qParent.pop();
				PGraph::TObj::TNodeI CurNodeNI = Graph->GetNI(curNodeID);
				float x = 1;
				for (int e = 0; e < CurNodeNI.GetOutDeg(); e++) {
					int childID = CurNodeNI.GetOutNId(e);
					PGraph::TObj::TNodeI ChildNI = Graph->GetNI(childID);
					x += spread.GetDat(childID) * (1.0 / ChildNI.GetInDeg());
				}
				//spread.AddDat(curNodeID) = spread.GetDat(curNodeID) + x;
				spread.AddDat(curNodeID) = x;

				for (int e = 0; e < CurNodeNI.GetInDeg(); e++) {
					int parentID = CurNodeNI.GetInNId(e);
					if (!Visited.GetDat(parentID)) {
						qParent.push(parentID);
						Visited.AddDat(parentID) = true;
					}
				}
			}
		}


		for (size_t i = 0; i < candidateSet.size(); i++)
		{
			TIntIntH item1;
			item1.AddDat(candidateSet[i]) = spread.GetDat(candidateSet[i]);
			qCELF.push(item1);
		}

		NScelf::QueueSortFirstK(qCELF, qCELF.size());
		//NScelf::printqueue(qCELF);
		//# Select the first node and remove from candidate list
		TIntIntH item2 = qCELF.front();
		int k = 0;
		seeds.push_back(item2.GetKey(0));
		infectedNodes = item2.GetDat(seeds[k]);
		qCELF.pop();
		cout << endl << "CELFPPPGL   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
		NStools::SaveToFile("CELFPPPGL", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
		NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELFPPPGL", "InfluenceSpraed");
		NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELFPPPGL", "RunningTime");
				//	# Find the next k - 1 nodes using the list - sorting procedure
		for (k = 1; k < SeedSize; k++) {
			bool check = false;
			int infectedNodesLocal;
			while (!check) {
				//# Recalculate spread of top node
				TIntIntH item3 = qCELF.front();
				int current = item3.GetKey(0);
				//Q[0] = (current, IC(g, S + [current], p, mc) - spread)
				//remove front of q (current) => then Evaluate the IC(s+current) => push cur and its influence to q
				qCELF.pop();
				seeds.push_back(current);
				infectedNodesLocal = NSIS::callInfluenceSpreadModel(Graph, seeds, Model, ICProb, MC);
				seeds.pop_back();
				TIntIntH item4;
				item4.AddDat(current) = infectedNodesLocal - infectedNodes;
				qCELF.push(item4);
				NScelf::QueueSortOthersK(qCELF, qCELF.size());
				TIntIntH item5 = qCELF.front();
				check = (item5.GetKey(0) == current);
			}//end while
			 //# Select the next node
			TIntIntH item6 = qCELF.front();
			seeds.push_back(item6.GetKey(0));
			infectedNodes += item6.GetDat(seeds[k]);
			qCELF.pop();
			cout << endl << "CELFPPPGL   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << ExeTm1.GetSecs()<<" Sec";
			NStools::SaveToFile("CELFPPPGL", k + 1, infectedNodes, Graph->GetNodes(), MC, Model, ExeTm1.GetSecs(),GC);
			NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "CELFPPPGL", "InfluenceSpraed");
			NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, ExeTm1.GetSecs(), Model, "CELFPPPGL", "RunningTime");
			}//end for
	}

}
