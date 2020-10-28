#pragma once
#include "stdafx.h"

namespace NSsimPath {
	vector<int> callSimPath1(const PNGraph& Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & R, const double & L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callSimPath1_v2(const vector<int> vertexSet, const PNGraph& Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & R, const double & L);
	vector<int> callSimPath2(const PNGraph& Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & R, const double & L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callSimPath2_v2(const vector<int> vertexSet, const PNGraph& Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & R, const double & L);
	vector<int> callSimPath3(const PNGraph& Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & R, const double & L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot);
	vector<int> callSimPath3_v2(const vector<int> vertexSet, const PNGraph& Graph, const GlobalConst & GC, const char * Model, const int & SeedSize, const double & ICProbb, const int & MCS, const double & R, const double & L);
	void DifferenceVector(vector<int>& S1, vector<int>& S2, vector<int>& S3);
	void SetIntersection(vector<int>& S1, vector<int>& S2, vector<int>& S3);
	void QueueSort(queue<TIntFltH>& q);
	void printqueue(queue<TIntFltH> q);
	bool SatisfyCondition(const int & x, const int & y, const stack<int>& Q, list<int>* D, const vector<int>& Wset);
	void PeekUPTopL(vector<int>& Uset, queue<TIntFltH> qCELF, const int & L);

	void RemoveXfromQCELF(queue<TIntFltH>& qCELF, const int & item);

	void qCELFUpdate(queue<TIntFltH>& qCELF, const int & item, const float & value);

	template<class PGraph>
	void SparceTrippleEdgeData(const PGraph& Graph, TIntPrIntH &NodNodDat, const int& data) {
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int curNodeiuID = NI.GetId();
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbourjvID = NI.GetOutNId(e);
				TIntPr x(curNodeiuID, neighbourjvID);
				//x.Val1 = curNodeiuID;
				//x.Val2 = neighbourjvID;
				NodNodDat.AddDat(x) = data;
			}
		}
		//for (size_t i = 0; i < NodNodDat.Len(); i++) {
		//	//TIntPr x = NodNodDat.GetKey(i);
		//	//cout << "\n" << x.Val1() << "\t" << x.Val2() << "\t" << NodNodDat[i].Val;
		//	cout << "\n" << NodNodDat.GetKey(i).Val1() << "\t" << NodNodDat.GetKey(i).Val2() << "\t" << NodNodDat[i].Val;
		//}
		//cout << "\n";
	}

	template<class PGraph>
	void get_vertex_cover1(const PGraph& Graph, vector<int>& VertexCover) {
		const int N = Graph->GetNodes();
		const int M = Graph->GetEdges();
		TIntIntH dv;
		TIntBoolH Visited;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			dv.AddDat(NI.GetId(), NI.GetOutDeg());
			Visited.AddDat(NI.GetId(), false);
		}
		dv.SortByDat(false);
		int checked = 0;
		int i = 0;
		while (checked < M) {
			int curNode = dv.GetKey(i).Val;
			VertexCover.push_back(curNode);
			i++;
			PGraph::TObj::TNodeI NI = Graph->GetNI(curNode);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbourID = NI.GetOutNId(e);
				//if (!Visited.GetDat(neighbourID))
				Visited.AddDat(neighbourID, true);
				checked++;
			}
		}
	}

	template<class PGraph>
	void get_vertex_cover2(const PGraph& Graph, vector<int>& VertexCover) {
		//Treat graph as UnDIRECTED :
		//Treat graph as UnDIRECTED :
		//Treat graph as UnDIRECTED :
		//Treat graph as UnDIRECTED :
		const int N = Graph->GetNodes();
		const int M = Graph->GetEdges();
		TIntIntH degreeV;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			degreeV.AddDat(NI.GetId(), NI.GetDeg());
		}
		degreeV.SortByDat(false);
		//MohsenTNT::PrintHaShTable(degreeV, N);
		int checkedEdge = 0;
		int i = 0;
		while (checkedEdge < M) {
			int curNode = degreeV.GetKey(i).Val;
			VertexCover.push_back(curNode);
			i++;
			checkedEdge += degreeV.GetDat(curNode);
		}
	}

	template<class PGraph>
	void get_vertex_cover3(const PGraph& Graph, vector<int>& VertexCover) {
		const int N = Graph->GetNodes();
		const int M = Graph->GetEdges();
		TIntPrIntH VisitedNND;
		NSsimPath::SparceTrippleEdgeData(Graph, VisitedNND, 0);
		////////////////////////////////
		TIntIntH dv;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			dv.AddDat(NI.GetId(), NI.GetDeg());
		dv.SortByDat(false);
		////////////////////////////////
		int checked = 0;
		int i = 0;
		while (checked < M) {
			int curNode = dv.GetKey(i).Val;
			VertexCover.push_back(curNode);
			i++;
			PGraph::TObj::TNodeI NI = Graph->GetNI(curNode);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int neighbourID = NI.GetOutNId(e);
				TIntPr x(curNode, neighbourID);
				if (!VisitedNND.GetDat(x)) {
					VisitedNND.AddDat(x) = 1;
					checked++;
				}
			}
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourID = NI.GetInNId(e);
				TIntPr x(neighbourID, curNode);
				if (!VisitedNND.GetDat(x)) {
					VisitedNND.AddDat(x) = 1;
					checked++;
				}
			}
		}
	}

	template<class PGraph>
	void Forward(stack<int>&Q, list<int> * D, float &spd, float &pp, const double& R, const vector<int>& Wset, const vector<int>& Uset, const PGraph& Graph, TIntPrIntH& spdW_, const int &myu) {
		int x = Q.top();
		int y;
		PGraph::TObj::TNodeI NIx = Graph->GetNI(x);
		int count = 0;
		while (true) {
			//# any suitable chid is ok
			for (; count < NIx.GetOutDeg(); count++) {
				y = NIx.GetOutNId(count);
				if (SatisfyCondition(x, y, Q, D, Wset) == true)
					break;
			}
			//# no such child:
			if (count == NIx.GetOutDeg())
				return;

			PGraph::TObj::TNodeI NIy = Graph->GetNI(y);
			float EgdeProb = 1.0 / NIy.GetInDeg();

			if (pp * EgdeProb < R)
				D[x].push_back(y);
			else {
				Q.push(y);
				pp = pp * EgdeProb;
				spd = spd + pp;
				D[x].push_back(y);
				x = Q.top();
				//for v in U :
					//if v not in Q :
						//spdW_u[v] = spdW_u[v] + pp
				for (size_t k = 0; k < Uset.size(); k++) {
					bool isExist = false;
					int v = Uset[k];
					stack<int> tempQ = Q;
					while (!tempQ.empty()) {
						if (v == tempQ.top()) {
							isExist = true;
							break;
						}
						tempQ.pop();
					}
					if (isExist == false) {
						TIntPr myu_v(myu, v);
						if (spdW_.IsKey(myu_v))
							spdW_.AddDat(myu_v) = spdW_.GetDat(myu_v) + pp;
					}
				}
				NIx = Graph->GetNI(x);
				count = 0;
			}
		}
	}

	template<class PGraph>
	float backtrack(int& u, const double& R, const vector<int>& Wset, const vector<int>& Uset, const PGraph& Graph, TIntPrIntH& spdW_, const int &myu) {
		//this function enumarates all simple path starting from u
		stack <int>Q;//current node on the path
		Q.push(u);
		float spd = 1;//
		float pp = 1;//weight of the current path
		list<int> * D = new list<int>[Graph->GetMxNId()];
		//D=NULL D[x]  maintains out-neighbors of x that have been seen so far
		int v;
		while (!Q.empty()) {
			//The subroutine FORWARD is called repeatedly which 
			//gives a new pruned path that is undiscovered until now
			NSsimPath::Forward(Q, D, spd, pp, R, Wset, Uset, Graph, spdW_,myu);
			//the last node is removed from the path 
			//the variable pp is set accordingly.
			//Since u no longer exists in the current path, D[u] is deleted
			u = Q.top();
			Q.pop();
			D[u].clear();
			if (!Q.empty()) {
				v = Q.top();
				PGraph::TObj::TNodeI NIu = Graph->GetNI(u);
				pp = pp / (1.0 / NIu.GetInDeg());
			}
		}
		delete []D;
		return spd;
	}

	template<class PGraph>
	float simpath_spread( vector<int>& seed, const double& R, const vector<int>& Uset, const PGraph& Graph, vector<int>& Vset, TIntPrIntH& spdW_) {
		int spread = 0;
		vector<int> Wset;
		DifferenceVector(Vset, seed, Wset);

		//if U is None or spdW_ is None :
		  //spdW_ = np.zeros(graph.node_num + 1)
		if (Uset.size() == 0 || spdW_.Len() == 0){
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 0);
		}
		for (size_t i = 0; i < seed.size(); i++) {
			int  u = seed[i];
			Wset.push_back(u);
			spread = spread + backtrack(u, R, Wset, Uset, Graph, spdW_,u); // estimate σW(u)
			Wset.pop_back();
		}
		return spread;
	}


	template<class PGraph>
	void GetSimPath1(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const int& SeedSize, const double& R, const double& L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot, const char* Model,const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		const int N = Graph->GetNodes();
		const int MxNId = Graph->GetMxNId();
		queue<TIntFltH> qCELF;

		//C = set(get_vertex_cover(graph))
		vector<int> Cset;
		get_vertex_cover3(Graph, Cset);

		//V = set(graph.nodes)
		vector<int> Vset(N, -1111);
		int i = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Vset[i++] = NI.GetId();

		//V_C = V.difference(C)
		vector<int> V_Cset;
		DifferenceVector(Vset, Cset, V_Cset);

		//spread = np.zeros(graph.node_num + 1)
		//# spread[x] is spd of S + x
		vector<float> spread(MxNId, 0.0);

		//spdV_ = np.ones((graph.node_num + 1, graph.node_num + 1))
		/*float **spdV_ = new float*[MxNId];
		for (size_t i = 0; i < MxNId; i++)
			spdV_[i] = new float[MxNId];
		for (size_t i = 0; i < MxNId; i++)
			for (size_t j = 0; j < MxNId; j++)
				spdV_[i][j] = 1;*/
		TIntPrIntH spdV_;
		NSsimPath::SparceTrippleEdgeData(Graph, spdV_, 1);


		//for u in C :
		//	U = V_C.intersection(set(graph.get_parents(u)))
		//	spread[u] = simpath_spread(set([u]), r, U, graph, spdV_)
		for (size_t j = 0; j < Cset.size(); j++) {
			int u = Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(u);
			vector<int> parents;
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourInID = NI.GetInNId(e);
				parents.push_back(neighbourInID);
			}
			vector<int> Uset;
			NSsimPath::SetIntersection(V_Cset, parents, Uset);
			vector <int> set_u;
			set_u.push_back(u);
			spread[u] = simpath_spread(set_u, R, Uset ,Graph,Vset,spdV_);//simpath_spread(set([u]), r, U, graph, spdV_)
			TIntFltH item;
			item.AddDat(u) = spread[u];
			qCELF.push(item);
		}

		//for v in V_C :
		//	v_children = graph.get_children(v)
		//	for child in v_children :
		//		spread[v] = spread[v] + spdV_[child][v] * graph.get_weight(v, child)
		//	spread[v] = spread[v] + 1
		for (size_t j = 0; j < V_Cset.size(); j++) {
			int v = V_Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(v);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int child = NI.GetOutNId(e);
				PGraph::TObj::TNodeI childNI = Graph->GetNI(child);
				//????????????TIntPr child_v(child,v);
				//????????????TIntPr child_v(child,v);

				//spread[v] = spread[v] + spdV_[child][v] * (1.0 / childNI.GetInDeg());
				TIntPr child_v(child, v);
				if (spdV_.IsKey(child_v))
					spread[v] = spread[v] + spdV_.GetDat(child_v) * (1.0 / childNI.GetInDeg());
				//spread[v] = spread[v] + spdV_[v][child] * (1.0 / childNI.GetInDeg());
				//spread[v] = spread[v] + spread[child] *(1.0 / childNI.GetInDeg());
			}
			spread[v] = spread[v] + 1;
			TIntFltH item;
			item.AddDat(v) = spread[v];
			qCELF.push(item);
		}

		NSsimPath::QueueSort(qCELF);

		vector<int> Wset;
		Wset = Vset;
		int spd = 0;

		//checked = np.zeros(graph.node_num + 1)
		vector<bool> checked(MxNId, false);
		vector<int> Uset;

		//while (seed.size() < SeedSize) {
		//	TIntFltH x = qCELF.front();
		//	seed.push_back(x.GetKey(0));
		//	qCELF.pop();
		//}
		//for (int i = 0; i < MxNId; i++)
		//	delete spdV_[i];//delete the A object allocations.
		//delete[] spdV_;//delete the array of pointers

		//return;


		while (seed.size() < SeedSize) {
			//U = celf.topn(l)
			Uset.clear();
			NSsimPath::PeekUPTopL(Uset, qCELF, L);

			//spdW_ = np.ones((graph.node_num + 1, graph.node_num + 1))
			/*float **spdW_ = new float*[MxNId];
			for (size_t i = 0; i < MxNId; i++)
				spdW_[i] = new float[MxNId];
			*/
			TIntPrIntH spdW_;
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 1);


			//spdV_x = np.zeros(graph.node_num + 1)
			vector<float> spdV_x(MxNId, 0.0);

			//simpath_spread(S, r, U, graph, spdW_ = spdW_)
			NSsimPath::simpath_spread(seed, R, Uset, Graph, Vset, spdW_);

			//for x in U :
			//for s in S :
			//spdV_x[x] = spdV_x[x] + spdW_[s][x]
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				for (size_t j = 0; j < seed.size(); j++) {
					int s = seed[j];
					TIntPr x_s(x, s);
					if (spdV_.IsKey(x_s))
						spdV_x[x] = spdV_x[x] + spdW_.GetDat(x_s);
				}
			}

			//for x in U :
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				if (checked[x] != false) {
					//S.add(x)
					seed.push_back(x);
					EndTime = ExeTmR.GetSecs();
					double exactRunTime = EndTime - ElapsedSimTimes;
					int infectedNodes=NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
					cout << endl << "SimPath1   K=" << seed.size() << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
					NStools::SaveToFile("SimPath1", seed.size(), infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
					NSplot::setPlottingValue(curveInfoISV, indexPlot, seed.size(), infectedNodes, Model, "SimPath1", "InfluenceSpraed");
					NSplot::setPlottingValue(curveInfoRTV, indexPlot, seed.size(), exactRunTime, Model, "SimPath1", "RunningTime");
					ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;

					//W = W.difference(set([x]))
					vector <int>xset;
					xset.push_back(x);
					NSsimPath::DifferenceVector(Wset, xset, Wset);

					//spd = spread[x]
					spd = spread[x];

					//checked = np.zeros(graph.node_num + 1)
					for (size_t i = 0; i < checked.size(); i++)
						checked[i] = false;

					//celf.remove(x)
					NSsimPath::RemoveXfromQCELF(qCELF, x);

					//break;
					break;
				}
				else {
					//spread[x] = backtrack(x, r, W, None, None, graph) + spdV_x[x]
					vector<int> NoneUset;
					//float *NoneSpdw_ = new float[MxNId];
					TIntPrIntH NoneSpdw_;
					spread[x] = NSsimPath::backtrack(x, R, Wset, NoneUset, Graph, NoneSpdw_,x) + spdV_x[x];
					
					//checked[x] = 1;
					checked[x] = true;

					//celf.update(x, spread[x] - spd)
					NSsimPath::qCELFUpdate(qCELF, x, spread[x] - spd);
					NSsimPath::QueueSort(qCELF);
				}
			}
			
		}


	}

	template<class PGraph>
	void GetSimPath1_v2(const vector<int> vertexSet, const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const int& SeedSize, const double& R, const double& L, const char* Model, const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		const int N = Graph->GetNodes();
		const int MxNId = Graph->GetMxNId();
		queue<TIntFltH> qCELF;

		//C = set(get_vertex_cover(graph))
		vector<int> Cset;
		get_vertex_cover3(Graph, Cset);

		//V = set(graph.nodes)
		vector<int> Vset(N, -1111);
		int i = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Vset[i++] = NI.GetId();

		//V_C = V.difference(C)
		vector<int> V_Cset;
		DifferenceVector(Vset, Cset, V_Cset);

		//spread = np.zeros(graph.node_num + 1)
		//# spread[x] is spd of S + x
		vector<float> spread(MxNId, 0.0);

		//spdV_ = np.ones((graph.node_num + 1, graph.node_num + 1))
		/*float **spdV_ = new float*[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		spdV_[i] = new float[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		for (size_t j = 0; j < MxNId; j++)
		spdV_[i][j] = 1;*/
		TIntPrIntH spdV_;
		NSsimPath::SparceTrippleEdgeData(Graph, spdV_, 1);


		//for u in C :
		//	U = V_C.intersection(set(graph.get_parents(u)))
		//	spread[u] = simpath_spread(set([u]), r, U, graph, spdV_)
		for (size_t j = 0; j < Cset.size(); j++) {
			int u = Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(u);
			vector<int> parents;
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourInID = NI.GetInNId(e);
				parents.push_back(neighbourInID);
			}
			vector<int> Uset;
			NSsimPath::SetIntersection(V_Cset, parents, Uset);
			vector <int> set_u;
			set_u.push_back(u);
			spread[u] = simpath_spread(set_u, R, Uset, Graph, Vset, spdV_);//simpath_spread(set([u]), r, U, graph, spdV_)
			TIntFltH item;
			item.AddDat(u) = spread[u];
			qCELF.push(item);
		}

		//for v in V_C :
		//	v_children = graph.get_children(v)
		//	for child in v_children :
		//		spread[v] = spread[v] + spdV_[child][v] * graph.get_weight(v, child)
		//	spread[v] = spread[v] + 1
		for (size_t j = 0; j < V_Cset.size(); j++) {
			int v = V_Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(v);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int child = NI.GetOutNId(e);
				PGraph::TObj::TNodeI childNI = Graph->GetNI(child);
				//????????????TIntPr child_v(child,v);
				//????????????TIntPr child_v(child,v);

				//spread[v] = spread[v] + spdV_[child][v] * (1.0 / childNI.GetInDeg());
				TIntPr child_v(child, v);
				if (spdV_.IsKey(child_v))
					spread[v] = spread[v] + spdV_.GetDat(child_v) * (1.0 / childNI.GetInDeg());
				//spread[v] = spread[v] + spdV_[v][child] * (1.0 / childNI.GetInDeg());
				//spread[v] = spread[v] + spread[child] *(1.0 / childNI.GetInDeg());
			}
			spread[v] = spread[v] + 1;
			TIntFltH item;
			item.AddDat(v) = spread[v];
			qCELF.push(item);
		}

		NSsimPath::QueueSort(qCELF);

		vector<int> Wset;
		Wset = Vset;
		int spd = 0;

		//checked = np.zeros(graph.node_num + 1)
		vector<bool> checked(MxNId, false);
		vector<int> Uset;

		//while (seed.size() < SeedSize) {
		//	TIntFltH x = qCELF.front();
		//	seed.push_back(x.GetKey(0));
		//	qCELF.pop();
		//}
		//for (int i = 0; i < MxNId; i++)
		//	delete spdV_[i];//delete the A object allocations.
		//delete[] spdV_;//delete the array of pointers

		//return;


		while (seed.size() < SeedSize) {
			//U = celf.topn(l)
			Uset.clear();
			NSsimPath::PeekUPTopL(Uset, qCELF, L);

			//spdW_ = np.ones((graph.node_num + 1, graph.node_num + 1))
			/*float **spdW_ = new float*[MxNId];
			for (size_t i = 0; i < MxNId; i++)
			spdW_[i] = new float[MxNId];
			*/
			TIntPrIntH spdW_;
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 1);


			//spdV_x = np.zeros(graph.node_num + 1)
			vector<float> spdV_x(MxNId, 0.0);

			//simpath_spread(S, r, U, graph, spdW_ = spdW_)
			NSsimPath::simpath_spread(seed, R, Uset, Graph, Vset, spdW_);

			//for x in U :
			//for s in S :
			//spdV_x[x] = spdV_x[x] + spdW_[s][x]
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				for (size_t j = 0; j < seed.size(); j++) {
					int s = seed[j];
					TIntPr x_s(x, s);
					if (spdV_.IsKey(x_s))
						spdV_x[x] = spdV_x[x] + spdW_.GetDat(x_s);
				}
			}

			//for x in U :
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				if (checked[x] != false) {
					//S.add(x)
					seed.push_back(x);
					EndTime = ExeTmR.GetSecs();
					double exactRunTime = EndTime - ElapsedSimTimes;
					int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
					NStools::SaveToFile("SimPath1_v2", seed.size(), infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
					ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;

					//W = W.difference(set([x]))
					vector <int>xset;
					xset.push_back(x);
					NSsimPath::DifferenceVector(Wset, xset, Wset);

					//spd = spread[x]
					spd = spread[x];

					//checked = np.zeros(graph.node_num + 1)
					for (size_t i = 0; i < checked.size(); i++)
						checked[i] = false;

					//celf.remove(x)
					NSsimPath::RemoveXfromQCELF(qCELF, x);

					//break;
					break;
				}
				else {
					//spread[x] = backtrack(x, r, W, None, None, graph) + spdV_x[x]
					vector<int> NoneUset;
					//float *NoneSpdw_ = new float[MxNId];
					TIntPrIntH NoneSpdw_;
					spread[x] = NSsimPath::backtrack(x, R, Wset, NoneUset, Graph, NoneSpdw_, x) + spdV_x[x];

					//checked[x] = 1;
					checked[x] = true;

					//celf.update(x, spread[x] - spd)
					NSsimPath::qCELFUpdate(qCELF, x, spread[x] - spd);
					NSsimPath::QueueSort(qCELF);
				}
			}
			
		}


	}

	template<class PGraph>
	void GetSimPath2(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const int& SeedSize, const double& R, const double& L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot, const char* Model, const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		const int N = Graph->GetNodes();
		const int MxNId = Graph->GetMxNId();
		queue<TIntFltH> qCELF;

		//C = set(get_vertex_cover(graph))
		vector<int> Cset;
		get_vertex_cover3(Graph, Cset);

		//V = set(graph.nodes)
		vector<int> Vset(N, -1111);
		int i = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Vset[i++] = NI.GetId();

		//V_C = V.difference(C)
		vector<int> V_Cset;
		DifferenceVector(Vset, Cset, V_Cset);

		//spread = np.zeros(graph.node_num + 1)
		//# spread[x] is spd of S + x
		vector<float> spread(MxNId, 0.0);

		//spdV_ = np.ones((graph.node_num + 1, graph.node_num + 1))
		/*float **spdV_ = new float*[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		spdV_[i] = new float[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		for (size_t j = 0; j < MxNId; j++)
		spdV_[i][j] = 1;*/
		TIntPrIntH spdV_;
		NSsimPath::SparceTrippleEdgeData(Graph, spdV_, 1);


		//for u in C :
		//	U = V_C.intersection(set(graph.get_parents(u)))
		//	spread[u] = simpath_spread(set([u]), r, U, graph, spdV_)
		for (size_t j = 0; j < Cset.size(); j++) {
			int u = Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(u);
			vector<int> parents;
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourInID = NI.GetInNId(e);
				parents.push_back(neighbourInID);
			}
			vector<int> Uset;
			NSsimPath::SetIntersection(V_Cset, parents, Uset);
			vector <int> set_u;
			set_u.push_back(u);
			spread[u] = simpath_spread(set_u, R, Uset, Graph, Vset, spdV_);//simpath_spread(set([u]), r, U, graph, spdV_)
			TIntFltH item;
			item.AddDat(u) = spread[u];
			qCELF.push(item);
		}

		//for v in V_C :
		//	v_children = graph.get_children(v)
		//	for child in v_children :
		//		spread[v] = spread[v] + spdV_[child][v] * graph.get_weight(v, child)
		//	spread[v] = spread[v] + 1
		for (size_t j = 0; j < V_Cset.size(); j++) {
			int v = V_Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(v);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int child = NI.GetOutNId(e);
				PGraph::TObj::TNodeI childNI = Graph->GetNI(child);
				//????????????TIntPr child_v(child,v);
				//????????????TIntPr child_v(child,v);

				//spread[v] = spread[v] + spdV_[child][v] * (1.0 / childNI.GetInDeg());
				TIntPr v_child(v, child);
				if (spdV_.IsKey(v_child))
					spread[v] = spread[v] + spdV_.GetDat(v_child) * (1.0 / childNI.GetInDeg());
					//spread[v] = spread[v] + spread[child] *(1.0 / childNI.GetInDeg());
			}
			spread[v] = spread[v] + 1;
			TIntFltH item;
			item.AddDat(v) = spread[v];
			qCELF.push(item);
		}

		NSsimPath::QueueSort(qCELF);

		vector<int> Wset;
		Wset = Vset;
		int spd = 0;

		//checked = np.zeros(graph.node_num + 1)
		vector<bool> checked(MxNId, false);
		vector<int> Uset;

		//while (seed.size() < SeedSize) {
		//	TIntFltH x = qCELF.front();
		//	seed.push_back(x.GetKey(0));
		//	qCELF.pop();
		//}
		//for (int i = 0; i < MxNId; i++)
		//	delete spdV_[i];//delete the A object allocations.
		//delete[] spdV_;//delete the array of pointers

		//return;


		while (seed.size() < SeedSize) {
			//U = celf.topn(l)
			Uset.clear();
			NSsimPath::PeekUPTopL(Uset, qCELF, L);

			//spdW_ = np.ones((graph.node_num + 1, graph.node_num + 1))
			/*float **spdW_ = new float*[MxNId];
			for (size_t i = 0; i < MxNId; i++)
			spdW_[i] = new float[MxNId];
			*/
			TIntPrIntH spdW_;
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 1);


			//spdV_x = np.zeros(graph.node_num + 1)
			vector<float> spdV_x(MxNId, 0.0);

			//simpath_spread(S, r, U, graph, spdW_ = spdW_)
			NSsimPath::simpath_spread(seed, R, Uset, Graph, Vset, spdW_);

			//for x in U :
			//for s in S :
			//spdV_x[x] = spdV_x[x] + spdW_[s][x]
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				for (size_t j = 0; j < seed.size(); j++) {
					int s = seed[j];
					TIntPr x_s(x, s);
					if (spdV_.IsKey(x_s))
						spdV_x[x] = spdV_x[x] + spdW_.GetDat(x_s);
				}
			}

			//for x in U :
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				if (checked[x] != false) {
					//S.add(x)
					seed.push_back(x);
					EndTime = ExeTmR.GetSecs();
					double exactRunTime = EndTime - ElapsedSimTimes;
					int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
					cout << endl << "SimPath2   K=" << seed.size() << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
					NStools::SaveToFile("SimPath2", seed.size(), infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
					NSplot::setPlottingValue(curveInfoISV, indexPlot, seed.size(), infectedNodes, Model, "SimPath2", "InfluenceSpraed");
					NSplot::setPlottingValue(curveInfoRTV, indexPlot, seed.size(), exactRunTime, Model, "SimPath2", "RunningTime");
					ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;

					//W = W.difference(set([x]))
					vector <int>xset;
					xset.push_back(x);
					NSsimPath::DifferenceVector(Wset, xset, Wset);

					//spd = spread[x]
					spd = spread[x];

					//checked = np.zeros(graph.node_num + 1)
					for (size_t i = 0; i < checked.size(); i++)
						checked[i] = false;

					//celf.remove(x)
					NSsimPath::RemoveXfromQCELF(qCELF, x);

					//break;
					break;
				}
				else {
					//spread[x] = backtrack(x, r, W, None, None, graph) + spdV_x[x]
					vector<int> NoneUset;
					//float *NoneSpdw_ = new float[MxNId];
					TIntPrIntH NoneSpdw_;
					spread[x] = NSsimPath::backtrack(x, R, Wset, NoneUset, Graph, NoneSpdw_, x) + spdV_x[x];

					//checked[x] = 1;
					checked[x] = true;

					//celf.update(x, spread[x] - spd)
					NSsimPath::qCELFUpdate(qCELF, x, spread[x] - spd);
					NSsimPath::QueueSort(qCELF);
				}
			}
			
		}


	}

	template<class PGraph>
	void GetSimPath2_v2(const vector<int> vertexSet, const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const int& SeedSize, const double& R, const double& L, const char* Model, const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		const int N = Graph->GetNodes();
		const int MxNId = Graph->GetMxNId();
		queue<TIntFltH> qCELF;

		//C = set(get_vertex_cover(graph))
		vector<int> Cset;
		get_vertex_cover3(Graph, Cset);

		//V = set(graph.nodes)
		vector<int> Vset(N, -1111);
		int i = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Vset[i++] = NI.GetId();

		//V_C = V.difference(C)
		vector<int> V_Cset;
		DifferenceVector(Vset, Cset, V_Cset);

		//spread = np.zeros(graph.node_num + 1)
		//# spread[x] is spd of S + x
		vector<float> spread(MxNId, 0.0);

		//spdV_ = np.ones((graph.node_num + 1, graph.node_num + 1))
		/*float **spdV_ = new float*[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		spdV_[i] = new float[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		for (size_t j = 0; j < MxNId; j++)
		spdV_[i][j] = 1;*/
		TIntPrIntH spdV_;
		NSsimPath::SparceTrippleEdgeData(Graph, spdV_, 1);


		//for u in C :
		//	U = V_C.intersection(set(graph.get_parents(u)))
		//	spread[u] = simpath_spread(set([u]), r, U, graph, spdV_)
		for (size_t j = 0; j < Cset.size(); j++) {
			int u = Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(u);
			vector<int> parents;
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourInID = NI.GetInNId(e);
				parents.push_back(neighbourInID);
			}
			vector<int> Uset;
			NSsimPath::SetIntersection(V_Cset, parents, Uset);
			vector <int> set_u;
			set_u.push_back(u);
			spread[u] = simpath_spread(set_u, R, Uset, Graph, Vset, spdV_);//simpath_spread(set([u]), r, U, graph, spdV_)
			TIntFltH item;
			item.AddDat(u) = spread[u];
			qCELF.push(item);
		}

		//for v in V_C :
		//	v_children = graph.get_children(v)
		//	for child in v_children :
		//		spread[v] = spread[v] + spdV_[child][v] * graph.get_weight(v, child)
		//	spread[v] = spread[v] + 1
		for (size_t j = 0; j < V_Cset.size(); j++) {
			int v = V_Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(v);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int child = NI.GetOutNId(e);
				PGraph::TObj::TNodeI childNI = Graph->GetNI(child);
				//????????????TIntPr child_v(child,v);
				//????????????TIntPr child_v(child,v);

				//spread[v] = spread[v] + spdV_[child][v] * (1.0 / childNI.GetInDeg());
				TIntPr v_child(v, child);
				if (spdV_.IsKey(v_child))
					spread[v] = spread[v] + spdV_.GetDat(v_child) * (1.0 / childNI.GetInDeg());
				//spread[v] = spread[v] + spread[child] *(1.0 / childNI.GetInDeg());
			}
			spread[v] = spread[v] + 1;
			TIntFltH item;
			item.AddDat(v) = spread[v];
			qCELF.push(item);
		}

		NSsimPath::QueueSort(qCELF);

		vector<int> Wset;
		Wset = Vset;
		int spd = 0;

		//checked = np.zeros(graph.node_num + 1)
		vector<bool> checked(MxNId, false);
		vector<int> Uset;

		//while (seed.size() < SeedSize) {
		//	TIntFltH x = qCELF.front();
		//	seed.push_back(x.GetKey(0));
		//	qCELF.pop();
		//}
		//for (int i = 0; i < MxNId; i++)
		//	delete spdV_[i];//delete the A object allocations.
		//delete[] spdV_;//delete the array of pointers

		//return;


		while (seed.size() < SeedSize) {
			//U = celf.topn(l)
			Uset.clear();
			NSsimPath::PeekUPTopL(Uset, qCELF, L);

			//spdW_ = np.ones((graph.node_num + 1, graph.node_num + 1))
			/*float **spdW_ = new float*[MxNId];
			for (size_t i = 0; i < MxNId; i++)
			spdW_[i] = new float[MxNId];
			*/
			TIntPrIntH spdW_;
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 1);


			//spdV_x = np.zeros(graph.node_num + 1)
			vector<float> spdV_x(MxNId, 0.0);

			//simpath_spread(S, r, U, graph, spdW_ = spdW_)
			NSsimPath::simpath_spread(seed, R, Uset, Graph, Vset, spdW_);

			//for x in U :
			//for s in S :
			//spdV_x[x] = spdV_x[x] + spdW_[s][x]
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				for (size_t j = 0; j < seed.size(); j++) {
					int s = seed[j];
					TIntPr x_s(x, s);
					if (spdV_.IsKey(x_s))
						spdV_x[x] = spdV_x[x] + spdW_.GetDat(x_s);
				}
			}

			//for x in U :
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				if (checked[x] != false) {
					//S.add(x)
					seed.push_back(x);
					EndTime = ExeTmR.GetSecs();
					double exactRunTime = EndTime - ElapsedSimTimes;
					int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
					NStools::SaveToFile("SimPath2_v2", seed.size(), infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
					ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
					//W = W.difference(set([x]))
					vector <int>xset;
					xset.push_back(x);
					NSsimPath::DifferenceVector(Wset, xset, Wset);

					//spd = spread[x]
					spd = spread[x];

					//checked = np.zeros(graph.node_num + 1)
					for (size_t i = 0; i < checked.size(); i++)
						checked[i] = false;

					//celf.remove(x)
					NSsimPath::RemoveXfromQCELF(qCELF, x);

					//break;
					break;
				}
				else {
					//spread[x] = backtrack(x, r, W, None, None, graph) + spdV_x[x]
					vector<int> NoneUset;
					//float *NoneSpdw_ = new float[MxNId];
					TIntPrIntH NoneSpdw_;
					spread[x] = NSsimPath::backtrack(x, R, Wset, NoneUset, Graph, NoneSpdw_, x) + spdV_x[x];

					//checked[x] = 1;
					checked[x] = true;

					//celf.update(x, spread[x] - spd)
					NSsimPath::qCELFUpdate(qCELF, x, spread[x] - spd);
					NSsimPath::QueueSort(qCELF);
				}
			}
			
		}


	}

	template<class PGraph>
	void GetSimPath3(const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const int& SeedSize, const double& R, const double& L, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot, const char* Model, const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		const int N = Graph->GetNodes();
		const int MxNId = Graph->GetMxNId();
		queue<TIntFltH> qCELF;

		//C = set(get_vertex_cover(graph))
		vector<int> Cset;
		get_vertex_cover3(Graph, Cset);

		//V = set(graph.nodes)
		vector<int> Vset(N, -1111);
		int i = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Vset[i++] = NI.GetId();

		//V_C = V.difference(C)
		vector<int> V_Cset;
		DifferenceVector(Vset, Cset, V_Cset);

		//spread = np.zeros(graph.node_num + 1)
		//# spread[x] is spd of S + x
		vector<float> spread(MxNId, 0.0);

		//spdV_ = np.ones((graph.node_num + 1, graph.node_num + 1))
		/*float **spdV_ = new float*[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		spdV_[i] = new float[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		for (size_t j = 0; j < MxNId; j++)
		spdV_[i][j] = 1;*/
		TIntPrIntH spdV_;
		NSsimPath::SparceTrippleEdgeData(Graph, spdV_, 1);


		//for u in C :
		//	U = V_C.intersection(set(graph.get_parents(u)))
		//	spread[u] = simpath_spread(set([u]), r, U, graph, spdV_)
		for (size_t j = 0; j < Cset.size(); j++) {
			int u = Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(u);
			vector<int> parents;
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourInID = NI.GetInNId(e);
				parents.push_back(neighbourInID);
			}
			vector<int> Uset;
			NSsimPath::SetIntersection(V_Cset, parents, Uset);
			vector <int> set_u;
			set_u.push_back(u);
			spread[u] = simpath_spread(set_u, R, Uset, Graph, Vset, spdV_);//simpath_spread(set([u]), r, U, graph, spdV_)
			TIntFltH item;
			item.AddDat(u) = spread[u];
			qCELF.push(item);
		}

		//for v in V_C :
		//	v_children = graph.get_children(v)
		//	for child in v_children :
		//		spread[v] = spread[v] + spdV_[child][v] * graph.get_weight(v, child)
		//	spread[v] = spread[v] + 1
		for (size_t j = 0; j < V_Cset.size(); j++) {
			int v = V_Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(v);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int child = NI.GetOutNId(e);
				PGraph::TObj::TNodeI childNI = Graph->GetNI(child);
				//????????????TIntPr child_v(child,v);
				//????????????TIntPr child_v(child,v);

				//spread[v] = spread[v] + spdV_[child][v] * (1.0 / childNI.GetInDeg());
				/*TIntPr child_v(child, v);
				if (spdV_.IsKey(child_v))
					spread[v] = spread[v] + spdV_.GetDat(child_v) * (1.0 / childNI.GetInDeg());
				*///spread[v] = spread[v] + spdV_[v][child] * (1.0 / childNI.GetInDeg());
				spread[v] = spread[v] + spread[child] *(1.0 / childNI.GetInDeg());
			}
			spread[v] = spread[v] + 1;
			TIntFltH item;
			item.AddDat(v) = spread[v];
			qCELF.push(item);
		}

		NSsimPath::QueueSort(qCELF);

		vector<int> Wset;
		Wset = Vset;
		int spd = 0;

		//checked = np.zeros(graph.node_num + 1)
		vector<bool> checked(MxNId, false);
		vector<int> Uset;

		//while (seed.size() < SeedSize) {
		//	TIntFltH x = qCELF.front();
		//	seed.push_back(x.GetKey(0));
		//	qCELF.pop();
		//}
		//for (int i = 0; i < MxNId; i++)
		//	delete spdV_[i];//delete the A object allocations.
		//delete[] spdV_;//delete the array of pointers

		//return;


		while (seed.size() < SeedSize) {
			//U = celf.topn(l)
			Uset.clear();
			NSsimPath::PeekUPTopL(Uset, qCELF, L);

			//spdW_ = np.ones((graph.node_num + 1, graph.node_num + 1))
			/*float **spdW_ = new float*[MxNId];
			for (size_t i = 0; i < MxNId; i++)
			spdW_[i] = new float[MxNId];
			*/
			TIntPrIntH spdW_;
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 1);


			//spdV_x = np.zeros(graph.node_num + 1)
			vector<float> spdV_x(MxNId, 0.0);

			//simpath_spread(S, r, U, graph, spdW_ = spdW_)
			NSsimPath::simpath_spread(seed, R, Uset, Graph, Vset, spdW_);

			//for x in U :
			//for s in S :
			//spdV_x[x] = spdV_x[x] + spdW_[s][x]
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				for (size_t j = 0; j < seed.size(); j++) {
					int s = seed[j];
					TIntPr x_s(x, s);
					if (spdV_.IsKey(x_s))
						spdV_x[x] = spdV_x[x] + spdW_.GetDat(x_s);
				}
			}

			//for x in U :
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				if (checked[x] != false) {
					//S.add(x)
					seed.push_back(x);
					EndTime = ExeTmR.GetSecs();
					double exactRunTime = EndTime - ElapsedSimTimes;
					int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
					cout << endl << "SimPath3   K=" << seed.size() << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime<<" Sec";
					NStools::SaveToFile("SimPath3", seed.size(), infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
					NSplot::setPlottingValue(curveInfoISV, indexPlot, seed.size(), infectedNodes, Model, "SimPath3", "InfluenceSpraed");
					NSplot::setPlottingValue(curveInfoRTV, indexPlot, seed.size(), exactRunTime, Model, "SimPath3", "RunningTime");
					ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;

					//W = W.difference(set([x]))
					vector <int>xset;
					xset.push_back(x);
					NSsimPath::DifferenceVector(Wset, xset, Wset);

					//spd = spread[x]
					spd = spread[x];

					//checked = np.zeros(graph.node_num + 1)
					for (size_t i = 0; i < checked.size(); i++)
						checked[i] = false;

					//celf.remove(x)
					NSsimPath::RemoveXfromQCELF(qCELF, x);

					//break;
					break;
				}
				else {
					//spread[x] = backtrack(x, r, W, None, None, graph) + spdV_x[x]
					vector<int> NoneUset;
					//float *NoneSpdw_ = new float[MxNId];
					TIntPrIntH NoneSpdw_;
					spread[x] = NSsimPath::backtrack(x, R, Wset, NoneUset, Graph, NoneSpdw_, x) + spdV_x[x];

					//checked[x] = 1;
					checked[x] = true;

					//celf.update(x, spread[x] - spd)
					NSsimPath::qCELFUpdate(qCELF, x, spread[x] - spd);
					NSsimPath::QueueSort(qCELF);
				}
			}
			
		}


	}

	template<class PGraph>
	void GetSimPath3_v2(const vector<int> vertexSet, const PGraph& Graph, const GlobalConst & GC, vector<int>& seed, const int& SeedSize, const double& R, const double& L, const char* Model, const double& ICProbb, const int& MCS)
	{
		TExeTm ExeTmR;
		double ElapsedSimTimes = 0, EndTime;
		const int N = Graph->GetNodes();
		const int MxNId = Graph->GetMxNId();
		queue<TIntFltH> qCELF;

		//C = set(get_vertex_cover(graph))
		vector<int> Cset;
		get_vertex_cover3(Graph, Cset);

		//V = set(graph.nodes)
		vector<int> Vset(N, -1111);
		int i = 0;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++)
			Vset[i++] = NI.GetId();

		//V_C = V.difference(C)
		vector<int> V_Cset;
		DifferenceVector(Vset, Cset, V_Cset);

		//spread = np.zeros(graph.node_num + 1)
		//# spread[x] is spd of S + x
		vector<float> spread(MxNId, 0.0);

		//spdV_ = np.ones((graph.node_num + 1, graph.node_num + 1))
		/*float **spdV_ = new float*[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		spdV_[i] = new float[MxNId];
		for (size_t i = 0; i < MxNId; i++)
		for (size_t j = 0; j < MxNId; j++)
		spdV_[i][j] = 1;*/
		TIntPrIntH spdV_;
		NSsimPath::SparceTrippleEdgeData(Graph, spdV_, 1);


		//for u in C :
		//	U = V_C.intersection(set(graph.get_parents(u)))
		//	spread[u] = simpath_spread(set([u]), r, U, graph, spdV_)
		for (size_t j = 0; j < Cset.size(); j++) {
			int u = Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(u);
			vector<int> parents;
			for (int e = 0; e < NI.GetInDeg(); e++) {
				int neighbourInID = NI.GetInNId(e);
				parents.push_back(neighbourInID);
			}
			vector<int> Uset;
			NSsimPath::SetIntersection(V_Cset, parents, Uset);
			vector <int> set_u;
			set_u.push_back(u);
			spread[u] = simpath_spread(set_u, R, Uset, Graph, Vset, spdV_);//simpath_spread(set([u]), r, U, graph, spdV_)
			TIntFltH item;
			item.AddDat(u) = spread[u];
			qCELF.push(item);
		}

		//for v in V_C :
		//	v_children = graph.get_children(v)
		//	for child in v_children :
		//		spread[v] = spread[v] + spdV_[child][v] * graph.get_weight(v, child)
		//	spread[v] = spread[v] + 1
		for (size_t j = 0; j < V_Cset.size(); j++) {
			int v = V_Cset[j];
			PGraph::TObj::TNodeI NI = Graph->GetNI(v);
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				int child = NI.GetOutNId(e);
				PGraph::TObj::TNodeI childNI = Graph->GetNI(child);
				//????????????TIntPr child_v(child,v);
				//????????????TIntPr child_v(child,v);

				//spread[v] = spread[v] + spdV_[child][v] * (1.0 / childNI.GetInDeg());
				/*TIntPr child_v(child, v);
				if (spdV_.IsKey(child_v))
				spread[v] = spread[v] + spdV_.GetDat(child_v) * (1.0 / childNI.GetInDeg());
				*///spread[v] = spread[v] + spdV_[v][child] * (1.0 / childNI.GetInDeg());
				spread[v] = spread[v] + spread[child] * (1.0 / childNI.GetInDeg());
			}
			spread[v] = spread[v] + 1;
			TIntFltH item;
			item.AddDat(v) = spread[v];
			qCELF.push(item);
		}

		NSsimPath::QueueSort(qCELF);

		vector<int> Wset;
		Wset = Vset;
		int spd = 0;

		//checked = np.zeros(graph.node_num + 1)
		vector<bool> checked(MxNId, false);
		vector<int> Uset;

		//while (seed.size() < SeedSize) {
		//	TIntFltH x = qCELF.front();
		//	seed.push_back(x.GetKey(0));
		//	qCELF.pop();
		//}
		//for (int i = 0; i < MxNId; i++)
		//	delete spdV_[i];//delete the A object allocations.
		//delete[] spdV_;//delete the array of pointers

		//return;


		while (seed.size() < SeedSize) {
			//U = celf.topn(l)
			Uset.clear();
			NSsimPath::PeekUPTopL(Uset, qCELF, L);

			//spdW_ = np.ones((graph.node_num + 1, graph.node_num + 1))
			/*float **spdW_ = new float*[MxNId];
			for (size_t i = 0; i < MxNId; i++)
			spdW_[i] = new float[MxNId];
			*/
			TIntPrIntH spdW_;
			NSsimPath::SparceTrippleEdgeData(Graph, spdW_, 1);


			//spdV_x = np.zeros(graph.node_num + 1)
			vector<float> spdV_x(MxNId, 0.0);

			//simpath_spread(S, r, U, graph, spdW_ = spdW_)
			NSsimPath::simpath_spread(seed, R, Uset, Graph, Vset, spdW_);

			//for x in U :
			//for s in S :
			//spdV_x[x] = spdV_x[x] + spdW_[s][x]
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				for (size_t j = 0; j < seed.size(); j++) {
					int s = seed[j];
					TIntPr x_s(x, s);
					if (spdV_.IsKey(x_s))
						spdV_x[x] = spdV_x[x] + spdW_.GetDat(x_s);
				}
			}

			//for x in U :
			for (size_t i = 0; i < Uset.size(); i++) {
				int x = Uset[i];
				if (checked[x] != false) {
					//S.add(x)
					seed.push_back(x);
					EndTime = ExeTmR.GetSecs();
					double exactRunTime = EndTime - ElapsedSimTimes;
					int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
					NStools::SaveToFile("SimPath3_v3", seed.size(), infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime,GC);
					ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;

					//W = W.difference(set([x]))
					vector <int>xset;
					xset.push_back(x);
					NSsimPath::DifferenceVector(Wset, xset, Wset);

					//spd = spread[x]
					spd = spread[x];

					//checked = np.zeros(graph.node_num + 1)
					for (size_t i = 0; i < checked.size(); i++)
						checked[i] = false;

					//celf.remove(x)
					NSsimPath::RemoveXfromQCELF(qCELF, x);

					//break;
					break;
				}
				else {
					//spread[x] = backtrack(x, r, W, None, None, graph) + spdV_x[x]
					vector<int> NoneUset;
					//float *NoneSpdw_ = new float[MxNId];
					TIntPrIntH NoneSpdw_;
					spread[x] = NSsimPath::backtrack(x, R, Wset, NoneUset, Graph, NoneSpdw_, x) + spdV_x[x];

					//checked[x] = 1;
					checked[x] = true;

					//celf.update(x, spread[x] - spd)
					NSsimPath::qCELFUpdate(qCELF, x, spread[x] - spd);
					NSsimPath::QueueSort(qCELF);
				}
			}
			
		}


	}

}
