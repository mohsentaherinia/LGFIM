#include "stdafx.h"



//vector<double> IRIE::od(MAX_NODE, 0);
//vector<int> IRIE::dd(MAX_NODE, 0);

TIntFltH newIRIEclass::AP;
TIntFltH newIRIEclass::dpV;
TIntBoolH newIRIEclass::Visited;
int newIRIEclass::k = 0;
bool newIRIEclass::saturated = false;




 void newIRIEclass::NewIRIEGlobal(const PNGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot)
{
	int recentSeedNode;
	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		AP.AddDat(NI.GetId(), 0.0);		//AP(u)=0
		dpV.AddDat(NI.GetId(), 1.0);		//AP(u)=0
		Visited.AddDat(NI.GetId(), false);		//AP(u)=0
	}
	k=0;
	saturated = false;
	for (; k < SeedSize; k++)
	{
		if (!k)
			recentSeedNode = initialInfluence(Graph, GC, seed, alpha, SeedSize, Model, ICProbb, MCS, curveInfoISV, curveInfoRTV, indexPlot);
		else
			recentSeedNode = residualInfluence(Graph, GC, seed, alpha, SeedSize, Model, ICProbb, MCS, curveInfoISV, curveInfoRTV, indexPlot);
		saturated = false;

		computeAP(Graph, seed, recentSeedNode, 320);

	}
}

TExeTm ExeTmR;

double ElapsedSimTimes = 0, EndTime;


int newIRIEclass::initialInfluence(const PNGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot)
{
	
	int NNodes = Graph->GetNodes();
	bool changed = true;
	int i = 0;
	double thetaa = 0.0001;
	double edgeProb;
	long count = 0;

	TIntFltH newdpV;

	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		newdpV.AddDat(NI.GetId(), 0.0);		
		dpV.AddDat(NI.GetId(), 1.0);		

	}

	while (count++ < NUM_SUBS_LOOP && changed)
	{
		changed = false;

		for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			newdpV.AddDat(NI.GetId(), 0.0);
		}

		float value;
		int curNodeiuID, neighbourjvID;
		for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			curNodeiuID = NI.GetId();
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				neighbourjvID = NI.GetOutNId(e);
				PNGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
				edgeProb = (float)1.0 / NI2.GetInDeg();  //edgeProb = exp(-e.w1);
				if (edgeProb > EPSs)
				{
					value = newdpV.GetDat(curNodeiuID);
					value += alpha * edgeProb * dpV.GetDat(neighbourjvID);
					newdpV.AddDat(curNodeiuID) = value;
				}
			}
			newdpV.AddDat(curNodeiuID) = newdpV.GetDat(curNodeiuID) + 1;
		}
		int curNodei;
		for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			curNodei = NI.GetId();
			if ((newdpV.GetDat(curNodei) < dpV.GetDat(curNodei) - thetaa - EPSs) ||
				(newdpV.GetDat(curNodei) > dpV.GetDat(curNodei) + thetaa + EPSs))
				changed = true;
			dpV.AddDat(curNodei) = newdpV.GetDat(curNodei);
		}
	}

	
	int u = newIRIEclass::IRIEGetMax(Graph, seed);//U = Argmax dv(u)
	seed.push_back(u);//S = S + u
	//Visited.AddDat(u) = true;//V = V - u
							 //////////
	EndTime = ExeTmR.GetSecs();
	double exactRunTime = EndTime - ElapsedSimTimes;
	int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
	cout << endl << "newIRIEGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
	NStools::SaveToFile("newIRIEGlobal", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime, GC);
	NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "newIRIEGlobal", "InfluenceSpraed");
	NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "newIRIEGlobal", "RunningTime");
	ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	return u;


}

int newIRIEclass::residualInfluence(const PNGraph& Graph, const GlobalConst & GC, vector<int>& seed, const double& alpha, const int& SeedSize, const char* Model, const double& ICProbb, const int& MCS, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot)
{
	int NNodes = Graph->GetNodes();
	bool changed = true;
	int i = 0;
	double thetaa = 0.0001;
	double edgeProb;
	long count = 0;

	TIntFltH newdpV, newnewdpV;


	
	//new implemented IRIE
	/*
	for (i = 0;i < n;i++) newdp[i] = 0;

	for (i = 0;i < n;i++) dd[i]=Graph::GetNeighbor(i);

	for (i = 0;i < n;i++) dp[i] = 1 - ap[i];

	while(count++ < NUM_LOOP && changed && !saturated)
	{
	changed = false;

	for (i = 0;i < n;i++)
	newdp[i]=0;

	//cout << count << endl;

	for (i = 0;i < n;i++)
	{
	for (int j = 0; j < dd[i]; j++)
	{
	Edge e = Graph::GetEdge(i,j);

	edgeProb = exp(-e.w1);
	if(edgeProb > EPS || ap[i] < 1 - EPS)
	newdp[i] += alpha * (1 - ap[i]) * edgeProb * dp[e.v];


	//newdp[i] += alpha * exp(-e.w1) * dp[e.v];
	}
	newdp[i] += 1 - ap[i];
	}

	//delta=0;

	for (i = 0;i < n; i++)
	if((newdp[i] < dp[i] - theta - EPS) || (newdp[i] > dp[i] + theta + EPS))
	changed = true;

	for (i = 0;i < n; i++)
	if(newdp[i] > n)
	saturated = true;

	for (i = 0;i < n; i++)
	dp[i] = newdp[i];
	}

	int max = GetMax();

	list[round] = max;
	d[round] = dp[max];
	seed[max] = true;

	return max;
	*/

	//IRIE-1 initial
	/*
	for (i = 0;i < n;i++) newdp[i] = (1 - ap[i]) * dp[i];

	for (i = 0;i < n;i++)
	{
	for (int j = 0; j < dd[i]; j++)
	{
	Edge e = Graph::GetEdge(i,j);

	edgeProb = exp(-e.w1);
	if(edgeProb > EPS || ap[i] < 1 - EPS)
	newnewdp[i] += alpha * (1 - ap[i]) * edgeProb * newdp[e.v];


	//newdp[i] += alpha * exp(-e.w1) * dp[e.v];
	}
	newnewdp[i] += 1 - ap[i];
	}

	for (i = 0;i < n; i++)
	if(newnewdp[i] > n)
	saturated = true;

	double max = -1000000.0;
	int mp = -1;
	for (int j=0; j<n; j++)
	{
	double tmp = newnewdp[j];
	if (tmp >max)
	{
	max = tmp;
	mp = j;
	}
	}

	list[round] = mp;
	d[round] = newnewdp[mp];
	seed[mp] = true;

	return mp;
	*/

	//IRIE-1 previous
	/*
	for (i = 0;i < n;i++) newdp[i] = (1 - ap[i]) * dp[i];

	for (i = 0;i < n;i++)
	{
	for (int j = 0; j < dd[i]; j++)
	{
	Edge e = Graph::GetEdge(i,j);

	edgeProb = exp(-e.w1);
	if(edgeProb > EPS || ap[i] < 1 - EPS)
	newnewdp[i] += alpha * (1 - ap[i]) * edgeProb * newdp[e.v];
	}
	newnewdp[i] += 1 - ap[i];
	}

	for (i = 0;i < n; i++)
	dp[i] = newnewdp[i];

	for (i = 0;i < n; i++)
	if(newnewdp[i] > n)
	saturated = true;

	double max = -1000000.0;
	int mp = -1;
	for (int j=0; j<n; j++)
	{
	double tmp = newnewdp[j];
	if (tmp >max)
	{
	max = tmp;
	mp = j;
	}
	}

	list[round] = mp;
	d[round] = newnewdp[mp];
	seed[mp] = true;

	return mp;
	*/

	//IRIE-5 initial
	/*
	for (i = 0;i < n;i++) newdp[i] = (1 - ap[i]) * dp[i];

	while(count++ < 5 && changed && !saturated)
	{
	changed = false;

	for (i = 0;i < n;i++)
	newnewdp[i] = 0;

	//cout << count << endl;

	for (i = 0;i < n;i++)
	{
	for (int j = 0; j < dd[i]; j++)
	{
	Edge e = Graph::GetEdge(i,j);

	edgeProb = exp(-e.w1);
	if(edgeProb > EPS || ap[i] < 1 - EPS)
	newnewdp[i] += alpha * (1 - ap[i]) * edgeProb * newdp[e.v];


	//newdp[i] += alpha * exp(-e.w1) * dp[e.v];
	}
	newnewdp[i] += 1 - ap[i];
	}

	for (i = 0;i < n; i++)
	if((newnewdp[i] < newdp[i] - theta) || (newnewdp[i] > newdp[i] + theta))
	changed = true;

	for (i = 0;i < n; i++)
	if(newnewdp[i] > n)
	saturated = true;

	for (i = 0;i < n; i++)
	newdp[i] = newnewdp[i];

	}

	double max = -1000000.0;
	int mp = -1;
	for (int j=0; j<n; j++)
	{
	double tmp = newnewdp[j];
	if (tmp >max)
	{
	max = tmp;
	mp = j;
	}
	}

	list[round] = mp;
	d[round] = newnewdp[mp];
	seed[mp] = true;

	return mp;
	*/

	//IRIE-5 previous

	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		int curNodet = NI.GetId();
		newdpV.AddDat(curNodet) = (1 - AP.GetDat(curNodet))*dpV.GetDat(curNodet);
	}


	while (count++ < NUM_SUBS_LOOP && changed)
	{
		changed = false;

		for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int curNodet = NI.GetId();
			newnewdpV.AddDat(NI.GetId(), 0.0);
		}

		float value;
		int curNodeiuID, neighbourjvID;
		for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			curNodeiuID = NI.GetId();
			for (int e = 0; e < NI.GetOutDeg(); e++) {
				neighbourjvID = NI.GetOutNId(e);
				PNGraph::TObj::TNodeI NI2 = Graph->GetNI(neighbourjvID);
				edgeProb = (float)1.0 / NI2.GetInDeg();  //edgeProb = exp(-e.w1);
				if (edgeProb > EPSs|| AP.GetDat(curNodeiuID)<1-EPSs)
				{
					value = newnewdpV.GetDat(curNodeiuID);
					value += alpha *  (1 - AP.GetDat(curNodeiuID)) * edgeProb * newdpV.GetDat(neighbourjvID);
					newnewdpV.AddDat(curNodeiuID) = value;
				}
			}
			newnewdpV.AddDat(curNodeiuID) = newnewdpV.GetDat(curNodeiuID) + 1 - AP.GetDat(curNodeiuID);
		}

		int curNodei;
		for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			curNodei = NI.GetId();
			if ((newnewdpV.GetDat(curNodei) < newdpV.GetDat(curNodei) - thetaa) ||
				(newnewdpV.GetDat(curNodei) > newdpV.GetDat(curNodei) + thetaa))
				changed = true;
			if (newnewdpV.GetDat(curNodei) > NNodes)
				saturated = true;
			newdpV.AddDat(curNodei) = newnewdpV.GetDat(curNodei);
		}
		
	}

	int curNodei;
	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		curNodei = NI.GetId();
		dpV.AddDat(curNodei) = newnewdpV.GetDat(curNodei);
	}

	int u = newIRIEclass::IRIEGetMax(Graph, seed);//U = Argmax dv(u)
	seed.push_back(u);//S = S + u
	EndTime = ExeTmR.GetSecs();
	double exactRunTime = EndTime - ElapsedSimTimes;
	int infectedNodes = NSIS::callInfluenceSpreadModel(Graph, seed, Model, ICProbb, MCS);
	cout << endl << "newIRIEGlobal   K=" << k + 1 << " " << infectedNodes << " is " << (double)infectedNodes / Graph->GetNodes() << "  %Percent RunTime=" << exactRunTime << " Sec";
	NStools::SaveToFile("newIRIEGlobal", k + 1, infectedNodes, Graph->GetNodes(), MCS, Model, exactRunTime, GC);
	NSplot::setPlottingValue(curveInfoISV, indexPlot, k + 1, infectedNodes, Model, "newIRIEGlobal", "InfluenceSpraed");
	NSplot::setPlottingValue(curveInfoRTV, indexPlot, k + 1, exactRunTime, Model, "newIRIEGlobal", "RunningTime");
	ElapsedSimTimes += ExeTmR.GetSecs() - EndTime;
	return u;

}

int newIRIEclass::IRIEGetMax(const PNGraph& Graph, vector<int>& seed)
{
	int MaxIndexKey = -9999999;
	float MaxValue = -99999999;
	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		int curNode = NI.GetId();
		if (dpV.GetDat(curNode) > MaxValue)
			if (!NStools::findNodeInVector(seed, curNode)) {
				MaxIndexKey = curNode;
				MaxValue = dpV.GetDat(curNode);
			}
	}
	return MaxIndexKey;
}

void newIRIEclass::computeAP(const PNGraph& Graph, const vector<int> &seeds, const int &src, const int& thetaa)
{
	
	//int i = 0, j, k;
	
	const int n = Graph->GetNodes();
	int* vertex = new int[n];   //??
	bool* Visited = new bool[n]; //??
	Heap h;
	HeapNode hNode;
	HeapNode curNode;
	int minNode;
	double minDist;
	double threshold = log((double)thetaa);

	initHeap(&h, n + 1);
	
	int z = 0;
	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		Visited[z] = false;
		vertex[z] = NI.GetId();
		z++;
	}
	for (size_t index2 = 0; index2 < seeds.size() - 1 && seeds[index2] != NULL; index2++)
		Visited[findinArray(vertex, seeds[index2], n)] = true;

	for (int i = 0; i < n; i++) {
		if (!Visited[i] && vertex[i] != src) { //if (!Visited[i]) {
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
				PNGraph::TObj::TNodeI NI3 = Graph->GetNI(neighbourID);
				double edgeProb = (float)1.0 / NI3.GetInDeg();

				if (h.elements[h.index[L]].value > minDist + edgeProb + EPss)
					decreaseKeyHeap(&h, h.index[L], minDist + edgeProb);
			}
		}
		AP.AddDat(vertex[minNode]) = AP.GetDat(vertex[minNode]) + exp(-minDist);
		if (AP.GetDat(vertex[minNode]) > 1 - EPss)
			AP.AddDat(vertex[minNode]) = 1;
		removeMinHeap(&h);
	}
	freeHeap(&h);
}

int newIRIEclass::findinArray(const int vertex[], int neighbourID, int n)
{
	for (size_t i = 0; i < n; i++)
	{
		if (vertex[i] == neighbourID)
			return i;
	}
	return -1;
}

vector<int> newIRIEclass::callNewIRIEGlobal(const PNGraph &Graph, const GlobalConst & GC, const char* Model, const int& SeedSize, const double& ICProbb, const int& MCS, const double& alpha, vector<NSplot::curveInfo>& curveInfoISV, vector<NSplot::curveInfo>& curveInfoRTV, const int indexPlot) {
	vector<int> seedIRIE;
	newIRIEclass::NewIRIEGlobal(Graph, GC, seedIRIE, alpha, SeedSize, Model, ICProbb, MCS, curveInfoISV, curveInfoRTV, indexPlot);
	NSplot::plotAllInfluenceSpread(curveInfoISV);
	NSplot::plotAllInfluenceSpread(curveInfoRTV);
	NStools::SaveSeedSetToFile("NewIRIEGlobal", GC, seedIRIE);
	return seedIRIE;
}