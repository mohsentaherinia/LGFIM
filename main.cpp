#include "stdafx.h"
PNGraph init(int argc, char* argv[], const GlobalConst &GC);
void CompareSavePlotAll(const PNGraph &Graph, const GlobalConst &GC, vector<NSplot::curveInfo> &curveInfoISV, vector<NSplot::curveInfo> &curveInfoRTV, int &indexPlot);
void simulation1(int argc, char* argv[]);
void simulation2Scalibility();
void simulation3Scalibility();
void startSim(const PNGraph &Graph, const GlobalConst &GC);
struct dsProperties {
	char* name;
	int coef1;
	float coef2;
};

int main(int argc, char* argv[]) {
	
	//simulation2Scalibility();
	//simulation3Scalibility();
	simulation1(argc, argv);
	//simulation4CompareWithCELF(argc, argv);
	system("pause");
	return 0;
}



void simulation2Scalibility() {
	char * CDalgName[4] = { "Girvan-Newman","Clauset-Newman-Moore","Infomap","CLPA" };
	char *diffusionModel[3] = { "LTPfull","ICPfull","WCPfull" };
	const int datasetCount = 9;
	dsProperties datasett[datasetCount] = {
		{"../0syn V=1000 E=10000.txt", 1 , 0.001},
		{ "../1syn V=2000 E=20000.txt", 2 , 0.001 },
		{ "../2syn V=4000 E=40000.txt", 4 , 0.001 },
		{ "../3syn V=8000 E=80000.txt", 8 , 0.001 },
		{ "../4syn V=16000 E=160000.txt", 16 , 0.001 },
		{ "../5syn V=32000 E=320000.txt", 32 , 0.001 },
		{ "../6syn V=64000 E=640000.txt", 64 , 0.001 },
		{ "../7syn V=128000 E=1280000.txt", 128 , 0.001 },
		{ "../8syn V=256000 E=2560000.txt" ,256 , 0.001 }
	};
	
	GlobalConst GC;
	GC._MCConst = 10000;				 //Must be 24 at least for OMP;
	GC._MCGreedyconst = 10000;		 //Must be 24 at least for OMP;
	GC._ICProbConsttemp = 0.01;
	GC._seedSizeConst = 50; 
	GC._NodeDivision = 1;
	GC._ComDetMethod = CDalgName[3];
	GC._SixDegree = 6;
	GC._simMode = 1;
	FILE *F = fopen("__Info.not", "at+");
	fprintf(F, "\n\n\nInfluence Maximization. build: %s, %s. Time: %s\nMCconst = %d\nMCGreedyconst = %d\nNodeDivision = %d\nseedSize = %d\nICProb = %Lf\n", __TIME__, __DATE__, TExeTm::GetCurTm(), GC._MCConst, GC._MCGreedyconst, GC._NodeDivision, GC._seedSizeConst, GC._ICProbConsttemp);
	fclose(F);


	for (size_t DS = 0; DS <datasetCount; DS++)
	{
		
		GC._inputFileName = datasett[DS].name;
		GC._coef1 = datasett[DS].coef1;
		GC._coef2 = datasett[DS].coef2;
		
		/*
		int V = pow(2, DS) * 1000;
		int E = V * 10;
		PNGraph Graph = TSnap::GenRndGnm<PNGraph>(V,E);
		{ TFOut FOut(GC._inputFileName); Graph->Save(FOut); }*/

		TFIn FIn(GC._inputFileName);
		PNGraph Graph = TNGraph::Load(FIn);
		

		FILE *F = fopen("__Info.not", "at+");
		fprintf(F, "\n\n\nDatasetName = %s\t#Nodes = %d\t#Edges = %d\tcoef1 = %d\tcoef2 = %Lf", GC._inputFileName, Graph->GetNodes(), Graph->GetEdges(), GC._coef1, GC._coef2);
		fclose(F);
		printf("Treat graph as \"DIRECTED Graph\": \n\n");
		NSMemoryConsumation::memoryUsage("Initialize", "", GC);
		for (size_t df = 2; df < 3; df++)
		{
			GC._DifModel = diffusionModel[df];
			startSim(Graph, GC);
		}
	}
}

void simulation3Scalibility() {
	char * CDalgName[4] = { "Girvan-Newman","Clauset-Newman-Moore","Infomap","CLPA" };
	char *diffusionModel[3] = { "LTPfull","ICPfull","WCPfull" };
	const int datasetCount = 9;
	dsProperties datasett[datasetCount] = {
		{"../0syn V=2000 E=2000.txt",4,0.001},
		{ "../1syn V=2000 E=4000.txt",4,0.001 },
		{ "../2syn V=2000 E=8000.txt",4,0.001 },
		{ "../3syn V=2000 E=16000.txt",4,0.001 },
		{ "../4syn V=2000 E=32000.txt",4,0.001 },
		{ "../5syn V=2000 E=64000.txt",4,0.001 },
		{ "../6syn V=2000 E=128000.txt",4,0.001 },
		{ "../7syn V=2000 E=256000.txt",4,0.001 },
		{ "../8syn V=2000 E=512000.txt",4,0.001 }
	};
	
	GlobalConst GC;
	GC._MCConst = 10000;				 //Must be 24 at least for OMP;
	GC._MCGreedyconst = 10000;		 //Must be 24 at least for OMP;
	GC._ICProbConsttemp = 0.01;
	GC._seedSizeConst = 50;
	GC._NodeDivision = 1;
	GC._ComDetMethod = CDalgName[3];
	GC._SixDegree = 6;
	GC._simMode = 2;
	FILE *F = fopen("__Info.not", "at+");
	fprintf(F, "\n\n\nInfluence Maximization. build: %s, %s. Time: %s\nMCconst = %d\nMCGreedyconst = %d\nNodeDivision = %d\nseedSize = %d\nICProb = %Lf\n", __TIME__, __DATE__, TExeTm::GetCurTm(), GC._MCConst, GC._MCGreedyconst, GC._NodeDivision, GC._seedSizeConst, GC._ICProbConsttemp);
	fclose(F);

	
	for (size_t DS = 0; DS <datasetCount; DS++)
	{

		GC._inputFileName = datasett[DS].name;
		GC._coef1 = datasett[DS].coef1;
		GC._coef2 = datasett[DS].coef2;
		
		/*
		int V = 2000;
		int E = pow(2, DS+1) * 1000;
		PNGraph Graph = TSnap::GenRndGnm<PNGraph>(V, E);
		{ TFOut FOut(GC._inputFileName); Graph->Save(FOut); }*/

		
		TFIn FIn(GC._inputFileName);
		PNGraph Graph = TNGraph::Load(FIn);

		FILE *F = fopen("__Info.not", "at+");
		fprintf(F, "\n\n\nDatasetName = %s\t#Nodes = %d\t#Edges = %d\tcoef1 = %d\tcoef2 = %Lf", GC._inputFileName, Graph->GetNodes(), Graph->GetEdges(), GC._coef1, GC._coef2);
		fclose(F);
		printf("Treat graph as \"DIRECTED Graph\": \n\n");
		NSMemoryConsumation::memoryUsage("Initialize", "", GC);
		for (size_t df = 2; df < 3; df++)
		{
			GC._DifModel = diffusionModel[df];
			startSim(Graph, GC);
		}
	}
}

void simulation1(int argc, char* argv[]) {
	char * CDalgName[4] = { "Girvan-Newman","Clauset-Newman-Moore","Infomap","CLPA" };
	char *diffusionModel[3] = { "LTPfull","ICPfull","WCPfull" };

	const int datasetCount = 16;
	dsProperties datasett[datasetCount] = {
		{ "../0NetHEPT.txt", 40 , 0.0005 },
		{ "../99out.advogato.txt", 10 , 0.0015 },
		{ "../11p2p-Gnutella09.txt", 40 , 0.0004 },
		{ "../2Cit-HepTh.txt", 20 , 0.0003 },
		{ "../99out.ego-gplus.txt", 20 , 0.0005 },
		{ "../3soc-Epinions1.txt", 100 , 0.00001 },
		{ "../4Slashdot0811.txt", 30 , 0.0005 },
		{ "../5web-Stanford.txt", 30 , 0.0001 },
		{ "../10out.amazon0601.txt", 100 , 0.00005 },
		{ "../99out.digg-friends.txt", 1200 , 0.00003 },
		{ "../6web-Google.txt", 200 , 0.00002 },
		{ "../99out.youtube-u-growth.txt", 800 , 0.000015 },
		{ "../7soc-pokec-relationships.txt", 500 , 0.000001 },
		{ "../9wiki-topcats.txt", 500 , 0.0001 },
		{ "../99out.flickr-growth.txt", 800 , 0.000015 },
		{ "../8soc-LiveJournal1.txt" , 1000 , 0.000001 }
	};

	GlobalConst GC;
	GC._MCConst = 10000;				 //Must be 24 at least for OMP;
	GC._MCGreedyconst = 10000;		 //Must be 24 at least for OMP;
	GC._ICProbConsttemp = 0.01;
	GC._seedSizeConst = 50;
	GC._NodeDivision = 1;
	GC._ComDetMethod = CDalgName[3];
	GC._SixDegree = 6;
	GC._simMode = 0;
	FILE *F = fopen("__Info.not", "at+");
	fprintf(F, "\n\n\nInfluence Maximization. build: %s, %s. Time: %s\nMCconst = %d\nMCGreedyconst = %d\nNodeDivision = %d\nseedSize = %d\nICProb = %Lf\n", __TIME__, __DATE__, TExeTm::GetCurTm(), GC._MCConst, GC._MCGreedyconst, GC._NodeDivision, GC._seedSizeConst, GC._ICProbConsttemp);
	fclose(F);


	for (size_t DS = 0; DS <datasetCount; DS++)
	{
		GC._inputFileName = datasett[DS].name;
		GC._coef1 = datasett[DS].coef1;
		GC._coef2 = datasett[DS].coef2;
		PNGraph Graph = init(argc, argv, GC);
		NSMemoryConsumation::memoryUsage("Initialize", "", GC);
		for (size_t df = 2; df < 3; df++)
		{
			GC._DifModel = diffusionModel[df];
			startSim(Graph, GC);
		}
	}
}


void startSim(const PNGraph &Graph, const GlobalConst &GC) {
	//double Q = NScommunity::CommunityDetectionRun(Graph, GC._inputFileName + 3, GC._ComDetMethod);

	////Uncomment if is exist a huge community in main comminity detection phase
	////comment above main comminity detection phase
	//PNGraph LargeSubG;
	//NScommunity::constructSubGraphOfLargeCommunityFromFile(Graph, LargeSubG);
	//double Q = NScommunity::CommunityDetectionRun(LargeSubG, GC._inputFileName + 3, GC._ComDetMethod);
	
	//char *fn = GC._inputFileName;
	TExeTm ExeTm;
	FILE *F = fopen("__Info.not", "at+");
	fprintf(F, "\n----------------------------For Diffusion Model \"%s\"----------------", GC._DifModel);
	fprintf(F, "\nStarting  Time = %s", TSecTm::GetCurTm().GetTmStr().CStr());
	fflush(F);
	int indexPlot = 0;
	const int n = 8;
	vector<NSplot::curveInfo> curveInfoISV(n), curveInfoRTV(n);
	for (size_t i = 0; i < n; i++)
	{
		strcpy(curveInfoISV[i].dataSet, GC._inputFileName + 4);
		strcpy(curveInfoRTV[i].dataSet, GC._inputFileName + 4);
	}
	CDIM(Graph, GC._DifModel, GC._ComDetMethod, GC,GC._coef1,GC._coef2, curveInfoISV, curveInfoRTV, indexPlot);
	CompareSavePlotAll(Graph, GC, curveInfoISV, curveInfoRTV, indexPlot);
	fprintf(F, "\nEnding Time = %s\nRunning Time = %s", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTm.GetTmStr());
	fclose(F);

}

PNGraph init(int argc, char* argv[], const GlobalConst &GC) {
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("Influence Maximization. build: %s, %s. Time: %s\t", __TIME__, __DATE__, TExeTm::GetCurTm()));
	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", GC._inputFileName, "Input un/directed graph");
	printf("Loading %s...", InFNm.CStr());
	PNGraph Graph = TSnap::LoadEdgeList<PNGraph>(InFNm);
	//PNGraph Graph = TSnap::GenRndGnm<PNGraph>(500, 5000); // generate a random graph

	TIntV NIdV;
	for (PNGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		if (NIdV.Len() < Graph->GetNodes() / GC._NodeDivision) {
			NIdV.Add(NI.GetId());
			/*PNGraph::TObj::TEdgeI Eij = Graph->GetEI(1, 2);
			Eij.GetSrcNId();
			Eij.GetDstNId();
			*/
		}
	}
	
	////uncomment if the orginal greph is undirected -- convert to directed
	//for (TNGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
	//	////printf("edge (%d, %d)\n", EI.GetSrcNId(), EI.GetDstNId());
	//	Graph->AddEdge(EI.GetDstNId(), EI.GetSrcNId());
	//}

	PNGraph SubG = TSnap::GetSubGraph(Graph, NIdV);
	printf("\nnodes:%d  edges:%d\n", SubG->GetNodes(), SubG->GetEdges());
	////
	//printf("Graph (%d, %d)\n", SubG->GetNodes(), SubG->GetEdges());
	//for (PNGraph::TObj::TNodeI NI = SubG->BegNI(); NI < SubG->EndNI(); NI++) {
	//	printf("  %d: ", NI.GetId());
	//	for (int e = 0; e < NI.GetOutDeg(); e++) {
	//		printf(" %d", NI.GetOutNId(e));
	//	}
	//	printf("\n");
	//}
	////
	///////Save Running Information//////////////////////
	FILE *F = fopen("__Info.not", "at+");
	fprintf(F, "\n\n\nDatasetName = %s\t#Nodes = %d\t#Edges = %d\tcoef1 = %d\tcoef2 = %Lf", GC._inputFileName,SubG->GetNodes(), SubG->GetEdges(),GC._coef1,GC._coef2);
	fclose(F);
	printf("Treat graph as \"DIRECTED Graph\": \n\n");
	////////////////////////////////////////////////////
	return SubG;
}

void CompareSavePlotAll(const PNGraph &Graph, const GlobalConst &GC, vector<NSplot::curveInfo> &curveInfoISV, vector<NSplot::curveInfo> &curveInfoRTV, int &indexPlot) {
	FILE *F = fopen("__Info.not", "at+");
	vector<int> targetSeedSet;
	NSMemoryConsumation::memStruct mStart, mFinish;

	/*TExeTm ExeTime0;
	printf("\nMyLAIMGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "MyLAIMGlobal", GC);
	targetSeedSet = NSMyLAIM::callMyLAIMGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++, GC._SixDegree, 0.001);
	mFinish = NSMemoryConsumation::memoryUsage("After", "MyLAIMGlobal", GC);
	printf("\nMyLAIMGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime0.GetTmStr());
	fprintf(F, "\nMyLAIMGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime0.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	
	TExeTm ExeTime1;
	printf("\nDegDisGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "DegDisGlobal", GC);
	targetSeedSet = NSDegreeDiscount::callDegDisGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, .01, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "DegDisGlobal", GC);
	printf("\nDegDisGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime1.GetTmStr());
	fprintf(F, "\nDegDisGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime1.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	
	TExeTm ExeTime2;
	printf("\nMaxDegreeRankGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "MaxDegreeRankGlobal", GC);
	targetSeedSet = NSbaseLineRank::callMaxDegreeRankGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "MaxDegreeRankGlobal", GC);
	printf("\nMaxDegreeRankGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime2.GetTmStr());
	fprintf(F, "\nMaxDegreeRankGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime2.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
	
	TExeTm ExeTime3;
	printf("\nPageRank2Global is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "PageRank2Global", GC);
	targetSeedSet = NSbaseLineRank::callPageRank2Global(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "PageRank2Global", GC);
	printf("\nPageRank2Global was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime3.GetTmStr());
	fprintf(F, "\nPageRank2Global Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime3.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);


	TExeTm ExeTime7;
	printf("\nIRIEPMIAGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "IRIEPMIAGlobal", GC);
	targetSeedSet = NSIRIE::callIRIEPMIAGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, 0.7, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "IRIEPMIAGlobal", GC);
	printf("\nIRIEPMIAGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime7.GetTmStr());
	fprintf(F, "\nIRIEPMIAGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime7.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
*/

		/*TExeTm ExeTime5;
		printf("\nIRGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
		mStart = NSMemoryConsumation::memoryUsage("Before", "IRGlobal", GC);
		targetSeedSet = NSIR::callIRGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, 0.7, curveInfoISV, curveInfoRTV, indexPlot++);
		mFinish = NSMemoryConsumation::memoryUsage("After", "IRGlobal", GC);
		printf("\nIRGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime5.GetTmStr());
		fprintf(F, "\nIRGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime5.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
*/
			

	
		
	
	TExeTm ExeTime8;
	printf("\nCELFGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "CELFGlobal", GC);
	targetSeedSet = NScelf::callCELFGlobal(Graph,GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "CELFGlobal", GC);
	printf("\nCELFGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime8.GetTmStr());
	fprintf(F, "\nCELFGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime8.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

//
//	//TExeTm ExeTime4;
//	//printf("\nKunduRankGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
//	//mStart = NSMemoryConsumation::memoryUsage("Before", "KunduRankGlobal", GC);
//	//targetSeedSet = NSbaseLineRank::callKunduRankGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
//	//mFinish = NSMemoryConsumation::memoryUsage("After", "KunduRankGlobal", GC);
//	//printf("\nKunduRankGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime4.GetTmStr());
//	//fprintf(F, "\nKunduRankGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime4.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
//
//
//	//TExeTm ExeTime5;
//	//printf("\nIRGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
//	//mStart = NSMemoryConsumation::memoryUsage("Before", "IRGlobal", GC);
//	//targetSeedSet = NSIR::callIRGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, 0.7, curveInfoISV, curveInfoRTV, indexPlot++);
//	//mFinish = NSMemoryConsumation::memoryUsage("After", "IRGlobal", GC);
//	//printf("\nIRGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime5.GetTmStr());
//	//fprintf(F, "\nIRGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime5.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
//
//
//	
//
//	/*TExeTm ExeTime7;
//	printf("\nIRIEPMIAGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
//	mStart = NSMemoryConsumation::memoryUsage("Before", "IRIEPMIAGlobal", GC);
//	targetSeedSet = NSIRIE::callIRIEPMIAGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, 0.7, curveInfoISV, curveInfoRTV, indexPlot++);
//	mFinish = NSMemoryConsumation::memoryUsage("After", "IRIEPMIAGlobal", GC);
//	printf("\nIRIEPMIAGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime7.GetTmStr());
//	fprintf(F, "\nIRIEPMIAGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime7.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
//	*/
//
//	//TExeTm ExeTime8;
//	//printf("\nCELFGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
//	//mStart = NSMemoryConsumation::memoryUsage("Before", "CELFGlobal", GC);
//	//targetSeedSet = NScelf::callCELFGlobal(Graph,GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
//	//mFinish = NSMemoryConsumation::memoryUsage("After", "CELFGlobal", GC);
//	//printf("\nCELFGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime8.GetTmStr());
//	//fprintf(F, "\nCELFGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime8.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
//
//
//	
//	//TExeTm ExeTime9;
//	//printf("\nCELFPPPGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
//	//mStart = NSMemoryConsumation::memoryUsage("Before", "CELFPPPGlobal", GC);
//	//targetSeedSet = NScelf::callCELFPPPGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
//	//mFinish = NSMemoryConsumation::memoryUsage("After", "CELFPPPGlobal", GC);
//	//printf("\nCELFPPPGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime9.GetTmStr());
//	//fprintf(F, "\nCELFPPPGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime9.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
//	/*
//
//	TExeTm ExeTime10;
//	printf("\nGreedyGlobal is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
//	mStart = NSMemoryConsumation::memoryUsage("Before", "GreedyGlobal", GC);
//	targetSeedSet = NSgreedy::callGreedyGlobal(Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
//	mFinish = NSMemoryConsumation::memoryUsage("After", "GreedyGlobal", GC);
//	printf("\nGreedyGlobal was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime10.GetTmStr());
//	fprintf(F, "\nGreedyGlobal Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime10.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
//*/

}


