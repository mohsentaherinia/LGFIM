#pragma once
#include "stdafx.h"


template<class PNGraph>
vector<int> MiningLocally(const PNGraph &SubG, const char* ModelName,const int &N, const GlobalConst &GC,const int &coef1,const float &coef2)
{
	vector<int> CurrentSeedSet;
	int local_i_k = (coef1 * GC._seedSizeConst * SubG->GetNodes() / N) + 1;
	if (local_i_k < 0)//overflow for big dataset
		local_i_k = 50;
	local_i_k = min(local_i_k, GC._seedSizeConst);
	

	if ((float)SubG->GetNodes() / (float)N <= coef2 || (SubG->GetNodes() <= local_i_k))
		return CurrentSeedSet;
	
	//cout << "\nKi=" << local_i_k;
	if (SubG->GetEdges()==0)
		return CurrentSeedSet;
	CurrentSeedSet = NSIRIE::callIRIEPMIALocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCConst, 0.7);
	//CurrentSeedSet = NSMyLAIM::callMyLAIMLocal(SubG,GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst,GC._SixDegree,0.001);
	/*CurrentSeedSet = NSDegreeDiscount::callDegDisLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst,0.01);
	CurrentSeedSet = NSIRIE::callIRIEPMIALocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCConst, 0.7);
	CurrentSeedSet = NSbaseLineRank::callMaxDegreeRankLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NSbaseLineRank::callPageRank2Local(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NSbaseLineRank::callPageRank1Local(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NSbaseLineRank::callKunduRankLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NSbaseLineRank::callMinDegreeRankLocal(SubG,GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NSIR::callIRLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst,0.7);
	CurrentSeedSet = NSgreedy::callGreedyLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NScelf::callCELFPPPLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	CurrentSeedSet = NScelf::callCELFLocal(SubG, GC, ModelName, local_i_k, GC._ICProbConsttemp, GC._MCGreedyconst);
	*/return CurrentSeedSet;
}

template<class PNGraph>
void CDIM(const PNGraph &Graph, const char* ModelName, const char * CDalgName, const GlobalConst &GC, const int &coef1, const float &coef2, vector<NSplot::curveInfo> &curveInfoISV, vector<NSplot::curveInfo> &curveInfoRTV, int &indexPlot) {
	int cmtyNumber;
	PNGraph *SubG;
	double Q = 0.0;
	int N = Graph->GetNodes();
	SubG = NScommunity::ObtainExtractedCommunityFromFile(Graph, GC._inputFileName, CDalgName, cmtyNumber, Q);
	TExeTm ExeTm0;
	printf("\nAll of the communities are exploring locally: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	vector<int> CandidateSet0; 
	vector<int> CandidateSet1;
	vector<int> CandidateSet2;
	vector<int> CandidateSet3;
	vector<int> CandidateSet4;
	vector<int> CandidateSet5;
	vector<int> CandidateSet6;
	vector<int> CandidateSet7;
		
#pragma omp parallel sections
	{
#pragma omp section
	{
		for (size_t i = 0; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet0.insert(CandidateSet0.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 1; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet1.insert(CandidateSet1.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 2; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet2.insert(CandidateSet2.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 3; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet3.insert(CandidateSet3.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 4; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet4.insert(CandidateSet4.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 5; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet5.insert(CandidateSet5.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 6; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet6.insert(CandidateSet6.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
#pragma omp section
	{
		for (size_t i = 7; i < cmtyNumber; i += 8)
		{
			vector<int> TempSeedSet = MiningLocally(SubG[i], ModelName, N, GC,coef1,coef2);
			CandidateSet7.insert(CandidateSet7.begin(), TempSeedSet.begin(), TempSeedSet.end());
		}
	}
	}
	
	vector<int> CandidateSet;
	CandidateSet.insert(CandidateSet.begin(), CandidateSet0.begin(), CandidateSet0.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet1.begin(), CandidateSet1.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet2.begin(), CandidateSet2.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet3.begin(), CandidateSet3.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet4.begin(), CandidateSet4.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet5.begin(), CandidateSet5.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet6.begin(), CandidateSet6.end());
	CandidateSet.insert(CandidateSet.begin(), CandidateSet7.begin(), CandidateSet7.end());
	if (CandidateSet.size() < GC._seedSizeConst) {
		FILE *F = fopen("__Info.not", "at+");
		fprintf(F, "\n\n\nthe size of the candidate set smaller than from \"k\"\n");
		fprintf(F, "the size of the candidate set smaller than from \"k\"\n");
		fprintf(F, "the size of the candidate set smaller than from \"k\"\n");
		fprintf(F, "Local Exploring failure...\n");
		fclose(F);
		cout << "the size of the candidate set smaller than from \"k\"\n";
		return;
	}
	printf("\nLocal exploration of communities was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTm0.GetTmStr());
	FILE *F = fopen("__Info.not", "at+");
	fprintf(F, "\nLocal exploration Running Time = %s", ExeTm0.GetTmStr()); fflush(F);

	vector<int> CurrentSeedSet;
	NSMemoryConsumation::memStruct mStart, mFinish;


	FILE *FRT = fopen("__Compare_RT_IRIEGLwithCELF.not", "at+");
	fprintf(FRT, "\n\n\t\tIRIE GLobal Of local\nDatasetName = %s\t#Nodes = %d\t#Edges = %d\tcoef1 = %d\tcoef2 = %Lf\n", GC._inputFileName, Graph->GetNodes(), Graph->GetEdges(), GC._coef1, GC._coef2);
	fclose(FRT);

	TExeTm ExeTime50;
	printf("\nIRIEPMIAGL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "IRIEPMIAGL", GC);
	CurrentSeedSet = NSIRIE::callIRIEPMIAGlobalOfLocal(CandidateSet, Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, 0.7, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "IRIEPMIAGL", GC);
	printf("\nIRIEPMIAGL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime50.GetTmStr());
	fprintf(F, "\nIRIEPMIAGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime50.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);


	/*TExeTm ExeTime500;
	printf("\nCofimCELF (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "CofimCELF", GC);
	CurrentSeedSet = NSCOFIM::callCofimCELFGlobal(SubG, cmtyNumber, Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "CofimCELF", GC);
	printf("\nCofimCELF (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime500.GetTmStr());
	fprintf(F, "\nCofimCELF Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime500.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	
	TExeTm ExeTime5000;
	printf("\nCofimGreedy (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "CofimGreedy", GC);
	CurrentSeedSet = NSCOFIM::callCofimGreedyGlobal(SubG, cmtyNumber, Graph, GC, GC._DifModel, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "CofimGreedy", GC);
	printf("\nCofimGreedy (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime5000.GetTmStr());
	fprintf(F, "\nCofimGreedy Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime5000.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
*/
	/*TExeTm ExeTime7;
	printf("\nCELFGL is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "CELFGL", GC);
	CurrentSeedSet = NScelf::callCELFGlobalOfLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "CELFGL", GC);
	printf("\nCELFGL (using CD) was finished  at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime7.GetTmStr());
	fprintf(F, "\nCELFGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime7.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
*/
	
	//TExeTm ExeTime111;
	//printf("\nMyLAIMGL1 is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	//mStart = NSMemoryConsumation::memoryUsage("Before", "MyLAIMGL1", GC);
	//CurrentSeedSet = NSMyLAIM::callMyLAIMGlobalOfLocal1(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++,GC._SixDegree,0.001);
	//mFinish = NSMemoryConsumation::memoryUsage("After", "MyLAIMGL1", GC);
	//printf("\nMyLAIMGL1 (using CD) was finished  at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime111.GetTmStr());
	//fprintf(F, "\nMyLAIMGL1 Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime111.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	/*TExeTm ExeTime112;
	printf("\nMyLAIMGL2 is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "MyLAIMGL2", GC);
	CurrentSeedSet = NSMyLAIM::callMyLAIMGlobalOfLocal2(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++, GC._SixDegree, 0.001);
	mFinish = NSMemoryConsumation::memoryUsage("After", "MyLAIMGL2", GC);
	printf("\nMyLAIMGL2 (using CD) was finished  at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime112.GetTmStr());
	fprintf(F, "\nMyLAIMGL2 Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime112.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
	*/
	/*TExeTm ExeTime2;
	printf("\nMaxDegreeRankGL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "MaxDegreeRankGL", GC);
	CurrentSeedSet = NSbaseLineRank::callMaxDegreeRankGlobalofLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "MaxDegreeRankGL", GC);
	printf("\nMaxDegreeRankGL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime2.GetTmStr());
	fprintf(F, "\nMaxDegreeRankGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime2.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
*/
	
	/*TExeTm ExeTime1;
	printf("\nDegDisGL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "DegDisGL", GC);
	CurrentSeedSet = NSDegreeDiscount::callDegDisGlobalOfLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst,0.000001, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "DegDisGL", GC);
	printf("\nDegDisGL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime1.GetTmStr());
	fprintf(F, "\nDegDisGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime1.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);
	*/

	/*
	TExeTm ExeTime2;
	printf("\nMaxDegreeRankGL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "MaxDegreeRankGL", GC);
	CurrentSeedSet = NSbaseLineRank::callMaxDegreeRankGlobalofLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "MaxDegreeRankGL", GC);
	printf("\nMaxDegreeRankGL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime2.GetTmStr());
	fprintf(F, "\nMaxDegreeRankGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime2.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	TExeTm ExeTime3;
	printf("\nPageRank2GL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "PageRank2GL", GC);
	CurrentSeedSet = NSbaseLineRank::callPageRank2GlobalofLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "PageRank2GL", GC);
	printf("\nPageRank2GL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime3.GetTmStr());
	fprintf(F, "\nPageRank2GL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime3.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);


	TExeTm ExeTime4;
	printf("\nKunduRankGL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "KunduRankGL", GC);
	CurrentSeedSet = NSbaseLineRank::callKunduRankGlobalofLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "KunduRankGL", GC);
	printf("\nKunduRankGL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime4.GetTmStr());
	fprintf(F, "\nKunduRankGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime4.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);


	TExeTm ExeTime5;
	printf("\nIRGL (using CD) is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "IRGL", GC);
	CurrentSeedSet = NSIR::callIRGlobalofLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCConst,0.7, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "IRGL", GC);
	printf("\nIRGL (using CD) was finished at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime5.GetTmStr());
	fprintf(F, "\nIRGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime5.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	
	TExeTm ExeTime6;
	printf("\nCELFPPPGL is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "CELFPPPGL", GC);
	CurrentSeedSet = NScelf::callCELFPPPGlobalOfLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "CELFPPPGL", GC);
	printf("\nCELFPPPGL (using CD) was finished  at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime6.GetTmStr());
	fprintf(F, "\nCELFPPPGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime6.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	
	TExeTm ExeTime8;	
	printf("\nGreedyGL is starting at now: %s\n", TSecTm::GetCurTm().GetTmStr().CStr());
	mStart = NSMemoryConsumation::memoryUsage("Before", "GreedyGL", GC);
	CurrentSeedSet = NSgreedy::callGreedyGlobalOfLocal(CandidateSet, Graph, GC, ModelName, GC._seedSizeConst, GC._ICProbConsttemp, GC._MCGreedyconst, curveInfoISV, curveInfoRTV, indexPlot++);
	mFinish = NSMemoryConsumation::memoryUsage("After", "GreedyGL", GC);
	printf("\nGreedyGL (using CD) was finished  at: %s (Run Time = %s)\n", TSecTm::GetCurTm().GetTmStr().CStr(), ExeTime8.GetTmStr());
	fprintf(F, "\nGreedyGL Running Time = %s, PeakWorkingSetSize = %f (MB), PrivateUsage = %f (MB)", ExeTime8.GetTmStr(), mFinish.PeakWorkingSetSize - mStart.PeakWorkingSetSize, mFinish.PrivateUsage - mStart.PrivateUsage); fflush(F);

	*/

	
	fclose(F);
}


