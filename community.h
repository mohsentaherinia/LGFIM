#pragma once
#include "stdafx.h"
namespace NScommunity {
	void SaveCommunityInfoTofile(const TCnComV & CmtyV, const string & DataSetName, const string & CmtyAlg);
	vector<int>* readCommunityInfoFromFile(const char * DataSetName, const char * CmtyAlg, int & cmtyCount);
	void readLargeCommunityNodesFromFile(TIntV & NIdV);

	template<class PGraph>
	double CommunityDetectionRun(const PGraph &GraphORG, const char*DataSetName, const string& CmtyAlg) {
		PUNGraph Graph = TSnap::ConvertGraph<PUNGraph, PNGraph>(GraphORG);
		TExeTm ExeTm;
		printf("\nThe \"%s\" Community Detection algorithm Start at: %s (%s).\nPlease Wait(TOO SLOW!!!)", CmtyAlg.c_str(), ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
		//Try
		TSnap::DelSelfEdges(Graph);
		TCnComV CmtyV;
		double Q = 0.0;
		if (CmtyAlg == "Girvan-Newman") {
				Q = TSnap::CommunityGirvanNewman(Graph, CmtyV);
		}
		else if (CmtyAlg == "Clauset-Newman-Moore") {
			FILE *F = fopen("_CommunityDtectionInfo.not", "at+");
			fprintf(F, "\n\nCommunity Detection algorithm CNM\nDataset Name = %s\nStarting Time of CNM = %s",DataSetName, TExeTm::GetCurTm()); fflush(F);
			Q = TSnap::CommunityCNM(Graph, CmtyV);
			fprintf(F, "\nEnd Time of CNM =%s ,Modularity=%Lf", TExeTm::GetCurTm(), Q); fflush(F);
			fclose(F);
			cout <<"\nCNM Q= " << Q << endl;
			/*Q = NSlpa::GetModularityLescovich(Graph, CmtyV);
			cout << "Lescovich CNM Q= " << Q << endl;
			Q = NSlpa::GetModularity2(Graph, CmtyV);
			cout << "Khodam CNM   Q= " << Q << endl;*/
		}
		else if (CmtyAlg == "Infomap") {
			Q = TSnap::Infomap(Graph, CmtyV);
		}
		else if (CmtyAlg == "CLPA") {
			Q = NSclpa::callCLPA(GraphORG, CmtyV, DataSetName);//GraphORG GraphORG GraphORG GraphORG GraphORG
		}
		else {
			Fail;
		}
		printf("\nrun time Community Detection Algorithm: %s (%s)\n", ExeTm.GetTmStr(), TSecTm::GetCurTm().GetTmStr().CStr());
		SaveCommunityInfoTofileDetails(GraphORG, CmtyV, DataSetName, CmtyAlg,Q);
		SaveCommunityInfoTofile(CmtyV, DataSetName, CmtyAlg);
		//Catch
		return Q;
	}


	//Return Subgraphs 
	//Read the file contain of the extracted community
	template<class PGraph>
	PNGraph* ObtainExtractedCommunityFromFile(const PGraph &Graph, const char*DataSetName, const char* CmtyAlg, int &CmtyCount, double &Q) {
		vector<int> *CmtyVec;
		CmtyVec = readCommunityInfoFromFile(DataSetName + 3, CmtyAlg, CmtyCount);

		PNGraph *SubG = new PNGraph[CmtyCount];
		for (int c = 0; c < CmtyCount; c++) {
			TIntV NIdV;
			//	printf("\nCommunity %d => ", c + 1);
			for (int i = 0; i < CmtyVec[c].size(); i++) {
				NIdV.Add(CmtyVec[c][i]);
			}
			SubG[c] = TSnap::GetSubGraph(Graph, NIdV);
			//printf("nodes:%d  edges:%d\n", SubG[c]->GetNodes(), SubG[c]->GetEdges());
		}
		return SubG;
	}

	template<class PGraph>
	void SaveCommunityInfoTofileDetails(const PGraph &Graph, const TCnComV &CmtyV, const string& DataSetName, const string& CmtyAlg, const double &Q) {
		char  Filename[80] = "Z_CD_Info_";
		strcat(Filename, DataSetName.c_str());
		strcat(Filename, "_");
		strcat(Filename, CmtyAlg.c_str());
		strcat(Filename, ".txt");

		FILE *F = fopen(Filename, "wt+");
		fprintf(F, "# Nodes: %d    Edges: %d\n", Graph->GetNodes(), Graph->GetEdges());
		fprintf(F, "# Algoritm: %s\n", CmtyAlg.c_str());
		if (CmtyAlg != "Infomap") {
			fprintf(F, "# Modularity: %f\n", Q);
		}
		else {
			fprintf(F, "# Average code length: %f\n", Q);
		}
		fprintf(F, "# Communities: %d\n", CmtyV.Len()); 
			for (int c = 0; c < CmtyV.Len(); c++)
				fprintf(F, "cmty # %d is composed of %d\n", c, CmtyV[c].Len());
		//fprintf(F, "%d\n",CmtyV[c].Len());
		fprintf(F, "# NId\tCommunityId\n");
		for (int c = 0; c < CmtyV.Len(); c++) {
			for (int i = 0; i < CmtyV[c].Len(); i++) {
				fprintf(F, "%d\t%d\n", CmtyV[c][i].Val, c);
			}
		}
		fclose(F);
	}

	template<class PGraph>
	void constructSubGraphOfLargeCommunityFromFile(const PGraph &Graph,PGraph &subGraph) {
		TIntV NIdV;
		NScommunity::readLargeCommunityNodesFromFile(NIdV);
		subGraph = TSnap::GetSubGraph(Graph, NIdV);
		printf("\nNew Sub Graghnodes:%d  edges:%d\n", subGraph->GetNodes(), subGraph->GetEdges());

	}


};