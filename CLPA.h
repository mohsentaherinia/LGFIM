#pragma once
#include "stdafx.h"
#define ConvergeRep 50
#define K_coeff 49
#define alphaConst 0.05
#define betaConst 0.05
#define mul 1


namespace NSclpa {
	int CummunityCapacity(const int & t, const int & N);

	int findMaximumLabelRandom(TIntIntH  FvL, vector<int> A);

	int findMaximumLabelRandom2(TIntIntH & FvL, vector<int>& A);

	template<class PGraph>
	void calcLabelFreqOfNbrsV(const PGraph&Graph, const int & v, TIntIntH& FvL, const TIntIntH &labelV) {
		PGraph::TObj::TNodeI NI = Graph->GetNI(v);
		for (int e = 0; e < NI.GetOutDeg(); e++) {
			int neighbourID = NI.GetOutNId(e);
			int label = labelV.GetDat(neighbourID);
			if (FvL.IsKey(label)) {
				FvL.AddDat(label) = FvL.GetDat(label) + 1; //++
			}
			else {
				FvL.AddDat(label, 1);
			}
		}
	}

	template<class PGraph>
	void internalFunction(const PGraph& Graph, const int& v, TIntIntH &labelV, TIntIntH& countLabelL, const int& cap, const int &P) {
		vector <int> A;
		TIntIntH FvL;
		calcLabelFreqOfNbrsV(Graph, v, FvL, labelV);
		for (size_t index = 0; index < FvL.Len(); index++) {
			int Label = FvL.GetKey(index);
			if (FvL.GetDat(Label) > 0)
				if (countLabelL.GetDat(Label) < cap)
					A.push_back(Label);
		}
		int Candidate_v_Label = findMaximumLabelRandom2(FvL, A);
		if (Candidate_v_Label == -1)
			return;
		double rnd = ((double)rand() / (RAND_MAX));
		if (rnd < P) {
			labelV.AddDat(v) = Candidate_v_Label;
			countLabelL.AddDat(Candidate_v_Label) = countLabelL.GetDat(Candidate_v_Label) + 1;
		}
		else
		{
			int Fvlv = 0;
			int current_v_label = labelV.GetDat(v);
			if (FvL.IsKey(current_v_label))
				Fvlv = FvL.GetDat(current_v_label);
			if (countLabelL.GetDat(Candidate_v_Label) > Fvlv) {
				labelV.AddDat(v) = Candidate_v_Label;
				countLabelL.AddDat(Candidate_v_Label) = countLabelL.GetDat(Candidate_v_Label) + 1;
			}
		}
	}

	template<class PGraph>
	void internalLoop(const PGraph& Graph, vector<int>& setVertex, TIntIntH &labelV, TIntIntH& countLabelL, const int& cap, const int &P) {
		const int div = 8;
#pragma omp parallel sections
		{
#pragma omp section
		{
			for (size_t index = 0 * setVertex.size() / div; index < 1 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 1 * setVertex.size() / div; index < 2 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 2 * setVertex.size() / div; index < 3 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 3 * setVertex.size() / div; index < 4 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 4 * setVertex.size() / div; index < 5 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 5 * setVertex.size() / div; index < 6 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 6 * setVertex.size() / div; index < 7 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}
#pragma omp section
		{
			for (size_t index = 7 * setVertex.size() / div; index < 8 * setVertex.size() / div; index++)
				internalFunction(Graph, setVertex[index], labelV, countLabelL, cap, P);
		}

		}
	}

	template<class PGraph>
	void combinationStep2(const PGraph& Graph, TCnComV &CmtyV, const double& alpha, const double &beta) {
		int N = Graph->GetNodes();
		set <int> cmtysFull, cmtysNotFull;//, cmtysCandidate;
		TIntIntH canSetFreqV;
		TIntIntH cmtyLabelNode;
		//traverse all node in all community and get their labels
		for (int i = 0; i < CmtyV.Len(); i++)
			for (size_t j = 0; j < CmtyV[i].Len(); j++)
				cmtyLabelNode.AddDat(CmtyV[i][j]);
		//traverse all community and classify them based on their lenght
		for (int c = 0; c < CmtyV.Len(); c++) {
			if (CmtyV[c].Len()>alpha*N)
				cmtysFull.insert(c);
			else
				cmtysNotFull.insert(c);
			//if (CmtyV[c].Len() < beta*N) {
			if (CmtyV[c].Len() < beta*N || CmtyV[c].Len() <= 10) {
				//cmtysCandidate.insert(c);
				canSetFreqV.AddDat(c) = CmtyV[c].Len();
			}
		}
		canSetFreqV.SortByDat(true);//asc

									//pick up one by one cmty from candidateSet and calculate label of nbr vertex . then select cmty with maximum labal that isExist in cmtynotfull set. then merge and modify all variable
		for (size_t i = 0; i < canSetFreqV.Len(); i++) {
			int candCmty = canSetFreqV.GetKey(i);
			if (CmtyV[candCmty].Len() == 0)
				continue;
			TIntIntH nbrLabelCountV;
			//calculate frequence community label in nbr of current community candCmty
			for (size_t k = 0; k < CmtyV[candCmty].Len(); k++)
			{
				int vertex = CmtyV[candCmty][k];
				PGraph::TObj::TNodeI NI = Graph->GetNI(vertex);
				for (int e = 0; e < NI.GetDeg(); e++) {
					int neighbourID = NI.GetNbrNId(e);
					int Nbrlabelcmty = cmtyLabelNode.GetDat(neighbourID);
					if (nbrLabelCountV.IsKey(Nbrlabelcmty))
						nbrLabelCountV.AddDat(Nbrlabelcmty) = nbrLabelCountV.GetDat(Nbrlabelcmty) + 1;//++
					else
						nbrLabelCountV.AddDat(Nbrlabelcmty) = 1;
				}
			}
			//find maximun frequence community label in nbr of current community candCmty and merge
			nbrLabelCountV.SortByDat(false);//decending
			for (size_t index = 0; index < nbrLabelCountV.Len(); index++) {
				int newlabelcmty = nbrLabelCountV.GetKey(index);
				if (newlabelcmty == candCmty)
					continue;
				else {
					bool isExist = cmtysNotFull.find(newlabelcmty) != cmtysNotFull.end();
					if (isExist == false)
						continue;
					//join candCmty , newlabelcmty
					//modify CmtyV ,cmtyLabelNode
					for (size_t cv = 0; cv < CmtyV[candCmty].Len(); cv++) {
						CmtyV[newlabelcmty].Add(CmtyV[candCmty][cv]);
						cmtyLabelNode.AddDat(CmtyV[candCmty][cv]) = newlabelcmty;
					}
					CmtyV[candCmty].Clr();
					cmtysNotFull.erase(candCmty);
					canSetFreqV.DelKey(candCmty);
					canSetFreqV.Defrag();
					//
					if (CmtyV[newlabelcmty].Len()>alpha*N) {
						cmtysFull.insert(newlabelcmty);
						cmtysNotFull.erase(newlabelcmty);
					}
					if (CmtyV[newlabelcmty].Len() > beta*N)
						if (canSetFreqV.IsKey(newlabelcmty)) {
							canSetFreqV.DelKey(newlabelcmty);
							canSetFreqV.Defrag();
						}
					//
					break;
				}
			}
		}
		// delete all empty cmty
		TCnComV CmtyVTemp;
		for (int i = 0; i < CmtyV.Len(); i++)
			if (CmtyV[i].Len() != 0)
				CmtyVTemp.Add(CmtyV[i]);
		CmtyV = CmtyVTemp;
		//		
	}

	template<class PGraph>
	double interCLPA(const PGraph& Graph, TCnComV &CmtyV, const char* datasetName) {
		FILE *F = fopen("_CommunityDtectionInfo.not", "at+");
		fprintf(F, "\n\nCommunity Detection algorithm = %s\nDataset Name = %s\nConvergeRep = %d\nK_coeff = %d\nMul = %d\nAlphaConst = %Lf\nBetaConst = %Lf\nStarting Time of traditional CLPA = %s","CLPA", datasetName, ConvergeRep,K_coeff, mul, alphaConst, betaConst,TExeTm::GetCurTm()); fflush(F);
		vector <int> setVertex;
		TIntIntH labelV, countLabelL;
		for (PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
			int v = NI.GetId();
			labelV.AddDat(v, v);
			countLabelL.AddDat(v, 1);
			setVertex.push_back(v);
		}
		////////////////////////////////////////////////////////////
		for (size_t t = 1; t <= ConvergeRep; t++) {
			int cap = CummunityCapacity(t, Graph->GetNodes());
			cout << "\n" << t << "\t" << cap;
			double P = 1.0 / pow(t, 2.0);
			for (size_t kkk = 1; kkk <= K_coeff*mul; kkk++) {
				internalLoop(Graph, setVertex, labelV, countLabelL, cap, P);
			}
		}
		double Q = 0;
		NSlpa::ExtractCommunity2(Graph, labelV, CmtyV);
		//////////////////////////////////////////////////
		Q = NSlpa::GetModularityLescovich(Graph, CmtyV);
		NScommunity::SaveCommunityInfoTofileDetails(Graph, CmtyV, datasetName, "CLPA-BeforeMerging", Q);
		NScommunity::SaveCommunityInfoTofile(CmtyV, datasetName, "CLPA-BeforeMerging");
		fprintf(F, "\nEnd Time of traditional CLPA =%s ,Modularity=%Lf",TExeTm::GetCurTm(),Q); fflush(F);
		int cmtyLenght;
		do {
			cout << "\ncombinationStep";
			fprintf(F, "\nStart Time of Combination =%s", TExeTm::GetCurTm()); fflush(F);
			cmtyLenght = CmtyV.Len();
			combinationStep2(Graph, CmtyV, alphaConst, betaConst);
			Q = NSlpa::GetModularityLescovich(Graph, CmtyV);
			NScommunity::SaveCommunityInfoTofileDetails(Graph, CmtyV, datasetName, "CLPA-AfterMerging", Q);
			NScommunity::SaveCommunityInfoTofile(CmtyV, datasetName, "CLPA-AfterMerging");
			fprintf(F, "\nEnd Time of Combination =%s ,Modularity=%Lf", TExeTm::GetCurTm(), Q); fflush(F);

		} while (cmtyLenght != CmtyV.Len());
		fclose(F);
		return Q;
	}

	template<class PGraph>
	double callCLPA(const PGraph& Graph, TCnComV &CmtyVG, const char* datasetName) {
		
		return NSclpa::interCLPA(Graph, CmtyVG, datasetName);
		//double Q = 0;
		//PNGraph SubG = Graph;
		//TCnComV tempCmtyV;
		//int kk = 1;
		//bool fullflag = true, notFullflag = true;
		//while (fullflag && notFullflag) {
		//	cout << "\n" << kk++;
		//	fullflag = false;
		//	notFullflag = false;
		//	tempCmtyV.Clr();
		//	NSclpa::interCLPA(SubG, tempCmtyV);
		//	int N = SubG->GetNodes();
		//	//create Graph with vertex of community that notFull
		//	TIntV NIdV;
		//	for (int c = 0; c < tempCmtyV.Len(); c++) {
		//		if (tempCmtyV[c].Len() <= alphaConst*N) { //Notfull
		//			notFullflag = true;
		//			for (size_t i = 0; i < tempCmtyV[c].Len(); i++) {
		//				NIdV.Add(tempCmtyV[c][i]);
		//			}
		//		}
		//		else { //is exist fuul community
		//			fullflag = true;
		//			CmtyVG.Add(tempCmtyV[c]);
		//		}
		//	}
		//	SubG = TSnap::GetSubGraph(SubG, NIdV);
		//	//
		//}
		//if (kk==2)
		//for (int c = 0; c < tempCmtyV.Len(); c++)
		//	CmtyVG.Add(tempCmtyV[c]);
		//return Q;
	}

}