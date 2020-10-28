#include "stdafx.h"

int NSIS::callInfluenceSpreadModel(const PNGraph &Graph, const vector<int>& seed, const string& Model, const double& ICProbb, const int& MCS) {
	int CountInfNode=-10000000;
	if (Model == "IC") {
		CountInfNode = NSICModel::ICinfluenceSpread(Graph, seed, ICProbb, MCS);
	}
	if (Model == "ICP") {
		CountInfNode = NSICModelP::ICinfluenceSpreadParallel4(Graph, seed, ICProbb, MCS);
	}
	if (Model == "ICPfull") {
		CountInfNode = NSICModelPF::ICinfluenceSpreadParallelFull(Graph, seed, ICProbb, MCS);
	}
	if (Model == "LT") {
		CountInfNode = NSLTModel::LTinfluenceSpread(Graph, seed, MCS);
	}
	if (Model == "LTP") {
		CountInfNode = NSLTModelP::LTinfluenceSpreadParallel4(Graph, seed, MCS);
	}
	if (Model == "LTPfull") {
		CountInfNode = NSLTModelPF::LTinfluenceSpreadParallelFull(Graph, seed, MCS);
	}
	if (Model == "WC") {
		CountInfNode = NSWCModel::WCinfluenceSpread(Graph, seed, MCS);
	}
	if (Model == "WCP") {
		CountInfNode = NSWCModelP::WCinfluenceSpreadParallel4(Graph, seed, MCS);
	}
	if (Model == "WCPfull") {
		CountInfNode = NSWCModelPF::WCinfluenceSpreadParallelFull(Graph, seed, MCS);
	}

	return CountInfNode;
}

bool NStools::findNodeInVector(vector<int>& seed, const int & Node)
	{
		for (vector<int>::iterator it = seed.begin(); it != seed.end(); ++it)
			if (*it == Node)
				return true;
		return false;
	}

void NStools::SaveToFile(char  *algName, const int& k, const int& infectedNodes, const int& GraphSize, const int& mc, const char* Model, const double& Runtime, const GlobalConst &GC)
	{
		char  Filename[80] = "_";
		strcat(Filename, GC._inputFileName + 4);
		strncpy_s(Filename, Filename, strlen(Filename) - 4);
		strcat(Filename, "_");
		strcat(Filename, Model);
		strcat(Filename, "_");
		strcat(Filename, algName);
		strcat(Filename, ".not");
		FILE *F = fopen(Filename, "at+");
		fprintf(F, "\n%s Algorithm\t%s DiffusionModel\t # %d MCSimulation\t ==>>> K= %d\tinfectedNOdes= %d\tPercent= %f\t%Lf Sec", algName, Model, mc, k, infectedNodes, (double)infectedNodes / GraphSize, Runtime);
		fclose(F);
	}


//void NStools::SaveToFile(char  *algName, const int& k, const int& infectedNodes, const int& GraphSize, const int& mc, const char* Model, const char* Runtime, const GlobalConst &GC)
//{
//	char  Filename[80] = "_";
//	strcat(Filename, GC._inputFileName + 4);
//	strncpy_s(Filename, Filename, strlen(Filename) - 4);
//	strcat(Filename, "_");
//	strcat(Filename, Model);
//	strcat(Filename, "_");
//	strcat(Filename, algName);
//	strcat(Filename, ".not");
//	FILE *F = fopen(Filename, "at+");
//	fprintf(F, "\n%s Algorithm\t%s DiffusionModel\t # %d MCSimulation\t ==>>> K= %d\tinfectedNOdes= %d\tPercent= %f\t%s Sec", algName, Model, mc, k, infectedNodes, (double)infectedNodes / GraphSize, Runtime);
//	fclose(F);
//}

void NStools::SaveSeedSetToFile(char  *algName, const GlobalConst &GC, const vector<int>& SeedSet)
{
	char  Filename[80] = "_";
	strcat(Filename, GC._inputFileName + 4);
	strncpy_s(Filename, Filename, strlen(Filename) - 4);
	strcat(Filename, "_");
	strcat(Filename, GC._DifModel);
	strcat(Filename, "__SeedSetGloal");
	strcat(Filename, ".not");
	FILE *F = fopen(Filename, "at+");
	fprintf(F, "\n\nAlgorithm = %s \tDiffusionModel = %s\n", algName, GC._DifModel);
	for (size_t i = 0; i < SeedSet.size(); i++)
		fprintf(F, "%d,\t", SeedSet[i]);
	fclose(F);
}

void NStools::SaveSeedSetToFileLocal(char  *algName, const GlobalConst &GC, const vector<int>& SeedSet)
{
	char  Filename[80] = "_";
	strcat(Filename, GC._inputFileName + 4);
	strncpy_s(Filename, Filename, strlen(Filename) - 4);
	strcat(Filename, "_");
	strcat(Filename, GC._DifModel);
	strcat(Filename, "__SeedSetLocal");
	strcat(Filename, ".not");
	FILE *F = fopen(Filename, "at+");
	fprintf(F, "\n\nAlgorithm = %s \tDiffusionModel = %s\n", algName, GC._DifModel);
	for (size_t i = 0; i < SeedSet.size(); i++)
		fprintf(F, "%d,\t", SeedSet[i]);
	fclose(F);
}

void NSplot::plotAllInfluenceSpread(vector<curveInfo> &curveInfo)
{
	char  Filename[80] = "_";
	strcat(Filename, curveInfo[0].dataSet);
	strncpy_s(Filename, Filename, strlen(Filename) - 4);
	strcat(Filename, "_");
	strcat(Filename, curveInfo[0].algName);
	strcat(Filename, "__");
	strcat(Filename, curveInfo[0].typeCurve);
	

	//char  title[80] = "DifMod=";
	//strcat(title, curveInfo[0].algName);
	//strcat(title, "  DS=");
	char  title[80]="";
	strcat(title, curveInfo[0].dataSet);
	title[strlen(title) - 4] = '\0';

	TGnuPlot Gp(Filename, title);
	for (size_t i = 0; i < curveInfo.size(); i++)
	{
		Gp.AddPlot(curveInfo[i].XY, gpwLinesPoints, curveInfo[i].modelName);
	}
	Gp.SetXYLabel("k (Seed Set Size)", curveInfo[0].typeCurve);

	Gp.SavePng(); //or Gp.SaveEps();
				  ///open Picture in MS Windows////////////////////
	char fname[80];
	strcpy(fname, Filename);
	strcat(fname, ".png");
	system(fname);
}

void NSplot::setPlottingValue(vector<curveInfo> &curveInfoX, const int & index, const int & xAxis, const float &yAxis, const char* algName, const char* modelName, const char* typeCurve) {
	strcpy(curveInfoX[index].algName, algName);
	strcpy(curveInfoX[index].modelName, modelName);
	strcpy(curveInfoX[index].typeCurve, typeCurve);
	TFltPr xy;
	xy.Val1 = xAxis;
	xy.Val2 = yAxis;
	curveInfoX[index].XY.Add(xy);
}
