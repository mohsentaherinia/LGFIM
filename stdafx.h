#pragma once

#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include <list>
#include <array>
#include <set>
#include <omp.h>
#include <fstream>
#include <string>
#include <algorithm>


struct GlobalConst {
	int _MCConst;
	int _MCGreedyconst;
	double _ICProbConsttemp;
	int _seedSizeConst;
	char* _inputFileName;
	int _NodeDivision;
	int _coef1;
	float _coef2;
	char* _DifModel;
	char* _ComDetMethod;
	int _SixDegree;
	int _simMode;
};





using namespace std;

#include "targetver.h"
#include "Snap.h"

#include "ICModel.h"
#include "ICModelParallel.h"
#include "ICModelParallelFull.h"
#include "LTModel.h"
#include "LTModelParallel.h"
#include "LTModelParallelFull.h"
#include "WCModel.h"
#include "WCModelParallel.h"
#include "WCModelParallelFull.h"

#include "Tools.h"
#include "Heap.h"


#include "BaseLineRank.h"
#include "DegreeDiscount.h"
#include "Greedy.h"
#include "Celf.h"
#include "IR.h"
#include "IRIE.h"
#include "SimPath.h"
#include "community.h"
#include "LPA.h"
#include "CLPA.h"
#include "CGA.h"
#include "m_IMCD.h"
#include "MyLAIM.h"
#include "MemoryConsumation.h"
#include "COFIM.h"
#include "NewIRIE.h"



