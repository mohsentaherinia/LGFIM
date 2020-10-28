#pragma once
#include "stdafx.h"
#include <stdio.h>
#include <tchar.h>
#include <psapi.h>
#include <iostream>
#include <windows.h>
#include <tlhelp32.h>
using namespace std;
namespace NSMemoryConsumation {

	struct memStruct {
		float PeakWorkingSetSize, PrivateUsage;
	};

	DWORD FindProcessId(const std::string & processName);
	int func(std::string processName);
	memStruct memoryUsage(const string & state, char * AlgName, const GlobalConst & GC);
	

};


