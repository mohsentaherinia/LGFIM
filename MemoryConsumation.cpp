#include "stdafx.h"
DWORD NSMemoryConsumation::FindProcessId(const std::string& processName)
{
	PROCESSENTRY32 processInfo;
	processInfo.dwSize = sizeof(processInfo);

	HANDLE processesSnapshot = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, NULL);
	if (processesSnapshot == INVALID_HANDLE_VALUE)
		return 0;

	Process32First(processesSnapshot, &processInfo);
	if (!processName.compare(processInfo.szExeFile))
	{
		CloseHandle(processesSnapshot);
		return processInfo.th32ProcessID;
	}

	while (Process32Next(processesSnapshot, &processInfo))
	{
		if (!processName.compare(processInfo.szExeFile))
		{
			CloseHandle(processesSnapshot);
			return processInfo.th32ProcessID;
		}
	}

	CloseHandle(processesSnapshot);
	return 0;
}

int NSMemoryConsumation::func(std::string processName)
{

	DWORD processID = NSMemoryConsumation::FindProcessId(processName);

	if (processID == 0)
		std::cout << "Could not find " << processName.c_str() << std::endl;
	//else
		//std::cout << "Process Name is = " << processName << "Process ID is " << processID << std::endl;
	return processID;
}

NSMemoryConsumation::memStruct NSMemoryConsumation::memoryUsage(const string&  state, char * AlgName, const GlobalConst &GC)
{
	memStruct memOutput;
	FILE *F = fopen("__Memory.not", "at+");
	const char wtext[] = "MYIM.exe";
	string nn = wtext;
	DWORD processID = NSMemoryConsumation::func(nn);
	HANDLE hProcess = OpenProcess(
		PROCESS_QUERY_INFORMATION | PROCESS_VM_READ | SYNCHRONIZE,
		FALSE,
		processID);
	PROCESS_MEMORY_COUNTERS_EX pmc;
	ZeroMemory(&pmc, sizeof(PROCESS_MEMORY_COUNTERS_EX));
	ZeroMemory(&pmc, sizeof(PROCESS_MEMORY_COUNTERS_EX));
	GetProcessMemoryInfo(hProcess, (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	memOutput.PeakWorkingSetSize = ((float)pmc.PeakWorkingSetSize / 1024.0) / 1024.0;
	memOutput.PrivateUsage = ((float)pmc.PrivateUsage / 1024.0) / 1024.0;
	
	//PeakWorkingSetSize is equal to VS Private Bytes Window
	if ((state=="Before")) {
		fprintf(F, "Allocated Memory \"Before\" Executing the algorithm  \"%s\" on Dataset \"%s\" DF=%s, K=%d\n", AlgName, GC._inputFileName, GC._DifModel, GC._seedSizeConst);
		fprintf(F, "Before:\n");
		fprintf(F, "PeakWorkingSetSize--(MB) : %f\n", memOutput.PeakWorkingSetSize);
		fprintf(F, "PrivateUsage--------(MB) : %f\n", memOutput.PrivateUsage);
	}
	if ((state== "After")) {
		fprintf(F, "After:\n");
		fprintf(F, "PeakWorkingSetSize--(MB) : %f\n", memOutput.PeakWorkingSetSize);
		fprintf(F, "PrivateUsage--------(MB) : %f\n\n", memOutput.PrivateUsage);
	}
	if ((state=="Initialize")) {
		fprintf(F, "\n\n\n************************************************************************************************\n");
		fprintf(F, "Allocated Memory After Initilaize Dataset = %s\n", GC._inputFileName);
		fprintf(F, "PeakWorkingSetSize--(MB) : %f\n", memOutput.PeakWorkingSetSize);
		fprintf(F, "PrivateUsage--------(MB) : %f\n\n", memOutput.PrivateUsage);
	}
	fflush(F);
	CloseHandle(hProcess);
	fclose(F);
	return memOutput;
	//system("pause");
}
