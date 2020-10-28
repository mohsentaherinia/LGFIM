#include "stdafx.h"

void NSCGA::integrate(const int& cm, const int& cl, TCnComV& CmtyV) {
	for (size_t i = 0; i < CmtyV[cl].Len(); i++)
		CmtyV[cm].Add(CmtyV[cl][i]);
	CmtyV[cl].Clr();
}

bool NSCGA::IsexistVertexInL(const int & u, const set<int>L[], const int & cm) {
	for (set<int>::iterator it = L[cm].begin(); it != L[cm].end(); ++it)
		if (*it == u)
			return true;
	return false;
}

void NSCGA::DifferenceVector(vector<int>& S1, vector<int>& S2, vector<int>& S3) {
	S3.resize(S1.size());
	vector<int>::iterator it;
	std::sort(S1.begin(), S1.end());
	std::sort(S2.begin(), S2.end());

	it = std::set_difference(S1.begin(), S1.end(), S2.begin(), S2.end(), S3.begin());
	S3.resize(it - S3.begin());
}

int NSCGA::getCmtyLAbelOfvertex(const TCnComV& CmtyV, const int& vertex) {
	for (int c = 0; c < CmtyV.Len(); c++) {
		for (int i = 0; i < CmtyV[c].Len(); i++) {
			if (CmtyV[c][i] == vertex)
				return c;
		}
	}
	return -1;
}

