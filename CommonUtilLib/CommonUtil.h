#pragma once

#include <vector>
using namespace std;

namespace Another_CommonUtil
{

	class CommonUtil
	{
	public:
		static void NormalizeScalar(const vector<double>& scalars, vector<double>& normalizeResult);
		static double FindMaxium(const vector<double>& scalars);
		static double FindMinum(vector<double> scalars);
		static void AcTan(const vector<double>& scalars, vector<double>& result);
		static void Filter(const vector<double>&scalars, vector<double>& result, int n);
	};
}