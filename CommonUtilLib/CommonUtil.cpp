#include "CommonUtil.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <algorithm>


namespace Another_CommonUtil
{
	/************************************************************************/
	/*              Normalize a scalar field to [0,1]                       */
	/************************************************************************/
	void CommonUtil::NormalizeScalar(const vector<double>& scalars, vector<double>& normalizeResult)
	{
		//assume the scalar field is not empty
		assert (scalars.size());
		
		vector<double> result;
		int cut = (int)scalars.size()/20;
		Filter(scalars,result,cut);

		//Normalize to [0,1]
		double maxium = FindMaxium(result) + 0.1;
		double minum = FindMinum(result) - 0.1;
		double len = maxium - minum;

		for (size_t i = 0; i<result.size(); i++)
		{
			double tmp = result[i];
			normalizeResult.push_back((tmp - minum)/len );
		}
		return;
	}
	void CommonUtil::AcTan(const vector<double>& scalars, vector<double>& result)
	{
		assert (scalars.size());

		for (size_t i = 0; i<scalars.size(); i++)
		{
			double tmp = atan( scalars[i] );
			result.push_back(tmp);
		}
		return;
	}


	void CommonUtil::Filter(const vector<double>&scalars, vector<double>& result, int n)
	{
		assert(scalars.size());
		vector<double> sorted;
		for (size_t i = 0; i < scalars.size(); i++)
		{
			sorted.push_back(scalars[i]);
		}

		sort(sorted.begin(), sorted.end());
		double tmpmin = sorted[n];
		double tmpmax = sorted[scalars.size()- n - 1];

		for (int i = 0; i <scalars.size(); i++)
		{
			double tmp = scalars[i];
			if (tmp > tmpmax)
			{
				result.push_back(tmpmax);
			}
			else if (tmp < tmpmin)
			{
				result.push_back((tmpmin));
			}
			else
			{
				result.push_back(tmp);
			}
		}
	}

	/************************************************************************/
	/* Find maxium entry in a scalar array                                  */
	/************************************************************************/
	double CommonUtil::FindMaxium(const vector<double>& scalars)
	{
		//assume the scalar field is not empty
		assert (scalars.size());

		//Find the maxium if not empty
		double tmpMax = 0;
		tmpMax = scalars[0];
		for (size_t i =0; i< scalars.size(); i++)
		{
			if (scalars[i]>tmpMax)
			{
				tmpMax = scalars[i];
			}
		}
		return tmpMax;
	}

	/************************************************************************/
	/* Find minum entry in a scalar array                                   */
	/************************************************************************/
	double CommonUtil::FindMinum(vector<double> scalars)
	{
		//assume the scalar field is not empty
		assert (scalars.size());

		//Find the maxium if not empty
		double tmpMin = 0;
		tmpMin = scalars[0];
		for (size_t i =0; i< scalars.size(); i++)
		{
			if (scalars[i]<tmpMin)
			{
				tmpMin = scalars[i];
			}
		}
		return tmpMin;
	}
}