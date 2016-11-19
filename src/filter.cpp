#include "filter.h"

#include <complex>

using namespace std;

template<class T>
vector<T> Filter<T>::computeSamples()
{
	vector<T> samples((M + 1)/2);

	//int cM = 1 + M/2;

	// Initialize samples with the values of Hr.
	for (unsigned int i = 0; i < samples.size(); ++i)
	{
		double f = 2.*i/M;
		double val = 1;

		if (lowpassOn && f >= lowpassF)
		{
			val = 0;
		}
		else if (highpassOn && f <= highpassF)
		{
			val = 0;
		}
		else if (notchOn)
		{
			double tmp = round(f/notchF);
			tmp = fabs(f - tmp*notchF);
			if (tmp <= 3./M*Fs/M) // Possibly turn the '3.' into a parameter.
			{
				val = 0;
			}
		}

		samples[/*2**/i] = val;
		//samples[2*i + 1] = 0;
	}

	return samples;
}

template class Filter<float>;
template class Filter<double>;
