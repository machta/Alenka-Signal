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

		if (lowpassOn && f >= lowpassF - 2/Fs*2)
		{
			val = 0;
		}
		else if (highpassOn && f <= highpassF + 1/Fs*2)
		{
			val = 0;
		}
		else if (notchOn)
		{
			double tmp = round(f/notchF);
			tmp = fabs(f - tmp*notchF);
			if (tmp <= notchWidth/M*Fs/M)
			{
				val = 0;
			}
		}

		samples[/*2**/i] = val;
		//samples[2*i + 1] = 0;
	}

	return samples;
}

template<class T>
void Filter<T>::printCoefficients(FILE* file, const vector<T>& coefficients)
{
	fprintf(file, "%lf\n%lf\n%lf\n", Fs, getLowpass(), getHighpass());
	for (unsigned int i = 0; i < M; ++i)
	{
		fprintf(file, "%lf\n", coefficients[i]);
	}
}

template class Filter<float>;
template class Filter<double>;
