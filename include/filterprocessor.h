#ifndef FILTERPROCESSOR_H
#define FILTERPROCESSOR_H

#include "openclcontext.h"
#include "openclprogram.h"

#include <clFFT.h>

#include <cassert>
#include <vector>

/**
 * @brief This class handles filtering of data blocks.
 */
template<class T>
class FilterProcessor
{
public:
	FilterProcessor(unsigned int blockLength, unsigned int blockChannels, OpenCLContext* context);
	~FilterProcessor();

	void process(cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue);
	void changeFilter(const std::vector<T>& coefficients)
	{
		M = coefficients.size();
		this->coefficients = coefficients;

		coefficientsChanged = true;
		samplesChanged = false;
	}
	void changeSampleFilter(int M, const std::vector<T>& samples)
	{
		assert((int)samples.size() == (M + 1)/2 && "Assure the right number of samples was provided.");

		this->M = M;
		this->samples = samples;

		samplesChanged = true;
		coefficientsChanged = false;
	}
	int getDelay()
	{
		return (M - 1)/2;
	}

private:
	unsigned int blockLength;
	unsigned int blockChannels;

	int M;
	bool coefficientsChanged = false;
	std::vector<T> coefficients;
	bool samplesChanged = false;
	std::vector<T> samples;

	cl_kernel filterKernel;
	cl_kernel zeroKernel;
	cl_mem filterBuffer;

	clfftPlanHandle fftPlan;
	clfftPlanHandle fftPlanBatch;
	clfftPlanHandle ifftPlanBatch;
	clfftPlanHandle ifftPlan;
};

#endif // FILTERPROCESSOR_H
