#ifndef ALENKASIGNAL_FILTERPROCESSOR_H
#define ALENKASIGNAL_FILTERPROCESSOR_H

#include <clFFT.h>

#include <cassert>
#include <vector>


namespace AlenkaSignal
{

class OpenCLContext;

enum class WindowFunction
{
	None, Hamming, Blackman
};

/**
 * @brief This class handles filtering of data blocks.
 */
template<class T>
class FilterProcessor
{
public:
	FilterProcessor(unsigned int blockLength, unsigned int blockChannels, OpenCLContext* context, WindowFunction windowFunction = WindowFunction::None);
	~FilterProcessor();

	void process(cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue);
	void changeFilter(const std::vector<T>& coefficients)
	{
		coefficientsChanged = true;
		M = coefficients.size();
		this->coefficients = coefficients;
	}
	void changeSampleFilter(int M, const std::vector<T>& samples);
	int delaySamples() const
	{
		return (M - 1)/2;
	}
	int discardSamples() const
	{
		return M - 1;
	}
	std::vector<T> getCoefficients() const
	{
		return coefficients;
	}

private:
	unsigned int blockLength;
	unsigned int blockChannels;
	WindowFunction windowFunction;

	int M;
	bool coefficientsChanged = false;
	std::vector<T> coefficients;

	cl_kernel filterKernel;
	cl_kernel zeroKernel;
	cl_mem filterBuffer;

	clfftPlanHandle fftPlan;
	clfftPlanHandle fftPlanBatch;
	clfftPlanHandle ifftPlanBatch;
};

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_FILTERPROCESSOR_H
