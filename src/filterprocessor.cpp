#include <AlenkaSignal/filterprocessor.h>

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/openclprogram.h>
#include <fasttransforms.h>

#include <cmath>
#include <complex>
#include <type_traits>

using namespace std;

namespace
{

const char* kernels =
#include "kernels.cl"
;

template<class T>
T hammingWindow(int n, int M)
{
	const T tmp = 2*M_PI*n/(M - 1);
	return 0.54 - 0.46*cos(tmp);
}

template<class T>
T blackmanWindow(int n, int M)
{
	const T a = 0.16, a0 = (1 - a)/2, a1 = 0.5, a2 = a/2, tmp = 2*M_PI*n/(M - 1);
	return a0 - a1*cos(tmp) + a2*cos(2*tmp);
}

} // namespace

namespace AlenkaSignal
{

template<class T>
FilterProcessor<T>::FilterProcessor(unsigned int blockLength, unsigned int channels, OpenCLContext* context, WindowFunction windowFunction)
	: blockLength(blockLength), blockChannels(channels), windowFunction(windowFunction)
{
	assert(blockLength%2 == 0);

	cl_int err;
	clfftStatus errFFT;

	clfftPrecision precision = CLFFT_SINGLE;
	string kernelsSource;

	if (is_same<double, T>::value)
	{
		precision = CLFFT_DOUBLE;
		kernelsSource = "#define float double\n#define float2 double2\n\n";
	}

	kernelsSource += kernels;
	OpenCLProgram program(kernelsSource, context);

	filterKernel = program.createKernel("filter");
	zeroKernel = program.createKernel("zero");

	cl_mem_flags flags = CL_MEM_READ_WRITE;
#ifdef NDEBUG
#if CL_1_2
	flags |= CL_MEM_HOST_WRITE_ONLY;
#endif
#endif

	filterBuffer = clCreateBuffer(context->getCLContext(), flags, (blockLength + 4)*sizeof(T), nullptr, &err);
	checkClErrorCode(err, "clCreateBuffer");

	// Construct the fft plans.
	size_t size = blockLength;
	size_t bufferDistance = size + 4;

	errFFT = clfftCreateDefaultPlan(&fftPlan, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	errFFT = clfftSetPlanPrecision(fftPlan, precision);
	checkClfftErrorCode(errFFT, "clfftSetPlanPrecision()");
	errFFT = clfftSetLayout(fftPlan, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
	checkClfftErrorCode(errFFT, "clfftSetLayout()");
	errFFT = clfftSetResultLocation(fftPlan, CLFFT_INPLACE);
	checkClfftErrorCode(errFFT, "clfftSetResultLocation()");
	errFFT = clfftSetPlanBatchSize(fftPlan, 1);
	checkClfftErrorCode(errFFT, "clfftSetPlanBatchSize()");
	//clfftSetPlanDistance(fftPlan, bufferDistance, bufferDistance/2);

	errFFT = clfftCreateDefaultPlan(&fftPlanBatch, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	errFFT = clfftSetPlanPrecision(fftPlanBatch, precision);
	checkClfftErrorCode(errFFT, "clfftSetPlanPrecision()");
	errFFT = clfftSetLayout(fftPlanBatch, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
	checkClfftErrorCode(errFFT, "clfftSetLayout()");
	errFFT = clfftSetResultLocation(fftPlanBatch, CLFFT_OUTOFPLACE);
	checkClfftErrorCode(errFFT, "clfftSetResultLocation()");
	errFFT = clfftSetPlanBatchSize(fftPlanBatch, blockChannels);
	checkClfftErrorCode(errFFT, "clfftSetPlanBatchSize()");
	errFFT = clfftSetPlanDistance(fftPlanBatch, bufferDistance, bufferDistance/2);
	checkClfftErrorCode(errFFT, "clfftSetPlanDistance()");

	errFFT = clfftCreateDefaultPlan(&ifftPlanBatch, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	errFFT = clfftSetPlanPrecision(ifftPlanBatch, precision);
	checkClfftErrorCode(errFFT, "clfftSetPlanPrecision()");
	errFFT = clfftSetLayout(ifftPlanBatch, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
	checkClfftErrorCode(errFFT, "clfftSetLayout()");
	errFFT = clfftSetResultLocation(ifftPlanBatch, CLFFT_INPLACE);
	checkClfftErrorCode(errFFT, "clfftSetResultLocation()");
	errFFT = clfftSetPlanBatchSize(ifftPlanBatch, blockChannels);
	checkClfftErrorCode(errFFT, "clfftSetPlanBatchSize()");
	errFFT = clfftSetPlanDistance(ifftPlanBatch, bufferDistance/2, bufferDistance);
	checkClfftErrorCode(errFFT, "clfftSetPlanDistance()");
}

template<class T>
FilterProcessor<T>::~FilterProcessor()
{
	cl_int err;
	err = clReleaseKernel(filterKernel);
	checkClErrorCode(err, "clReleaseKernel()");
	err = clReleaseKernel(zeroKernel);
	checkClErrorCode(err, "clReleaseKernel()");
	err = clReleaseMemObject(filterBuffer);
	checkClErrorCode(err, "clReleaseMemObject()");

	clfftStatus errFFT;
	errFFT = clfftDestroyPlan(&fftPlan);
	checkClfftErrorCode(errFFT, "clfftDestroyPlan()");
	errFFT = clfftDestroyPlan(&fftPlanBatch);
	checkClfftErrorCode(errFFT, "clfftDestroyPlan()");
	errFFT = clfftDestroyPlan(&ifftPlanBatch);
	checkClfftErrorCode(errFFT, "clfftDestroyPlan()");
}

template<class T>
void FilterProcessor<T>::process(cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue)
{
	cl_int err;
	clfftStatus errFFT;

#ifndef NDEBUG
	{
		size_t inSize;
		err = clGetMemObjectInfo(inBuffer, CL_MEM_SIZE, sizeof(size_t), &inSize, nullptr);
		checkClErrorCode(err, "clGetMemObjectInfo");

		size_t outSize;
		err = clGetMemObjectInfo(outBuffer, CL_MEM_SIZE, sizeof(size_t), &outSize, nullptr);
		checkClErrorCode(err, "clGetMemObjectInfo");

		assert(inSize >= (blockLength + 4)*blockChannels*sizeof(T) && "The inBuffer is too small.");
		assert(outSize >= (blockLength + 4)*blockChannels*sizeof(T) && "The inBuffer is too small.");
	}
#endif

	if (coefficientsChanged)
	{
		err = clEnqueueWriteBuffer(queue, filterBuffer, CL_TRUE, 0, M*sizeof(T), coefficients.data(), 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueWriteBuffer()");

		//printBuffer("before_filterBuffer.txt", filterBuffer, queue);

		// This section is disabled because of a bug in the implementation of clEnqueueFillBuffer().
//#if CL_1_2
//		float zero = 0;
//		err = clEnqueueFillBuffer(queue, filterBuffer, &zero, sizeof(zero), 0, width + 4, 0, nullptr, nullptr);
//		checkClErrorCode(err, "clEnqueueFillBuffer()");
//#else
		err = clSetKernelArg(zeroKernel, 0, sizeof(cl_mem), &filterBuffer);
		checkClErrorCode(err, "clSetKernelArg()");

		size_t globalWorkSize = blockLength;
		size_t globalWorkOffset = M;
		err = clEnqueueNDRangeKernel(queue, zeroKernel, 1, &globalWorkOffset, &globalWorkSize, nullptr, 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueNDRangeKernel()");
//#endif

		//printBuffer("after_filterBuffer_zero.txt", filterBuffer, queue);

		errFFT = clfftEnqueueTransform(fftPlan, CLFFT_FORWARD, 1, &queue, 0, nullptr, nullptr, &filterBuffer, nullptr, nullptr);
		checkClfftErrorCode(errFFT, "clfftEnqueueTransform");

		//printBuffer("after_filterBuffer.txt", filterBuffer, queue);
	}

	// TODO: apply a window function (on the device)

	//OpenCLContext::printBuffer("before_fft.txt", inBuffer, queue);

	// FFT.
	errFFT = clfftEnqueueTransform(fftPlanBatch, CLFFT_FORWARD, 1, &queue, 0, nullptr, nullptr, &inBuffer, &outBuffer, nullptr);
	checkClfftErrorCode(errFFT, "clfftEnqueueTransform");

	//OpenCLContext::printBuffer("after_fft.txt", outBuffer, queue);

	// Multiply.
	err = clSetKernelArg(filterKernel, 0, sizeof(cl_mem), &outBuffer);
	checkClErrorCode(err, "clSetKernelArg()");

	err = clSetKernelArg(filterKernel, 1, sizeof(cl_mem), &filterBuffer);
	checkClErrorCode(err, "clSetKernelArg()");

	size_t globalWorkSize[2] = {blockChannels, blockLength/2};

	err = clEnqueueNDRangeKernel(queue, filterKernel, 2, nullptr, globalWorkSize, nullptr, 0, nullptr, nullptr);
	checkClErrorCode(err, "clEnqueueNDRangeKernel()");

	//OpenCLContext::printBuffer("after_multiply.txt", outBuffer, queue);

	// IFFT.
	errFFT = clfftEnqueueTransform(ifftPlanBatch, CLFFT_BACKWARD, 1, &queue, 0, nullptr, nullptr, &outBuffer, nullptr, nullptr);
	checkClfftErrorCode(errFFT, "clfftEnqueueTransform");

	//OpenCLContext::printBuffer("after_ifft.txt", outBuffer, queue);
}

template<class T>
void FilterProcessor<T>::changeSampleFilter(int M, const std::vector<T>& samples)
{
	assert((int)samples.size() == (M + 1)/2 && "Assure the right number of samples was provided.");

	coefficientsChanged = true;
	this->M = M;
	//this->samples = samples;

	int cM = 1 + M/2;

	alglib::complex_1d_array inArray;
	inArray.setlength(cM);
	inArray[cM - 1].x = 0;
	inArray[cM - 1].y = 0;

	for (unsigned int i = 0; i < /*cM*/samples.size(); ++i)
	{
		assert(i < /*(int)*/samples.size());

		inArray[i].x = samples[i];
		inArray[i].y = 0;
	}

	// Multiply Hr by exp(...) to make the frequency response H. (eq. 10.2.35)
	for (int i = 0; i < cM; ++i)
	{
		/*complex<T> tmp(0, 1);
		tmp *= -2*M_PI*i*(M - 1)/2/M;
		tmp = exp(tmp);

		complex<T> tmp2(coefficients[2*i], coefficients[2*i + 1]);
		tmp *= tmp2;

		coefficients[2*i] = tmp.real();
		coefficients[2*i + 1] = tmp.imag();*/

		complex<double> tmp(0, 1);
		tmp *= -2*M_PI*i*(M - 1)/2/M;
		tmp = exp(tmp);

		inArray[i] *= alglib::complex(tmp.real(), tmp.imag());
	}

	// Compute the iFFT of H to make the FIR filter coefficients h. (eq. 10.2.33)
	alglib::real_1d_array outArray;
	outArray.setlength(M);
	alglib::fftr1dinv(inArray, M, outArray);

	coefficients.resize(M);
	for (int i = 0; i < M; i++)
		coefficients[i] = outArray[i];

	// Try to improve filter characteristics by applying a window function.
	if (windowFunction == WindowFunction::Hamming)
	{
		for (int i = 0; i < M; ++i)
		{
			coefficients[i] *= hammingWindow<T>(i, M);
		}
	}
	else if (windowFunction == WindowFunction::Blackman)
	{
		for (int i = 0; i < M; ++i)
		{
			coefficients[i] *= blackmanWindow<T>(i, M);
		}
	}
}

template class FilterProcessor<float>;
template class FilterProcessor<double>;

} // namespace AlenkaSignal
