#include "filterprocessor.h"

#include <cmath>
#include <complex>
#include <type_traits>

using namespace std;

namespace
{

const char* kernels =
#include "kernels.cl"
;

void printBuffer(FILE* file, float* data, int n)
{
	return;
	for (int i = 0; i < n; ++i)
	{
		fprintf(file, "%f\n", data[i]);
	}
}

void printBuffer(const std::string& filePath, cl_mem buffer, cl_command_queue queue)
{
	return;
#ifndef NDEBUG
	FILE* file = fopen(filePath.c_str(), "w");

	cl_int err;

	size_t size;
	err = clGetMemObjectInfo(buffer, CL_MEM_SIZE, sizeof(size_t), &size, nullptr);
	checkClErrorCode(err, "clGetMemObjectInfo");

	float* tmp = new float[size/sizeof(float)];

	err = clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, size, tmp, 0, nullptr, nullptr);
	checkClErrorCode(err, "clEnqueueReadBuffer");

	printBuffer(file, tmp, size/sizeof(float));

	delete[] tmp;

	fclose(file);
#endif
}

} // namespace

template<class T, bool test>
FilterProcessor<T, test>::FilterProcessor(unsigned int blockLength, unsigned int channels, OpenCLContext* context)
	: blockLength(blockLength), blockChannels(channels)
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

	filterBuffer = clCreateBuffer(context->getCLContext(), test ? CL_MEM_READ_WRITE : flags, blockLength*sizeof(T), nullptr, &err);
	checkClErrorCode(err, "clCreateBuffer");

	// Construct the fft plans.
	size_t size = blockLength;
	size_t bufferDistance = size;

	errFFT = clfftCreateDefaultPlan(&fftPlan, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	clfftSetPlanPrecision(fftPlan, precision);
	clfftSetLayout(fftPlan, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
	clfftSetResultLocation(fftPlan, CLFFT_INPLACE);
	clfftSetPlanBatchSize(fftPlan, 1);
	//clfftSetPlanDistance(fftPlan, bufferDistance, bufferDistance/2);

	errFFT = clfftCreateDefaultPlan(&fftPlanBatch, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	clfftSetPlanPrecision(fftPlanBatch, precision);
	clfftSetLayout(fftPlanBatch, CLFFT_REAL, CLFFT_HERMITIAN_INTERLEAVED);
	clfftSetResultLocation(fftPlanBatch, CLFFT_OUTOFPLACE);
	clfftSetPlanBatchSize(fftPlanBatch, blockChannels);
	clfftSetPlanDistance(fftPlanBatch, bufferDistance, bufferDistance/2);

	errFFT = clfftCreateDefaultPlan(&ifftPlanBatch, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	clfftSetPlanPrecision(ifftPlanBatch, precision);
	clfftSetLayout(ifftPlanBatch, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
	clfftSetResultLocation(ifftPlanBatch, CLFFT_INPLACE);
	clfftSetPlanBatchSize(ifftPlanBatch, blockChannels);
	clfftSetPlanDistance(ifftPlanBatch, bufferDistance/2, bufferDistance);

	errFFT = clfftCreateDefaultPlan(&ifftPlan, context->getCLContext(), CLFFT_1D, &size);
	checkClfftErrorCode(errFFT, "clfftCreateDefaultPlan()");
	clfftSetPlanPrecision(ifftPlan, precision);
	clfftSetLayout(ifftPlan, CLFFT_HERMITIAN_INTERLEAVED, CLFFT_REAL);
	clfftSetResultLocation(ifftPlan, CLFFT_INPLACE);
	clfftSetPlanBatchSize(ifftPlan, 1);
	//clfftSetPlanDistance(ifftPlan, bufferDistance, bufferDistance/2);
}

template<class T, bool test>
FilterProcessor<T, test>::~FilterProcessor()
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
	errFFT = clfftDestroyPlan(&ifftPlan);
	checkClfftErrorCode(errFFT, "clfftDestroyPlan()");
}

template<class T, bool test>
void FilterProcessor<T, test>::process(cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue)
{
	cl_int err;
	clfftStatus errFFT;

	if (coefficientsChanged)
	{
		assert(samplesChanged == false);

		err = clEnqueueWriteBuffer(queue, filterBuffer, CL_FALSE, 0, M*sizeof(T), coefficients.data(), 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueWriteBuffer()");
	}

	if (samplesChanged)
	{
		assert(coefficientsChanged == false);

		int cM = 1 + M/2;
		coefficients.insert(coefficients.begin(), 2*cM, 0);

		for (unsigned int i = 0; i < /*cM*/samples.size(); ++i)
		{
			assert(i < /*(int)*/samples.size());

			coefficients[2*i] = samples[i];
		}

		// Multiply Hr by exp(...) to make the frequency response H. (eq. 10.2.35)
		for (int i = 0; i < cM; ++i)
		{
			complex<T> tmp(0, 1);
			tmp *= -2*M_PI*i*(M - 1)/2/M;
			tmp = exp(tmp);

			complex<T> tmp2(coefficients[2*i], coefficients[2*i + 1]);
			tmp *= tmp2;

			coefficients[2*i] = tmp.real();
			coefficients[2*i + 1] = tmp.imag();
		}

		size_t size = M;
		clfftSetPlanLength(ifftPlan, CLFFT_1D, &size);
		clfftSetPlanScale(ifftPlan, CLFFT_BACKWARD, 1./size);

		err = clEnqueueWriteBuffer(queue, filterBuffer, CL_TRUE, 0, 2*cM*sizeof(T), coefficients.data(), 0, nullptr, nullptr);

		// Compute the iFFT of H to make the FIR filter coefficients h. (eq. 10.2.33)
		errFFT = clfftEnqueueTransform(ifftPlan, CLFFT_BACKWARD, 1, &queue, 0, nullptr, nullptr, &filterBuffer, nullptr, nullptr);
		checkClfftErrorCode(errFFT, "clfftEnqueueTransform()");

		if (test)
		{
			err = clEnqueueReadBuffer(queue, filterBuffer, CL_TRUE, 0, 2*cM*sizeof(T), coefficients.data(), 0, nullptr, nullptr);
			checkClErrorCode(err, "clEnqueueReadBuffer()");
		}
	}

	if (coefficientsChanged || samplesChanged)
	{
		coefficientsChanged = samplesChanged = false;

		printBuffer("before_filterBuffer.txt", filterBuffer, queue);

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

		printBuffer("after_filterBuffer_zero.txt", filterBuffer, queue);

		errFFT = clfftEnqueueTransform(fftPlan, CLFFT_FORWARD, 1, &queue, 0, nullptr, nullptr, &filterBuffer, nullptr, nullptr);
		checkClfftErrorCode(errFFT, "clfftEnqueueTransform");

		printBuffer("after_filterBuffer.txt", filterBuffer, queue);
	}

	// TODO: apply a window function (on the device)

	printBuffer("before_fft.txt", inBuffer, queue);

	// FFT.
	errFFT = clfftEnqueueTransform(fftPlanBatch, CLFFT_FORWARD, 1, &queue, 0, nullptr, nullptr, &inBuffer, &outBuffer, nullptr);
	checkClfftErrorCode(errFFT, "clfftEnqueueTransform");

	printBuffer("after_fft.txt", outBuffer, queue);

	// Multiply.
	err = clSetKernelArg(filterKernel, 0, sizeof(cl_mem), &outBuffer);
	checkClErrorCode(err, "clSetKernelArg()");

	err = clSetKernelArg(filterKernel, 1, sizeof(cl_mem), &filterBuffer);
	checkClErrorCode(err, "clSetKernelArg()");

	size_t globalWorkSize[2] = {blockChannels, blockLength/2};

	err = clEnqueueNDRangeKernel(queue, filterKernel, 2, nullptr, globalWorkSize, nullptr, 0, nullptr, nullptr);
	checkClErrorCode(err, "clEnqueueNDRangeKernel()");

	printBuffer("after_multiply.txt", outBuffer, queue);

	// IFFT.
	errFFT = clfftEnqueueTransform(ifftPlanBatch, CLFFT_BACKWARD, 1, &queue, 0, nullptr, nullptr, &outBuffer, nullptr, nullptr);
	checkClfftErrorCode(errFFT, "clfftEnqueueTransform");

	printBuffer("after_ifft.txt", outBuffer, queue);
}

template class FilterProcessor<float, false>;
template class FilterProcessor<float, true>;
template class FilterProcessor<double, false>;
template class FilterProcessor<double, true>;
