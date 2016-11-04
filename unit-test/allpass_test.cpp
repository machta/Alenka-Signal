#include "gtest/gtest.h"

#include "filterprocessor.h"

#include <fstream>
#include <string>
#include <memory>

using namespace std;

namespace
{


} // namespace

TEST(allpass_test, sample)
{
/*	clfftStatus errFFT;
	clfftSetupData setupData;

	errFFT = clfftInitSetupData(&setupData);
	//checkClfftErrorCode(errFFT, "clfftInitSetupData()");

	errFFT = clfftSetup(&setupData);
	//checkClfftErrorCode(errFFT, "clfftSetup()");*/


	int n = 20;

	OpenCLContext context(0, 0);
	FilterProcessor<float> processor(n, 1, &context);

	vector<float> signal;
	vector<float> output(n);
	for (int i = 1; i <= n; i++)
		signal.push_back(i);

	cl_command_queue queue = clCreateCommandQueue(context.getCLContext(), context.getCLDevice(), 0, nullptr);

	cl_mem_flags flags = CL_MEM_READ_WRITE;
	cl_mem inBuffer = clCreateBuffer(context.getCLContext(), flags | CL_MEM_COPY_HOST_PTR, n*sizeof(float), signal.data(), nullptr);
	cl_mem outBuffer = clCreateBuffer(context.getCLContext(), flags, n*sizeof(float), nullptr, nullptr);

	processor.changeSampleFilter(5, vector<float>{1, 1, 1});

	processor.process(inBuffer, outBuffer, queue);

	clEnqueueReadBuffer(queue, outBuffer, CL_TRUE, 0, n*sizeof(float), output.data(), 0, nullptr, nullptr);


//	/*clfftStatus*/ errFFT = clfftTeardown();
	//checkClfftErrorCode(errFFT, "clfftTeardown()");
}

TEST(allpass_test, coefficient)
{
	int n = 20;

	OpenCLContext context(0, 0);
	FilterProcessor<float> processor(n, 1, &context);

	vector<float> signal;
	vector<float> output(n);
	for (int i = 1; i <= n; i++)
		signal.push_back(i);

	cl_command_queue queue = clCreateCommandQueue(context.getCLContext(), context.getCLDevice(), 0, nullptr);

	cl_mem_flags flags = CL_MEM_READ_WRITE;
	cl_mem inBuffer = clCreateBuffer(context.getCLContext(), flags | CL_MEM_COPY_HOST_PTR, n*sizeof(float), signal.data(), nullptr);
	cl_mem outBuffer = clCreateBuffer(context.getCLContext(), flags, n*sizeof(float), nullptr, nullptr);

	processor.changeFilter(vector<float>{0, 0, 1, 0, 0});

	processor.process(inBuffer, outBuffer, queue);

	clEnqueueReadBuffer(queue, outBuffer, CL_TRUE, 0, n*sizeof(float), output.data(), 0, nullptr, nullptr);
}
