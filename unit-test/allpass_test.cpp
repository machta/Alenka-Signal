#include "gtest/gtest.h"

#include "filterprocessor.h"

#include <functional>

using namespace std;

namespace
{

void compareFloat(float a, float b)
{
	//EXPECT_FLOAT_EQ(a, b);
	EXPECT_NEAR(a, b, 0.00001);
}

void compareDouble(double a, double b)
{
	//EXPECT_DOUBLE_EQ(a, b);
	EXPECT_NEAR(a, b, 0.000001);
}

template<class T>
void test(function<void(FilterProcessor<T>*)> change, function<void(T, T)> compare)
{
	clfftStatus errFFT;
	clfftSetupData setupData;

	errFFT = clfftInitSetupData(&setupData);
	checkClfftErrorCode(errFFT, "clfftInitSetupData()");

	errFFT = clfftSetup(&setupData);
	checkClfftErrorCode(errFFT, "clfftSetup()");

	{
		cl_int err;
		int n = 20;

		OpenCLContext context(OPENCL_PLATFORM, OPENCL_DEVICE);
		FilterProcessor<T> processor(n, 1, &context);

		change(&processor);

		vector<T> signal(processor.discardSamples() - processor.delaySamples());
		vector<T> output(n);
		for (int i = 1; i <= n - processor.discardSamples() + processor.delaySamples(); i++)
			signal.push_back(i);

		cl_command_queue queue = clCreateCommandQueue(context.getCLContext(), context.getCLDevice(), 0, &err);
		checkClErrorCode(err, "clCreateCommandQueue");

		cl_mem_flags flags = CL_MEM_READ_WRITE;

		cl_mem inBuffer = clCreateBuffer(context.getCLContext(), flags | CL_MEM_COPY_HOST_PTR, (n + 4)*sizeof(T), signal.data(), &err);
		checkClErrorCode(err, "clCreateBuffer");

		cl_mem outBuffer = clCreateBuffer(context.getCLContext(), flags, (n + 4)*sizeof(T), nullptr, &err);
		checkClErrorCode(err, "clCreateBuffer");

		processor.process(inBuffer, outBuffer, queue);

		err = clEnqueueReadBuffer(queue, outBuffer, CL_TRUE, 0, n*sizeof(T), output.data(), 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueReadBuffer");

		for (int i = 0; i < n - processor.discardSamples(); i++)
			compare(output[i + processor.discardSamples()], signal[i + processor.discardSamples() - processor.delaySamples()]);

		err = clReleaseCommandQueue(queue);
		checkClErrorCode(err, "clReleaseCommandQueue");

		err = clReleaseMemObject(inBuffer);
		checkClErrorCode(err, "clReleaseMemObject");

		err = clReleaseMemObject(outBuffer);
		checkClErrorCode(err, "clReleaseMemObject");
	}

	errFFT = clfftTeardown();
	checkClfftErrorCode(errFFT, "clfftTeardown()");
}

} // namespace

TEST(allpass_test, sample_float)
{
	auto f = [] (FilterProcessor<float>* p) { p->changeSampleFilter(5, vector<float>{1, 1, 1}); };
	test<float>(f, &compareFloat);
}

TEST(allpass_test, coefficient_float)
{
	auto f = [] (FilterProcessor<float>* p) { p->changeFilter(vector<float>{0, 0, 1, 0, 0}); };
	test<float>(f, &compareFloat);
}

TEST(allpass_test, sample_double)
{
	auto f = [] (FilterProcessor<double>* p) { p->changeSampleFilter(5, vector<double>{1, 1, 1}); };
	test<double>(f, &compareDouble);
}

TEST(allpass_test, coefficient_double)
{
	auto f = [] (FilterProcessor<double>* p) { p->changeFilter(vector<double>{0, 0, 1, 0, 0}); };
	test<double>(f, &compareDouble);
}
