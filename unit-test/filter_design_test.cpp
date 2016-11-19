#include "gtest/gtest.h"

#include "filterprocessor.h"
#include "filter.h"

#include <functional>

using namespace std;

namespace
{

double answerLowpassD[8] = {-0.105969883127822, 0.0293291419087276, 0.220670858091272, 0.355969883127822, 0.355969883127822, 0.220670858091272, 0.0293291419087275, -0.105969883127822};
float answerLowpassF[8]  = {-0.105969883127822, 0.0293291419087276, 0.220670858091272, 0.355969883127822, 0.355969883127822, 0.220670858091272, 0.0293291419087275, -0.105969883127822};

void compareFloat(float a, float b)
{
	//EXPECT_FLOAT_EQ(a, b);
	EXPECT_NEAR(a, b, 0.0000001);
}

void compareDouble(double a, double b)
{
	//EXPECT_DOUBLE_EQ(a, b);
	EXPECT_NEAR(a, b, 0.000000001);
}

template<class T>
void test(function<void(T, T)> compare, T* answer)
{
	int n = 20;

	OpenCLContext context(0, 0);
	FilterProcessor<T, true> processor(n, 1, &context);

	vector<T> signal;
	vector<T> output(n);
	for (int i = 1; i <= n; i++)
	signal.push_back(i);

	cl_command_queue queue = clCreateCommandQueue(context.getCLContext(), context.getCLDevice(), 0, nullptr);

	cl_mem_flags flags = CL_MEM_READ_WRITE;
	cl_mem inBuffer = clCreateBuffer(context.getCLContext(), flags | CL_MEM_COPY_HOST_PTR, n*sizeof(T), signal.data(), nullptr);
	cl_mem outBuffer = clCreateBuffer(context.getCLContext(), flags, n*sizeof(T), nullptr, nullptr);

	Filter<T> filter(8, 200);
	filter.lowpass(true);
	filter.setLowpass(50);
	processor.changeSampleFilter(8, filter.computeSamples());

	processor.process(inBuffer, outBuffer, queue);

	clEnqueueReadBuffer(queue, outBuffer, CL_TRUE, 0, n*sizeof(T), output.data(), 0, nullptr, nullptr);

	auto res = processor.getCoefficients();
	for (int i = 0; i < 8; ++i)
	{
		compare(res[i], answer[i]);
	}
}

} // namespace

TEST(filter_design_test, simple_float)
{
	test<float>(&compareFloat, answerLowpassF);
}

TEST(filter_design_test, simple_double)
{
	test<double>(&compareDouble, answerLowpassD);
}
