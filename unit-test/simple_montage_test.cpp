#include "gtest/gtest.h"

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/montage.h>
#include <AlenkaSignal/montageprocessor.h>

#include <functional>

using namespace std;
using namespace AlenkaSignal;

namespace
{

void compareFloat(float a, float b)
{
	EXPECT_FLOAT_EQ(a, b);
}

void compareDouble(double a, double b)
{
	EXPECT_DOUBLE_EQ(a, b);
}

template<class T>
void test(function<void(T, T)> compare)
{
	int n = 20;
	int inChannels = 3;
	int offset = 5;

	OpenCLContext context(OPENCL_PLATFORM, OPENCL_DEVICE);
	cl_command_queue queue = clCreateCommandQueue(context.getCLContext(), context.getCLDevice(), 0, nullptr);
	MontageProcessor<T> processor(offset, n - offset, inChannels);

	string src = "out = in(0);";
	string msg;

	bool res = Montage<T>::test(src, &context, &msg);
	cerr << msg << endl;
	ASSERT_TRUE(res);
	Montage<T> m1(src, &context);

	Montage<T> m2("out = in(1);", &context);
	Montage<T> m3("out = in(0) + in(1);", &context);
	Montage<T> m4("out = in(2)*3.14;", &context);
	Montage<T> m5("out = -1;", &context);
	vector<Montage<T>*> montage = {&m1, &m2, &m3, &m4, &m5};

	vector<T> signal;
	for (int j = 0; j < inChannels; j++)
	for (int i = 1; i <= n; i++)
	signal.push_back(10*pow(10,j) + i);

	vector<T> output((n - offset)*montage.size());

	cl_mem_flags flags = CL_MEM_READ_WRITE;
	cl_mem inBuffer = clCreateBuffer(context.getCLContext(), flags | CL_MEM_COPY_HOST_PTR, n*inChannels*sizeof(T), signal.data(), nullptr);
	cl_mem outBuffer = clCreateBuffer(context.getCLContext(), flags, (n - offset)*montage.size()*sizeof(T), nullptr, nullptr);

	processor.process(montage, inBuffer, outBuffer, queue);

	clEnqueueReadBuffer(queue, outBuffer, CL_TRUE, 0, (n - offset)*montage.size()*sizeof(T), output.data(), 0, nullptr, nullptr);

	for (int i = 0; i < n - offset; i++)
	{
		compare(output[(n - offset)*0 + i], signal[offset + i]);
	}

	for (int i = 0; i < n - offset; i++)
	{
		compare(output[(n - offset)*1 + i], signal[n + offset + i]);
	}

	for (int i = 0; i < n - offset; i++)
	{
		compare(output[(n - offset)*2 + i], signal[offset + i] + signal[n + offset + i]);
	}

	for (int i = 0; i < n - offset; i++)
	{
		compare(output[(n - offset)*3 + i], signal[2*n + offset + i]*3.14);
	}

	for (int i = 0; i < n - offset; i++)
	{
		compare(output[(n - offset)*4 + i], -1);
	}
}

} // namespace

TEST(simple_montage_test, float)
{
	test<float>(&compareFloat);
}

TEST(simple_montage_test, double)
{
	test<double>(&compareDouble);
}
