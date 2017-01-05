#include "gtest/gtest.h"

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/montage.h>

using namespace std;
using namespace AlenkaSignal;

namespace
{

template<class T>
void test()
{
	OpenCLContext context(OPENCL_PLATFORM, OPENCL_DEVICE);

	ASSERT_TRUE(Montage<T>::test("", &context));
	ASSERT_TRUE(Montage<T>::test("out=in(1);", &context));
	ASSERT_TRUE(Montage<T>::test("  out = in(    1     )  ;    ", &context));
	ASSERT_TRUE(Montage<T>::test("out=in(00010);", &context));
	ASSERT_TRUE(Montage<T>::test("int i = 5; out=in(i);", &context));
	ASSERT_TRUE(Montage<T>::test("float sample = in(55); out = 2*sample*2;", &context));
	ASSERT_TRUE(Montage<T>::test("float sample = in(22); out = 2.2*sample*2.2;", &context));
	ASSERT_TRUE(Montage<T>::test("float sample = in(00); out = 2.2*sample*2;", &context));
	ASSERT_TRUE(Montage<T>::test("float sample = in(-1); out = 2*sample*2.2;", &context));
	ASSERT_TRUE(Montage<T>::test("out=in(10 - 2*3);", &context));
	ASSERT_TRUE(Montage<T>::test("in(10 - 2*3);", &context));

	ASSERT_FALSE(Montage<T>::test("fdfsdfssf", &context));
	ASSERT_FALSE(Montage<T>::test("out=in(1)", &context));
	ASSERT_FALSE(Montage<T>::test("out=in();", &context));
	ASSERT_FALSE(Montage<T>::test("out=in(1, 1)", &context));
	ASSERT_FALSE(Montage<T>::test("outt=in(1);", &context));
	ASSERT_FALSE(Montage<T>::test("out=IN(1);", &context));
	ASSERT_FALSE(Montage<T>::test("in(10 3);", &context));
}

} // namespace

TEST(montage_test_test, float)
{
	test<float>();
}

TEST(montage_test_test, double)
{
	test<double>();
}
