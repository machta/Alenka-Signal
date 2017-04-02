#include <AlenkaSignal/montage.h>

#include <AlenkaSignal/openclcontext.h>

#include <regex>

using namespace std;

namespace
{

template<class T>
string buildSource(const string& source, const string& headerSource)
{
	// TODO: add proper indentation
	string src;

	src += "#define float ";
	src += is_same<float, T>::value ? "float" : "double";
	src += "\n";

	src += R"(
#define PARA __global float* _input_, int _inputRowLength_, int _inputRowOffset_, int _channelsInFile_
#define PASS _input_, _inputRowLength_, _inputRowOffset_, _channelsInFile_

float in(int i, PARA)
{
	return 0 <= i && i < _channelsInFile_ ? _input_[_inputRowLength_*i + _inputRowOffset_ + get_global_id(0)] : /*NAN*/0; // The NAN value makes the signal line disappear, which makes it apparent that the user made a mistake. It caused problems during compilation on some platforms, so I replaced it.
}

#define in(a_) in(a_, PASS))";

	src += headerSource;

	src += "\n\n__kernel void montage(__global float* _input_, __global float* _output_, int _inputRowLength_, int _inputRowOffset_, int _outputRowLength_, int _outputRowIndex_, int _channelsInFile_)";

	src += "\n{\n\tfloat out = 0;\n\n{\n";
	src += source;
	src += "\n}\n\n\t";

	src += "_output_[_outputRowLength_*_outputRowIndex_ + get_global_id(0)] = out;\n}\n";

	//cerr << src << endl;
	return src;
}

// This test matches the most frequently used code: "out = int(1);" and similar.
const regex re(R"(\s*out\s*=\s*in\s*\(\s*\d+\s*\)\s*;\s*)");

bool testRegex(const std::string& source)
{
	return regex_match(source, re);
}

} // namespace

namespace AlenkaSignal
{

template<class T>
Montage<T>::Montage(const string& source, OpenCLContext* context, const string& headerSource)
	: program(OpenCLProgram(buildSource<T>(source, headerSource), context))
{
	//logToFile("Constructing montage with " << source.size() << " tracks.");
}

template<class T>
Montage<T>::~Montage()
{
	cl_int err;
	err = clReleaseKernel(kernel);
	checkClErrorCode(err, "clReleaseKernel()");
}

template<class T>
bool Montage<T>::test(const std::string& source, OpenCLContext* context, string* errorMessage, const string& headerSource)
{
	//logToFile("Testing montage code.");

	if (testRegex(source))
		return true;

	// Use the OpenCL compiler to test the source.
	OpenCLProgram program(buildSource<T>(source, headerSource), context);

	if (program.compilationSuccessful())
	{
		cl_kernel kernel = program.createKernel("montage");

		cl_int err = clReleaseKernel(kernel);
		checkClErrorCode(err, "clReleaseKernel()");

		return true;
	}
	else
	{
		if (errorMessage != nullptr)
		{
			*errorMessage = "Compilation failed:\n" + program.getCompilationLog();
		}
		return false;
	}
}

template<class T>
string Montage<T>::stripComments(const string& code)
{
	const static std::regex commentre(R"((/\*([^*]|(\*+[^*/]))*\*+/)|(//.*))");
	return regex_replace(code, commentre, string(""));
}

template class Montage<float>;
template class Montage<double>;

} // namespace AlenkaSignal
