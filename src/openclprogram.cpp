#include <AlenkaSignal/openclprogram.h>

#include <AlenkaSignal/openclcontext.h>

#include <cstdlib>
#include <cassert>

using namespace std;

namespace AlenkaSignal
{

OpenCLProgram::OpenCLProgram(const string& source, OpenCLContext* context) : context(context)
{
	cl_int err;

	const char* sourcePointer = source.c_str();
	size_t size = source.size();

	program = clCreateProgramWithSource(context->getCLContext(), 1, &sourcePointer, &size, &err);
	checkClErrorCode(err, "clCreateProgramWithSource()");

	err = clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);

	if (err == CL_SUCCESS)
	{
		invalid = false;
	}
	else
	{
		cl_build_status status;

		cl_int err2 = clGetProgramBuildInfo(program, context->getCLDevice(), CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, nullptr);
		checkClErrorCode(err2, "clGetProgramBuildInfo()");

		assert(status != CL_BUILD_IN_PROGRESS);

		invalid = status == CL_BUILD_ERROR;

		if (invalid)
		{
			string log = getCompilationLog();
			//logToFileAndConsole(log);
		}
		else
		{
			checkClErrorCode(err, "clBuildProgram()");
		}
	}
}

OpenCLProgram::~OpenCLProgram()
{
	cl_int err = clReleaseProgram(program);
	checkClErrorCode(err, "clReleaseProgram()");
}

cl_kernel OpenCLProgram::createKernel(const string& kernelName)
{
	if (compilationSuccessful() == false)
	{
		throw std::runtime_error("Cannot create kernel object from an OpenCLProgram that failed to compile.");
	}

	cl_int err;

	cl_kernel kernel = clCreateKernel(program, kernelName.c_str(), &err);
	checkClErrorCode(err, "clCreateKernel()");

	return kernel;
}

string OpenCLProgram::getCompilationLog() const
{
	size_t logLength;

	cl_int err = clGetProgramBuildInfo(program, context->getCLDevice(), CL_PROGRAM_BUILD_LOG, 0, nullptr, &logLength);
	checkClErrorCode(err, "clGetProgramBuildInfo()");

	char* tmp = new char[logLength + 1];
	tmp[logLength] = 0;

	err = clGetProgramBuildInfo(program, context->getCLDevice(), CL_PROGRAM_BUILD_LOG, logLength, tmp, nullptr);
	checkClErrorCode(err, "clGetProgramBuildInfo()");

	string str(tmp);

	delete[] tmp;

	return str;
}

} // namespace AlenkaSignal
