#ifndef OPENCLPROGRAM_H
#define OPENCLPROGRAM_H

#include "openclcontext.h"

#include <cstdio>
#include <string>
#include <stdexcept>

/**
 * @brief A wrapper for cl_program.
 */
class OpenCLProgram
{
public:
	/**
	 * @brief OpenCLProgram constructor.
	 * @param source The source string.
	 */
	OpenCLProgram(const std::string& source, OpenCLContext* context);
	~OpenCLProgram();

	/**
	 * @brief Returns a kernel object.
	 * @param kernelName The name of the kernel function.
	 *
	 * The returned kernel object is independent of this class and the caller
	 * takes its ownership.
	 */
	cl_kernel createKernel(const std::string& kernelName)
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

	/**
	 * @brief Returns cl_program compilation status.
	 * @return True if there was no error during compilation.
	 */
	bool compilationSuccessful() const
	{
		return !invalid;
	}

	/**
	 * @brief Returns a string with the compilation output (errors and warnings).
	 */
	std::string getCompilationLog() const;

private:
	cl_program program;

	bool invalid;
	OpenCLContext* context;
};

#endif // OPENCLPROGRAM_H
