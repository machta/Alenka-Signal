#ifndef ALENKASIGNAL_OPENCLPROGRAM_H
#define ALENKASIGNAL_OPENCLPROGRAM_H

#include <CL/cl_gl.h>

#include <cstdio>
#include <string>
#include <stdexcept>

namespace AlenkaSignal
{

class OpenCLContext;

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
	cl_kernel createKernel(const std::string& kernelName);

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

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_OPENCLPROGRAM_H