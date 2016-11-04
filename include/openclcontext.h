/**
 * @brief The header with the OpenCLContext class definition.
 *
 * @file
 */

#ifndef OPENCLCONTEXT_H
#define OPENCLCONTEXT_H

#include "error.h"

#include <CL/cl_gl.h>
#include <clFFT.h>

#include <string>
#include <sstream>

/**
 * @brief Simplified error code test for OpenCL functions.
 * @param val_ The error code.
 */
#define checkClErrorCode(val_, message_) if((val_) != CL_SUCCESS) { std::stringstream ss; ss << message_; OpenCLContext::CCEC(val_, ss.str(), __FILE__, __LINE__); }

/**
 * @brief Simplified error code test for clFFT functions
 * @param val_ The error code.
 */
#define checkClfftErrorCode(val_, message_) if((val_) != CLFFT_SUCCESS) { std::stringstream ss; ss << message_; OpenCLContext::CFCEC(val_, ss.str(), __FILE__, __LINE__); }

/**
 * @brief A wrapper for cl_context.
 *
 * CL_DEVICE_TYPE_ALL is used universally in the whole class (e.g. when calling clGetDeviceIDs).
 */
class OpenCLContext
{
public:
	/**
	 * @brief OpenCLContext constructor.
	 * @param platform Used as an index to an array returned by clGetPlatformIDs().
	 * @param device Used as an index to an array returned by clGetDeviceIDs().
	 * @param shareCurrentGLContext If true, use current context during creation.
	 *
	 * The current OpenGL context is needed to setup proper communication between OpenGL and OpenCL.
	 * This is the only platform dependent code in the whole program and
	 * will probably need to be modified when the code is ported to other platforms.
	 */
	OpenCLContext(unsigned int platform, unsigned int device, bool shareCurrentGLContext = false);

	~OpenCLContext();

	/**
	 * @brief Returns the underlying OpenCL object.
	 */
	cl_context getCLContext() const
	{
		return context;
	}

	/**
	 * @brief Returns the platform id resolved during construction.
	 */
	cl_platform_id getCLPlatform() const
	{
		return platformId;
	}

	/**
	 * @brief Returns the device id resolved during construction.
	 */
	cl_device_id getCLDevice() const
	{
		return deviceId;
	}

	/**
	 * @brief Returns a human-readable string with info about the selected platform.
	 *
	 * clGetPlatformInfo() is used to retrieve this info.
	 */
	std::string getPlatformInfo() const;

	/**
	 * @brief Returns a human-readable string with info about the selected device.
	 *
	 * clGetDeviceInfo() is used to retrieve this info.
	 */
	std::string getDeviceInfo() const;

	/**
	 * @brief A convenience function for using a barrier.
	 */
	static void enqueueBarrier(cl_command_queue commandQueue, cl_event event)
	{
		cl_int err;
#if CL_1_2
		err = clEnqueueBarrierWithWaitList(commandQueue, 1, &event, nullptr);
		checkClErrorCode(err, "clEnqueueBarrierWithWaitList()");
#else
		err = clEnqueueWaitForEvents(commandQueue, 1, &event);
		checkClErrorCode(err, "clEnqueueWaitForEvents()");
#endif

		err = clReleaseEvent(event);
		checkClErrorCode(err, "clReleaseEvent()");
	}

	static void CCEC(cl_int val, std::string message, const char* file, int line);
	static void CFCEC(clfftStatus val, std::string message, const char* file, int line);

private:
	cl_context context;
	cl_platform_id platformId;
	cl_device_id deviceId;
};

#endif // OPENCLCONTEXT_H
