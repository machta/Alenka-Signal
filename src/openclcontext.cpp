#include <AlenkaSignal/openclcontext.h>

#include <algorithm>
#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>

#if defined WIN_BUILD
#include <windows.h>
#elif defined UNIX_BUILD
#include <GL/glx.h>
#endif

using namespace std;

namespace
{
template <typename T>
string errorCodeToString(T val)
{
	using namespace std;

	stringstream ss;

	ss << dec << val << "(0x" << hex << val << dec << ")";

	return ss.str();
}

#define CASE(a_) case a_: return #a_

string clErrorCodeToString(cl_int code)
{
	using namespace std;

	switch (code)
	{
		CASE(CL_SUCCESS);
		CASE(CL_DEVICE_NOT_FOUND);
		CASE(CL_DEVICE_NOT_AVAILABLE);
		CASE(CL_COMPILER_NOT_AVAILABLE);
		CASE(CL_MEM_OBJECT_ALLOCATION_FAILURE);
		CASE(CL_OUT_OF_RESOURCES);
		CASE(CL_OUT_OF_HOST_MEMORY);
		CASE(CL_PROFILING_INFO_NOT_AVAILABLE);
		CASE(CL_MEM_COPY_OVERLAP);
		CASE(CL_IMAGE_FORMAT_MISMATCH);
		CASE(CL_IMAGE_FORMAT_NOT_SUPPORTED);
		CASE(CL_BUILD_PROGRAM_FAILURE);
		CASE(CL_MAP_FAILURE);
		CASE(CL_MISALIGNED_SUB_BUFFER_OFFSET);
		CASE(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST);
#if CL_1_2
		CASE(CL_COMPILE_PROGRAM_FAILURE);
		CASE(CL_LINKER_NOT_AVAILABLE);
		CASE(CL_LINK_PROGRAM_FAILURE);
		CASE(CL_DEVICE_PARTITION_FAILED);
		CASE(CL_KERNEL_ARG_INFO_NOT_AVAILABLE);
#endif
		CASE(CL_INVALID_VALUE);
		CASE(CL_INVALID_DEVICE_TYPE);
		CASE(CL_INVALID_PLATFORM);
		CASE(CL_INVALID_DEVICE);
		CASE(CL_INVALID_CONTEXT);
		CASE(CL_INVALID_QUEUE_PROPERTIES);
		CASE(CL_INVALID_COMMAND_QUEUE);
		CASE(CL_INVALID_HOST_PTR);
		CASE(CL_INVALID_MEM_OBJECT);
		CASE(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);
		CASE(CL_INVALID_IMAGE_SIZE);
		CASE(CL_INVALID_SAMPLER);
		CASE(CL_INVALID_BINARY);
		CASE(CL_INVALID_BUILD_OPTIONS);
		CASE(CL_INVALID_PROGRAM);
		CASE(CL_INVALID_PROGRAM_EXECUTABLE);
		CASE(CL_INVALID_KERNEL_NAME);
		CASE(CL_INVALID_KERNEL_DEFINITION);
		CASE(CL_INVALID_KERNEL);
		CASE(CL_INVALID_ARG_INDEX);
		CASE(CL_INVALID_ARG_VALUE);
		CASE(CL_INVALID_ARG_SIZE);
		CASE(CL_INVALID_KERNEL_ARGS);
		CASE(CL_INVALID_WORK_DIMENSION);
		CASE(CL_INVALID_WORK_GROUP_SIZE);
		CASE(CL_INVALID_WORK_ITEM_SIZE);
		CASE(CL_INVALID_GLOBAL_OFFSET);
		CASE(CL_INVALID_EVENT_WAIT_LIST);
		CASE(CL_INVALID_EVENT);
		CASE(CL_INVALID_OPERATION);
		CASE(CL_INVALID_GL_OBJECT);
		CASE(CL_INVALID_BUFFER_SIZE);
		CASE(CL_INVALID_MIP_LEVEL);
		CASE(CL_INVALID_GLOBAL_WORK_SIZE);
		CASE(CL_INVALID_PROPERTY);
#if CL_1_2
		CASE(CL_INVALID_IMAGE_DESCRIPTOR);
		CASE(CL_INVALID_COMPILER_OPTIONS);
		CASE(CL_INVALID_LINKER_OPTIONS);
		CASE(CL_INVALID_DEVICE_PARTITION_COUNT);
#endif
	}

	return "unknown code " + errorCodeToString(code);
}

string clfftErrorCodeToString(clfftStatus code)
{
	switch (code)
	{
		CASE(CLFFT_BUGCHECK);
		CASE(CLFFT_NOTIMPLEMENTED);
		CASE(CLFFT_TRANSPOSED_NOTIMPLEMENTED);
		CASE(CLFFT_FILE_NOT_FOUND);
		CASE(CLFFT_FILE_CREATE_FAILURE);
		CASE(CLFFT_VERSION_MISMATCH);
		CASE(CLFFT_INVALID_PLAN);
		CASE(CLFFT_DEVICE_NO_DOUBLE);
		CASE(CLFFT_DEVICE_MISMATCH);
		CASE(CLFFT_ENDSTATUS);
	default:
		return clErrorCodeToString(code);
	}
}

#undef CASE

} // namespace

namespace AlenkaSignal
{

OpenCLContext::OpenCLContext(unsigned int platform, unsigned int device, bool shareCurrentGLContext)
{
	cl_int err;

	// Retrieve the platform and device ids.
	cl_uint pCount = platform + 1;
	cl_platform_id* platforms = new cl_platform_id[pCount];

	err = clGetPlatformIDs(pCount, platforms, &pCount);
	checkClErrorCode(err, "clGetPlatformIDs()");

	if (platform >= pCount)
	{
		stringstream ss;
		ss << "Platform ID " << platform << " too high.";
		throw runtime_error(ss.str());
	}

	platformId = platforms[platform];

	cl_uint dCount = device + 1;
	cl_device_id* devices = new cl_device_id[dCount];

	err = clGetDeviceIDs(platformId, CL_DEVICE_TYPE_ALL, dCount, devices, &dCount);
	checkClErrorCode(err, "clGetDeviceIDs()");

	if (device >= dCount)
	{
		stringstream ss;
		ss << "Device ID " << device << " too high.";
		throw runtime_error(ss.str());
	}

	deviceId = devices[device];

	// Create the context.
	vector<cl_context_properties> properties {CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(platformId)};

	if (shareCurrentGLContext)
	{
		assert(false && "This shouldn't be used yet.");
#if defined WIN_BUILD
		properties.push_back(CL_GL_CONTEXT_KHR);
		properties.push_back(reinterpret_cast<cl_context_properties>(wglGetCurrentContext()));

		properties.push_back(CL_WGL_HDC_KHR);
		properties.push_back(reinterpret_cast<cl_context_properties>(wglGetCurrentDC()));
#elif defined UNIX_BUILD
		properties.push_back(CL_GL_CONTEXT_KHR);
		properties.push_back(reinterpret_cast<cl_context_properties>(glXGetCurrentContext()));

		properties.push_back(CL_GLX_DISPLAY_KHR);
		properties.push_back(reinterpret_cast<cl_context_properties>(glXGetCurrentDisplay()));
#endif
	}

	properties.push_back(0);

	context = clCreateContext(properties.data(), 1, &deviceId, nullptr, nullptr, &err);
	checkClErrorCode(err, "clCreateContext()");

	delete[] platforms;
	delete[] devices;

	//logToFile("OpenCLContext " << this << " created.");
}

OpenCLContext::~OpenCLContext()
{
	cl_int err = clReleaseContext(context);
	checkClErrorCode(err, "clReleaseContext()");

	//logToFile("OpenCLContext " << this << " destroyed.");
}

string OpenCLContext::getPlatformInfo() const
{
	cl_int err;

	cl_uint platformCount;

	err = clGetPlatformIDs(0, nullptr, &platformCount);
	checkClErrorCode(err, "clGetPlatformIDs()");

	cl_platform_id* platformIDs = new cl_platform_id[platformCount];

	err = clGetPlatformIDs(platformCount, platformIDs, nullptr);
	checkClErrorCode(err, "clGetPlatformIDs()");

	// Find the maximum size needed for the values.
	size_t size, maxSize = 0;

	for (cl_uint i = 0; i < platformCount; ++i)
	{
		err = clGetPlatformInfo(platformIDs[i], CL_PLATFORM_VERSION, 0, nullptr, &size);
		checkClErrorCode(err, "clGetPlatformInfo()");
		maxSize = max(maxSize, size);

		err = clGetPlatformInfo(platformIDs[i], CL_PLATFORM_NAME, 0, nullptr, &size);
		checkClErrorCode(err, "clGetPlatformInfo()");
		maxSize = max(maxSize, size);

		err = clGetPlatformInfo(platformIDs[i], CL_PLATFORM_VENDOR, 0, nullptr, &size);
		checkClErrorCode(err, "clGetPlatformInfo()");
		maxSize = max(maxSize, size);

		err = clGetPlatformInfo(platformIDs[i], CL_PLATFORM_EXTENSIONS, 0, nullptr, &size);
		checkClErrorCode(err, "clGetPlatformInfo()");
		maxSize = max(maxSize, size);
	}

	// Build the string.
	char* tmp = new char[maxSize];
	string str;

	str += "Available platforms (" + to_string(platformCount) + "):";
	for (cl_uint i = 0; i < platformCount; ++i)
	{
		if (platformIDs[i] == getCLPlatform())
			str += "\n * ";
		else
			str += "\n   ";

		err = clGetPlatformInfo(platformIDs[i], CL_PLATFORM_NAME, maxSize, tmp, nullptr);
		checkClErrorCode(err, "clGetPlatformInfo()");
		str += tmp;
	}

	str += "\n\nSelected platform: ";
	err = clGetPlatformInfo(getCLPlatform(), CL_PLATFORM_NAME, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetPlatformInfo()");
	str += tmp;

	str += "\nVersion: ";
	err = clGetPlatformInfo(getCLPlatform(), CL_PLATFORM_VERSION, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetPlatformInfo()");
	str += tmp;

	str += "\nVendor: ";
	err = clGetPlatformInfo(getCLPlatform(), CL_PLATFORM_VENDOR, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetPlatformInfo()");
	str += tmp;

	str += "\nExtensions: ";
	err = clGetPlatformInfo(getCLPlatform(), CL_PLATFORM_EXTENSIONS, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetPlatformInfo()");
	str += tmp;

	delete[] tmp;

	return str;
}

string OpenCLContext::getDeviceInfo() const
{
	cl_int err;

	cl_uint deviceCount;

	err = clGetDeviceIDs(getCLPlatform(), CL_DEVICE_TYPE_ALL, 0, nullptr, &deviceCount);
	checkClErrorCode(err, "clGetDeviceIDs()");

	cl_device_id* deviceIDs = new cl_device_id[deviceCount];

	err = clGetDeviceIDs(getCLPlatform(), CL_DEVICE_TYPE_ALL, deviceCount, deviceIDs, nullptr);
	checkClErrorCode(err, "clGetDeviceIDs()");

	// Find the maximum size needed for the value.
	size_t size, maxSize = 0;

	for (cl_uint i = 0; i < deviceCount; ++i)
	{
		err = clGetDeviceInfo(deviceIDs[i], CL_DEVICE_VERSION, 0, nullptr, &size);
		checkClErrorCode(err, "clGetDeviceInfo()");
		maxSize = max(maxSize, size);

		err = clGetDeviceInfo(deviceIDs[i], CL_DEVICE_NAME, 0, nullptr, &size);
		checkClErrorCode(err, "clGetDeviceInfo()");
		maxSize = max(maxSize, size);

		err = clGetDeviceInfo(deviceIDs[i], CL_DEVICE_VENDOR, 0, nullptr, &size);
		checkClErrorCode(err, "clGetDeviceInfo()");
		maxSize = max(maxSize, size);

		err = clGetDeviceInfo(deviceIDs[i], CL_DEVICE_EXTENSIONS, 0, nullptr, &size);
		checkClErrorCode(err, "clGetDeviceInfo()");
		maxSize = max(maxSize, size);
	}

	// Build the string.
	string str;
	char* tmp = new char[maxSize];

	str += "Available devices (" + to_string(deviceCount) + "):";
	for (cl_uint i = 0; i < deviceCount; ++i)
	{
		if (deviceIDs[i] == getCLDevice())
			str += "\n * ";
		else
			str += "\n   ";

		err = clGetDeviceInfo(deviceIDs[i], CL_DEVICE_NAME, maxSize, tmp, nullptr);
		checkClErrorCode(err, "clGetDeviceInfo()");
		str += tmp;
	}

	str += "\n\nSelected device: ";
	err = clGetDeviceInfo(getCLDevice(), CL_DEVICE_NAME, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetDeviceInfo()");
	str += tmp;

	str += "\nVersion: ";
	err = clGetDeviceInfo(getCLDevice(), CL_DEVICE_VERSION, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetDeviceInfo()");
	str += tmp;

	str += "\nVendor: ";
	err = clGetDeviceInfo(getCLDevice(), CL_DEVICE_VENDOR, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetDeviceInfo()");
	str += tmp;

	str += "\nExtensions: ";
	err = clGetDeviceInfo(getCLDevice(), CL_DEVICE_EXTENSIONS, maxSize, tmp, nullptr);
	checkClErrorCode(err, "clGetDeviceInfo()");
	str += tmp;

	str += "\nGlobal memory size: ";
	cl_ulong memSize;
	err = clGetDeviceInfo(getCLDevice(), CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &memSize, nullptr);
	checkClErrorCode(err, "clGetDeviceInfo()");
	stringstream ss;
	double memGigs = memSize/(1000.*1000*1000);
	ss << memGigs << " GB";
	str += ss.str();

	delete[] tmp;

	return str;
}

void OpenCLContext::CCEC(cl_int val, string message, const char* file, int line)
{
	std::stringstream ss;

	ss << "Unexpected error code: ";
	ss << clErrorCodeToString(val);
	ss << ", required CL_SUCCESS. ";

	ss << message << " " << file << ":" << line;

	throw std::runtime_error(ss.str());
}

void OpenCLContext::CFCEC(clfftStatus val, string message, const char* file, int line)
{
	std::stringstream ss;

	ss << "Unexpected error code: " << clfftErrorCodeToString(val);
	ss << ", required CLFFT_SUCCESS. ";

	ss << message << " " << file << ":" << line;

	throw std::runtime_error(ss.str());
}

void OpenCLContext::printBuffer(FILE* file, float* data, int n)
{
#ifndef NDEBUG
	for (int i = 0; i < n; ++i)
	{
		float tmp = data[i];
		if (std::isnan(tmp)/* || tmp < -1000*1000*1000 || tmp > 1000*1000*1000*/)
			tmp = 111111111;
		fprintf(file, "%f\n", tmp);
	}
#else
	(void)file;
	(void)data;
	(void)n;
#endif
}

void OpenCLContext::printBuffer(FILE* file, cl_mem buffer, cl_command_queue queue)
{
#ifndef NDEBUG
	cl_int err;

	size_t size;
	err = clGetMemObjectInfo(buffer, CL_MEM_SIZE, sizeof(size_t), &size, nullptr);
	checkClErrorCode(err, "clGetMemObjectInfo");

	float* tmp = new float[size/sizeof(float)];

	err = clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, size, tmp, 0, nullptr, nullptr);
	checkClErrorCode(err, "clEnqueueReadBuffer");

	printBuffer(file, tmp, size/sizeof(float));

	delete[] tmp;
#else
	(void)file;
	(void)buffer;
	(void)queue;
#endif
}

void OpenCLContext::printBuffer(const string& filePath, float* data, int n)
{
#ifndef NDEBUG
	FILE* file = fopen(filePath.c_str(), "w");
	printBuffer(file, data, n);
	fclose(file);
#else
	(void)filePath;
	(void)data;
	(void)n;
#endif
}

void OpenCLContext::printBuffer(const string& filePath, cl_mem buffer, cl_command_queue queue)
{
#ifndef NDEBUG
	FILE* file = fopen(filePath.c_str(), "w");
	printBuffer(file, buffer, queue);
	fclose(file);
#else
	(void)filePath;
	(void)buffer;
	(void)queue;
#endif
}

void OpenCLContext::printBufferDouble(FILE* file, double* data, int n)
{
#ifndef NDEBUG
	for (int i = 0; i < n; ++i)
	{
		double tmp = data[i];
		if (std::isnan(tmp)/* || tmp < -1000*1000*1000 || tmp > 1000*1000*1000*/)
			tmp = 111111111;
		fprintf(file, "%f\n", tmp);
	}
#else
	(void)file;
	(void)data;
	(void)n;
#endif
}

void OpenCLContext::printBufferDouble(FILE* file, cl_mem buffer, cl_command_queue queue)
{
#ifndef NDEBUG
	cl_int err;

	size_t size;
	err = clGetMemObjectInfo(buffer, CL_MEM_SIZE, sizeof(size_t), &size, nullptr);
	checkClErrorCode(err, "clGetMemObjectInfo");

	double* tmp = new double[size/sizeof(double)];

	err = clEnqueueReadBuffer(queue, buffer, CL_TRUE, 0, size, tmp, 0, nullptr, nullptr);
	checkClErrorCode(err, "clEnqueueReadBuffer");

	printBufferDouble(file, tmp, size/sizeof(double));

	delete[] tmp;
#else
	(void)file;
	(void)buffer;
	(void)queue;
#endif
}

void OpenCLContext::printBufferDouble(const string& filePath, double* data, int n)
{
#ifndef NDEBUG
	FILE* file = fopen(filePath.c_str(), "w");
	printBufferDouble(file, data, n);
	fclose(file);
#else
	(void)filePath;
	(void)data;
	(void)n;
#endif
}

void OpenCLContext::printBufferDouble(const string& filePath, cl_mem buffer, cl_command_queue queue)
{
#ifndef NDEBUG
	FILE* file = fopen(filePath.c_str(), "w");
	printBufferDouble(file, buffer, queue);
	fclose(file);
#else
	(void)filePath;
	(void)buffer;
	(void)queue;
#endif
}

} // namespace AlenkaSignal
