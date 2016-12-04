#ifndef MONTAGEPROCESSOR_H
#define MONTAGEPROCESSOR_H

#include "montage.h"
#include "openclprogram.h"

#include <CL/cl_gl.h>

#include <vector>

/**
 * @brief This class handles computation of montages.
 */
template<class T>
class MontageProcessor
{
public:
	/**
	 * @brief MontageProcessor constructor.
	 * @param offset Skip this many samples at the beginning of the input buffer.
	 */
	MontageProcessor(unsigned int offset, unsigned int blockLength, int channelsInFile);
	~MontageProcessor();

	/**
	 * @brief Enqueues all commands required for montage computation to queue.
	 */
	void process(const std::vector<Montage<T>*>& montage, cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue);
	
	/**
	 * @brief Returns the number of tracks of the montage. 
	 */
	unsigned int getNumberOfRows() const // TODO: remove if not used
	{
		return numberOfRows;
	}

private:
	cl_int inputRowLength;
	cl_int inputRowOffset;
	cl_int outputRowLength;
	cl_int channelsInFile;
	cl_kernel montageKernel = nullptr;
	unsigned int numberOfRows;

	void releaseMontage()
	{
		cl_int err;

		if (montageKernel != nullptr)
		{
			err = clReleaseKernel(montageKernel);
			checkClErrorCode(err, "clReleaseKernel()");
		}
	}
};

#endif // MONTAGEPROCESSOR_H
