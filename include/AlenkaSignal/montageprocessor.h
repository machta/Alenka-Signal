#ifndef ALENKASIGNAL_MONTAGEPROCESSOR_H
#define ALENKASIGNAL_MONTAGEPROCESSOR_H

#include <CL/cl_gl.h>

#include <vector>

namespace AlenkaSignal
{

template <class T>
class Montage;

/**
 * @brief This class handles computation of montages.
 */
template<class T>
class MontageProcessor
{
	cl_int inputRowLength, outputRowLength, inputRowCount, outputCopyCount;

public:
	/**
	 * @brief MontageProcessor constructor.
	 * @param offset Skip this many samples at the beginning of the input buffer.
	 */
	MontageProcessor(unsigned int inputRowLength, unsigned int outputRowLength, int inputRowCount, int outputCopyCount = 1) :
		inputRowLength(inputRowLength), outputRowLength(outputRowLength), inputRowCount(inputRowCount), outputCopyCount(outputCopyCount) {}

	/**
	 * @brief Enqueues all commands required for montage computation to queue.
	 */
	void process(const std::vector<Montage<T>*>& montage, cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue, int inputRowOffset = 0);
};

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_MONTAGEPROCESSOR_H
