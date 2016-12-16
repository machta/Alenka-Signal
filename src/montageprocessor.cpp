#include <AlenkaSignal/montageprocessor.h>

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/montage.h>

using namespace std;

namespace AlenkaSignal
{

template<class T>
MontageProcessor<T>::MontageProcessor(unsigned int offset, unsigned int blockLength, int channelsInFile) :
	inputRowLength(offset + blockLength), inputRowOffset(offset), outputRowLength(blockLength), channelsInFile(channelsInFile)
{
}

template<class T>
MontageProcessor<T>::~MontageProcessor()
{
}

template<class T>
void MontageProcessor<T>::process(const vector<Montage<T>*>& montage, cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue)
{
	cl_int err;

	for (unsigned int i = 0; i < montage.size(); i++)
	{
		cl_kernel montageKernel = montage[i]->getKernel();

		err = clSetKernelArg(montageKernel, 0, sizeof(cl_mem), &inBuffer);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, 1, sizeof(cl_mem), &outBuffer);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, 2, sizeof(cl_int), &inputRowLength);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, 3, sizeof(cl_int), &inputRowOffset);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, 4, sizeof(cl_int), &outputRowLength);
		checkClErrorCode(err, "clSetKernelArg()");

		cl_int index = i;
		err = clSetKernelArg(montageKernel, 5, sizeof(cl_int), &index);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, 6, sizeof(cl_int), &channelsInFile);
		checkClErrorCode(err, "clSetKernelArg()");

		size_t globalWorkSize = outputRowLength;

		err = clEnqueueNDRangeKernel(queue, montageKernel, 1, nullptr, &globalWorkSize, nullptr, 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueNDRangeKernel()");
	}
}

template class MontageProcessor<float>;
template class MontageProcessor<double>;

} // namespace AlenkaSignal
