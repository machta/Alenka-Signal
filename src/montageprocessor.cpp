#include <AlenkaSignal/montageprocessor.h>

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/montage.h>

using namespace std;

namespace AlenkaSignal
{

template<class T>
void MontageProcessor<T>::process(const vector<Montage<T>*>& montage, cl_mem inBuffer, cl_mem outBuffer, cl_command_queue queue, int inputRowOffset)
{
	cl_int err;

	for (unsigned int i = 0; i < montage.size(); i++)
	{
		cl_kernel montageKernel = montage[i]->getKernel();
		int pi = 0;

		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_mem), &inBuffer);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_mem), &outBuffer);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_int), &inputRowLength);
		checkClErrorCode(err, "clSetKernelArg()");

		cl_int offset = inputRowOffset;
		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_int), &offset);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_int), &inputRowCount);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_int), &outputRowLength);
		checkClErrorCode(err, "clSetKernelArg()");

		cl_int index = i;
		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_int), &index);
		checkClErrorCode(err, "clSetKernelArg()");

		err = clSetKernelArg(montageKernel, pi++, sizeof(cl_int), &outputCopyCount);
		checkClErrorCode(err, "clSetKernelArg()");

		size_t globalWorkSize = outputRowLength;

		err = clEnqueueNDRangeKernel(queue, montageKernel, 1, nullptr, &globalWorkSize, nullptr, 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueNDRangeKernel()");
	}
}

template class MontageProcessor<float>;
template class MontageProcessor<double>;

} // namespace AlenkaSignal
