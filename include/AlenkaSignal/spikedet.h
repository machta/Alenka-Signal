#ifndef ALENKASIGNAL_SPIKEDET_H
#define ALENKASIGNAL_SPIKEDET_H

#include <CSpikeDetector.h>

#include <vector>
#include <atomic>

namespace AlenkaSignal
{

template<class T>
class SpikedetDataLoader
{
public:
	virtual ~SpikedetDataLoader() {}

	virtual void readSignal(T* data, int64_t firstSample, int64_t lastSample) = 0;
	virtual int64_t sampleCount() = 0;
	virtual int channelCount() = 0;
};

template<class T>
class Spikedet
{
	const int fs;
	int channelCount;
	int progressComplete;
	std::atomic<int> progressCurrent;
	std::atomic<bool> cancelComputation;
	bool originalDecimation;

	DETECTOR_SETTINGS settings;

public:
	Spikedet(int fs, int channelCount, bool originalDecimation, DETECTOR_SETTINGS settings);
	~Spikedet();

	void runAnalysis(SpikedetDataLoader<T>* loader, CDetectorOutput*& out, CDischarges*& discharges);

	/**
	 * @brief progressPercentage is used to query the completion status of the analysis.
	 *
	 * This method is thread safe. It can be called from a thread other than the one
	 * that launched the detector via runAnalysis.
	 * @return Returns the percentage towards completion of the operation.
	 */
	int progressPercentage() const
	{
		return 100*progressCurrent/progressComplete;
	}

	/**
	 * @brief Use cancel to tell the detector to quit the computation at the earliest oppurtunity.
	 *
	 * Thread safe.
	 */
	void cancel()
	{
		cancelComputation = true;
	}

	static DETECTOR_SETTINGS defaultSettings()
	{
		return DETECTOR_SETTINGS(10, 60, 3.65, 3.65, 0, 5, 4, 300, 50, 0.005, 0.12, 200);
	}

private:

};

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_SPIKEDET_H
