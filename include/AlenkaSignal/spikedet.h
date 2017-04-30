#ifndef ALENKASIGNAL_SPIKEDET_H
#define ALENKASIGNAL_SPIKEDET_H

#include "../../spikedet/src/Definitions.h"
#include "../../spikedet/src/spikedetoutput.h"
#include "../../spikedet/src/CSettingsModel.h"

#include <vector>
#include <atomic>

class CSpikeDetector;

namespace AlenkaSignal
{

class SpikedetDataLoader
{
public:
	virtual ~SpikedetDataLoader() {}

	virtual void readSignal(SIGNALTYPE* data, int64_t firstSample, int64_t lastSample) = 0;
	virtual int64_t sampleCount() = 0;
	virtual int channelCount() = 0;
};

class Spikedet
{
	const int fs;
	int channelCount;
	bool original;
	DETECTOR_SETTINGS settings;
	std::atomic<int> progressCurrent;
	CSpikeDetector* detector = nullptr;

public:
	Spikedet(int fs, int channelCount, bool original, DETECTOR_SETTINGS settings);
	~Spikedet();

	void runAnalysis(SpikedetDataLoader* loader, CDetectorOutput* out, CDischarges* discharges);

	/**
	 * @brief progressPercentage is used to query the completion status of the analysis.
	 *
	 * This method is thread safe. It can be called from a thread other than the one
	 * that launched the detector via runAnalysis.
	 * @return Returns the percentage towards completion of the operation.
	 */
	int progressPercentage() const
	{
		return progressCurrent;
	}

	/**
	 * @brief Use cancel to tell the detector to quit the computation at the earliest oppurtunity.
	 *
	 * Thread safe.
	 */
	void cancel();

	static DETECTOR_SETTINGS defaultSettings()
	{
		return DETECTOR_SETTINGS(10, 60, 3.65, 3.65, 0, 5, 4, 300, 50, 0.005, 0.12, 200);
	}
};

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_SPIKEDET_H
