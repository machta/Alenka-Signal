#include "../include/AlenkaSignal/spikedet.h"

#include "../../spikedet/src/CSpikeDetector.h"

using namespace std;
using namespace AlenkaSignal;

namespace
{

class InputModel : public CInputModel
{
	SpikedetDataLoader* loader;
	vector<SIGNALTYPE> buffer;

public:
	InputModel(int fs, SpikedetDataLoader* loader) : loader(loader)
	{
		m_fs = fs;
		m_countSamples = loader->sampleCount();
		m_channels.resize(loader->channelCount());
	}

	virtual void OpenFile(const char* fileName) override { assert(0); }
	virtual void OpenFile(const wchar_t* fileName) override { assert(0); }
	virtual void CloseFile() override { assert(0); }
	virtual bool IsOpen() const override { assert(0); return true; }
	virtual bool IsEnd() const override { assert(0); return false; }

	virtual wxVector<SIGNALTYPE>* GetSegment(const int& start, const int& end) override
	{
		int channelCount = GetCountChannels();
		int size = end - start;

		buffer.resize(channelCount*size);
		loader->readSignal(buffer.data(), start, end - 1);

		auto channels = new wxVector<SIGNALTYPE>[channelCount];

		for (int i = 0; i < channelCount; ++i)
		{
			channels[i].resize(size);

			SIGNALTYPE* channelsPointer = channels[i].data();
			SIGNALTYPE* bufferPointer = buffer.data() + i*size;

			for (int j = 0; j < size; ++j)
				channelsPointer[j] = bufferPointer[j];
		}

		return channels;
	}
};

} // namespace

namespace AlenkaSignal
{

Spikedet::Spikedet(int fs, int channelCount, bool originalDecimation, DETECTOR_SETTINGS settings)
	: fs(fs), channelCount(channelCount), originalDecimation(originalDecimation), settings(settings)
{

}

Spikedet::~Spikedet()
{

}

void Spikedet::runAnalysis(SpikedetDataLoader* loader, CDetectorOutput* out, CDischarges* discharges)
{
	unique_ptr<CInputModel> model(new InputModel(fs, loader));

	unique_ptr<CSpikeDetector> detector(new CSpikeDetector(nullptr, model.get(), &settings, out, discharges));

	detector->Run();
}

} // namespace AlenkaSignal
