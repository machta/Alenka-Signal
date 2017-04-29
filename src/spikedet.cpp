#include "../include/AlenkaSignal/spikedet.h"

using namespace std;
using namespace AlenkaSignal;

namespace
{

class InputModel : public CInputModel
{
	SpikedetDataLoader<SIGNALTYPE>* loader;
	vector<SIGNALTYPE> buffer;

public:
	InputModel(SpikedetDataLoader<SIGNALTYPE>* loader) : loader(loader)
	{
		m_fs = loader->channelCount();
		m_countSamples = loader->sampleCount();
		m_channels.resize(loader->channelCount());
	}

	virtual void OpenFile(const char* fileName) override {}
	virtual void OpenFile(const wchar_t* fileName) override {}
	virtual void CloseFile() override {}
	virtual bool IsOpen() const override { return true; }
	virtual bool IsEnd() const override { return false; }

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

template<class T>
Spikedet<T>::Spikedet(int fs, int channelCount, bool originalDecimation, DETECTOR_SETTINGS settings)
	: fs(fs), channelCount(channelCount), originalDecimation(originalDecimation), settings(settings)
{

}

template<class T>
Spikedet<T>::~Spikedet()
{

}

template<class T>
void Spikedet<T>::runAnalysis(SpikedetDataLoader<T>* loader, CDetectorOutput*& out, CDischarges*& discharges)
{
	unique_ptr<CInputModel> model(new InputModel(loader));

	unique_ptr<CSpikeDetector> detector(new CSpikeDetector(nullptr, model.get(), &settings, out, discharges));

	detector->Run();
}

template class Spikedet<float>;
//template class NewSpikedet<double>;

} // namespace AlenkaSignal
