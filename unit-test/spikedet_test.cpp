#include "gtest/gtest.h"

#include "spikedet.h"
#include "edf.h"

#include <functional>

using namespace std;

namespace
{

template<class T>
class Loader : public SpikedetDataLoader<T>
{
public:
	DataFile* file;

	Loader(const string& filePath)
	{
		file = new EDF(filePath);
	}
	virtual ~Loader() override
	{
		delete file;
	}

	virtual void readSignal(T* data, int64_t firstSample, int64_t lastSample) override
	{
		file->readSignal(data, firstSample, lastSample);
	}
	virtual int64_t sampleCount() override
	{
		return file->getSamplesRecorded();
	}
	virtual int channelCount() override
	{
		return file->getChannelCount();
	}
};

template<class T>
void test(const string& file, DETECTOR_SETTINGS settings)
{
	OpenCLContext context(0, 0);

	Loader<T> loader(file);

	Spikedet<T> det(loader.file->getSamplingFrequency(), loader.channelCount(), settings, &context);
	CDetectorOutput* out = new CDetectorOutput;
	CDischarges* dis = new CDischarges(loader.channelCount());
	det.runAnalysis(&loader, out, dis);

	delete out;
	delete dis;
}

} // namespace

class spikedet_test : public ::testing::Test
{
protected:
	spikedet_test()
	{}

	virtual ~spikedet_test()
	{}

	string path = "unit-test/data/";
	DETECTOR_SETTINGS defaultSettings;
};

TEST_F(spikedet_test, IED_P001_default_float)
{
	test<float>(path + "IED_P001.edf", defaultSettings);
}

TEST_F(spikedet_test, IED_P002_default_float)
{
	test<float>(path + "IED_P002.edf", defaultSettings);
}
