#include <gtest/gtest.h>

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/spikedet.h>
#include <AlenkaFile/edf.h>

#include <functional>
#include <memory>

using namespace std;
using namespace AlenkaSignal;
using namespace AlenkaFile;

namespace
{

const string PATH = "unit-test/data/";

class DataFileLoader : public SpikedetDataLoader
{
public:
	DataFile* file;

	DataFileLoader(DataFile* file) : file(file) {}

	virtual void readSignal(SIGNALTYPE* data, int64_t firstSample, int64_t lastSample) override
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

class VectorLoader : public SpikedetDataLoader
{
public:
	vector<SIGNALTYPE> signal;
	int channels;
	int length;

	VectorLoader(vector<SIGNALTYPE> signal, int channels) : signal(signal), channels(channels), length(signal.size()/channels)
	{
		assert(signal.size() == channels*length);
	}

	virtual void readSignal(SIGNALTYPE* data, int64_t firstSample, int64_t lastSample) override
	{
		int64_t len = lastSample - firstSample + 1;
		for (int j = 0; j < channels; j++)
			for (int64_t i = firstSample; i <= lastSample; i++)
			{
				SIGNALTYPE sample = i < 0 || i >= length ? 0 : signal.at(j*length + i);

				data[j*len + i - firstSample] = sample;
			}
	}
	virtual int64_t sampleCount() override
	{
		return length;
	}
	virtual int channelCount() override
	{
		return channels;
	}
};

int test(SpikedetDataLoader* loader, double fs, bool original = false, DETECTOR_SETTINGS settings = Spikedet::defaultSettings())
{
	unique_ptr<CDetectorOutput> out(new CDetectorOutput);
	unique_ptr<CDischarges> dis(new CDischarges(loader->channelCount()));

	Spikedet det(fs, loader->channelCount(), original, settings);
	det.runAnalysis(loader, out.get(), dis.get());

	int spikes = static_cast<int>(out->m_pos.size());
	EXPECT_EQ(spikes, out->m_pos.size());
	return spikes;
}

void printException(function<void (void)> fun)
{
	try
	{
		fun();
	}
	catch (exception& e)
	{
		cerr << "Caught an std exception: " << e.what() << endl;
		throw;
	}
	catch (...)
	{
		cerr << "Caught an exception." << endl;
		throw;
	}
}

const int spikeCounts[7] = {894, 497, 382, 22, 1367, 648, 439};
const double ALLOWED_ERROR = 1;
const double ALLOWED_ERROR_ORIGINAL = 0;

double relativeError(int res, int sol)
{
	return abs(res - sol)*100./sol;
}

} // namespace

// TODO: Add test that repeatedly runs the analysis using lingle instance of Spikedet class.

TEST(spikedet_test, IED_P001_default)
{
	EDF file(PATH + "IED_P001.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[0]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[0]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P002_default)
{
	EDF file(PATH + "IED_P002.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[1]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[1]), ALLOWED_ERROR_ORIGINAL*2);
}

TEST(spikedet_test, IED_P003_default)
{
	EDF file(PATH + "IED_P003.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[2]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[2]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P004_default)
{
	EDF file(PATH + "IED_P004.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[3]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[3]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P005_default)
{
	EDF file(PATH + "IED_P005.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[4]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[4]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P006_default)
{
	EDF file(PATH + "IED_P006.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[5]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[5]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P007_default)
{
	EDF file(PATH + "IED_P007.edf");
	DataFileLoader loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency()); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[6]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test(&loader, file.getSamplingFrequency(), true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[6]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, index_bug)
{
	// This tests the bug when computing segment indices.
	EDF file(PATH + "edfsample.edf");
	DataFileLoader loader(&file);
	EXPECT_NO_THROW(printException([&] () { test(&loader, file.getSamplingFrequency()); }));
}

TEST(spikedet_test, zeroChannel_bug0)
{
	// This tests the strange case when you get nan values and it causes an exception.
	EDF file(PATH + "zeroChannel.edf");
	DataFileLoader loader(&file);
	EXPECT_NO_THROW(printException([&] () { test(&loader, file.getSamplingFrequency()); }));
}

TEST(spikedet_test, zeroChannel_bug1)
{
	// This test does not appear to be effective in reproducing the bug, but I will keep it all the same.
	int len = 100000;
	vector<SIGNALTYPE> signal(3*len);
	for (int i = len; i < 2*len; i++)
		signal[i] = i;

	VectorLoader loader(signal, 3);
	EXPECT_NO_THROW(printException([&] () { test(&loader, 200); }));
}
