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
const DETECTOR_SETTINGS defaultSettings;

template<class T>
class DataFileLoader : public SpikedetDataLoader<T>
{
public:
	DataFile* file;

	DataFileLoader(DataFile* file) : file(file) {}

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
class VectorLoader : public SpikedetDataLoader<T>
{
public:
	vector<T> signal;
	int channels;
	int length;

	VectorLoader(vector<T> signal, int channels) : signal(signal), channels(channels), length(signal.size()/channels)
	{
		assert(signal.size() == channels*length);
	}

	virtual void readSignal(T* data, int64_t firstSample, int64_t lastSample) override
	{
		int64_t len = lastSample - firstSample + 1;
		for (int j = 0; j < channels; j++)
			for (int64_t i = firstSample; i <= lastSample; i++)
			{
				T sample = i < 0 || i >= length ? 0 : signal.at(j*length + i);

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

template<class T>
int test(SpikedetDataLoader<T>* loader, double fs, DETECTOR_SETTINGS settings, bool originalDecimation = false)
{
	OpenCLContext::clfftInit();

	int spikes;
	{
		OpenCLContext context(OPENCL_PLATFORM, OPENCL_DEVICE);

		Spikedet<T> det(fs, loader->channelCount(), originalDecimation, settings, &context);
		CDetectorOutput* out;
		CDischarges* dis;

		out = new CDetectorOutput;
		dis = new CDischarges(loader->channelCount());
		det.runAnalysis(loader, out, dis);
		spikes = static_cast<int>(out->m_pos.size());
		delete out;
		delete dis;

		out = new CDetectorOutput;
		dis = new CDischarges(loader->channelCount());
		det.runAnalysis(loader, out, dis);
		EXPECT_EQ(spikes, out->m_pos.size());
		delete out;
		delete dis;
	}

	OpenCLContext::clfftDeinit();

	cerr << "Spikes detected: " << spikes << endl;
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

const int spikeCounts[7] = {894, 497, 382, 22, 1366, 646, 438};
const double ALLOWED_ERROR = 1;
const double ALLOWED_ERROR_ORIGINAL = 0.35;

double relativeError(int res, int sol)
{
	return abs(res - sol)*100./sol;
}

} // namespace

TEST(spikedet_test, IED_P001_default_float)
{
	EDF file(PATH + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[0]), 1.5);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[0]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P001_default_double)
{
	EDF file(PATH + "IED_P001.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[0]), 1.5);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[0]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P002_default_float)
{
	EDF file(PATH + "IED_P002.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[1]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[1]), ALLOWED_ERROR_ORIGINAL*2);
}

TEST(spikedet_test, IED_P002_default_double)
{
	EDF file(PATH + "IED_P002.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[1]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[1]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P003_default_float)
{
	EDF file(PATH + "IED_P003.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[2]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[2]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P003_default_double)
{
	EDF file(PATH + "IED_P003.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[2]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[2]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P004_default_float)
{
	EDF file(PATH + "IED_P004.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[3]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[3]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P004_default_double)
{
	EDF file(PATH + "IED_P004.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[3]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[3]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P005_default_float)
{
	EDF file(PATH + "IED_P005.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[4]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[4]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P005_default_double)
{
	EDF file(PATH + "IED_P005.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[4]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[4]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P006_default_double)
{
	EDF file(PATH + "IED_P006.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[5]), ALLOWED_ERROR);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[5]), ALLOWED_ERROR_ORIGINAL);
}

TEST(spikedet_test, IED_P006_default_float)
{
	EDF file(PATH + "IED_P006.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[5]), 1);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[5]), 1);
}

TEST(spikedet_test, IED_P007_default_float)
{
	EDF file(PATH + "IED_P007.edf");
	DataFileLoader<float> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[6]), 1);

	EXPECT_NO_THROW(printException([&] () { spikes = test<float>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[6]), 1);
}

TEST(spikedet_test, IED_P007_default_double)
{
	EDF file(PATH + "IED_P007.edf");
	DataFileLoader<double> loader(&file);
	int spikes;

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[6]), 1);

	EXPECT_NO_THROW(printException([&] () { spikes = test<double>(&loader, file.getSamplingFrequency(), defaultSettings, true); }));
	EXPECT_LE(relativeError(spikes, spikeCounts[6]), 1);
}

TEST(spikedet_test, index_bug_float)
{
	// This tests the bug when computing segment indices.
	EDF file(PATH + "edfsample.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST(spikedet_test, index_bug_double)
{
	// This tests the bug when computing segment indices.
	EDF file(PATH + "edfsample.edf");
	DataFileLoader<double> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST(spikedet_test, zeroChannel_bug0_float)
{
	// This tests the strange case when you get nan values and it causes an exception.
	EDF file(PATH + "zeroChannel.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST(spikedet_test, zeroChannel_bug0_double)
{
	// This tests the strange case when you get nan values and it causes an exception.
	EDF file(PATH + "zeroChannel.edf");
	DataFileLoader<double> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<double>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST(spikedet_test, zeroChannel_bug1_float)
{
	// This test does not appear to be effective in reproducing the bug, but I will keep it all the same.
	int len = 100000;
	vector<float> signal(3*len);
	for (int i = len; i < 2*len; i++)
		signal[i] = i;

	VectorLoader<float> loader(signal, 3);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, 200, defaultSettings); }));
}

TEST(spikedet_test, zeroChannel_bug1_double)
{
	// This test does not appear to be effective in reproducing the bug, but I will keep it all the same.
	int len = 100000;
	vector<double> signal(3*len);
	for (int i = len; i < 2*len; i++)
		signal[i] = i;

	VectorLoader<double> loader(signal, 3);
	EXPECT_NO_THROW(printException([&] () { test<double>(&loader, 200, defaultSettings); }));
}
