#include "gtest/gtest.h"

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

template<class T>
class DataFileLoader : public SpikedetDataLoader<T>
{
public:
	DataFile* file;

	DataFileLoader(DataFile* file) : file(file)
	{}
	virtual ~DataFileLoader() override
	{}

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
	virtual ~VectorLoader() override
	{}

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
void test(SpikedetDataLoader<T>* loader, double fs, DETECTOR_SETTINGS settings)
{
	clfftStatus errFFT;
	clfftSetupData setupData;

	errFFT = clfftInitSetupData(&setupData);
	checkClfftErrorCode(errFFT, "clfftInitSetupData()");

	errFFT = clfftSetup(&setupData);
	checkClfftErrorCode(errFFT, "clfftSetup()");

	{
		OpenCLContext context(OPENCL_PLATFORM, OPENCL_DEVICE);

		Spikedet<T> det(fs, loader->channelCount(), settings, &context);
		CDetectorOutput* out;
		CDischarges* dis;

		out = new CDetectorOutput;
		dis = new CDischarges(loader->channelCount());
		det.runAnalysis(loader, out, dis);
		unsigned int spikes = out->m_pos.size();
		cerr << "Spikes detected: " << spikes << endl;
		delete out;
		delete dis;

		out = new CDetectorOutput;
		dis = new CDischarges(loader->channelCount());
		det.runAnalysis(loader, out, dis);
		EXPECT_EQ(spikes, out->m_pos.size());
		delete out;
		delete dis;
	}

	errFFT = clfftTeardown();
	checkClfftErrorCode(errFFT, "clfftTeardown()");
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
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P002_default_float)
{
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P003_default_float)
{
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P004_default_float)
{
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P005_default_float)
{
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P006_default_float)
{
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P007_default_float)
{
	EDF file(path + "IED_P001.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

// TODO: make double version of the above tests

TEST_F(spikedet_test, index_bug)
{
	// This tests the bug when computing segment indices.
	EDF file(path + "edfsample.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, zeroChannel_bug0)
{
	// This tests the strange case when you get nan values and it causes an exception.
	EDF file(path + "zeroChannel.edf");
	DataFileLoader<float> loader(&file);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, file.getSamplingFrequency(), defaultSettings); }));
}

TEST_F(spikedet_test, zeroChannel_bug1)
{
	// This test does not appear to be effective, but I will keep it all the same.
	int len = 100000;
	vector<float> signal(3*len);
	for (int i = len; i < 2*len; i++)
		signal[i] = i;

	VectorLoader<float> loader(signal, 3);
	EXPECT_NO_THROW(printException([&] () { test<float>(&loader, 200, defaultSettings); }));
}
