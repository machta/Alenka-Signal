#include "gtest/gtest.h"

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/spikedet.h>
#include <Alenka-File/edf.h>

#include <functional>
#include <memory>

using namespace std;
using namespace AlenkaSignal;
using namespace AlenkaFile;

namespace
{

template<class T>
class Loader : public SpikedetDataLoader<T>
{
public:
	DataFile* file;

	Loader(DataFile* file) : file(file)
	{}
	virtual ~Loader() override
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
void test(DataFile* file, DETECTOR_SETTINGS settings)
{
	clfftStatus errFFT;
	clfftSetupData setupData;

	errFFT = clfftInitSetupData(&setupData);
	checkClfftErrorCode(errFFT, "clfftInitSetupData()");

	errFFT = clfftSetup(&setupData);
	checkClfftErrorCode(errFFT, "clfftSetup()");

	{
		OpenCLContext context(OPENCL_PLATFORM, OPENCL_DEVICE);

		Loader<T> loader(file);

		Spikedet<T> det(loader.file->getSamplingFrequency(), loader.channelCount(), settings, &context);
		CDetectorOutput* out = new CDetectorOutput;
		CDischarges* dis = new CDischarges(loader.channelCount());
		det.runAnalysis(&loader, out, dis);

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
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P002_default_float)
{
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P003_default_float)
{
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P004_default_float)
{
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P005_default_float)
{
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P006_default_float)
{
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

TEST_F(spikedet_test, IED_P007_default_float)
{
	unique_ptr<DataFile> file(new EDF(path + "IED_P001.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}

// TODO: make double version of the above tests

TEST_F(spikedet_test, bug)
{
	unique_ptr<DataFile> file(new EDF(path + "edfsample.edf"));
	EXPECT_NO_THROW(printException([&] () { test<float>(file.get(), defaultSettings); }));
}
