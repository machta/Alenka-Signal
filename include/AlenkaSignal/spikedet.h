#ifndef ALENKASIGNAL_SPIKEDET_H
#define ALENKASIGNAL_SPIKEDET_H

#include <CL/cl_gl.h>

#include <cstdio>
#include <cmath>
#include <vector>
#include <atomic>

#define wxVector std::vector

#define SIGNALTYPE float

namespace AlenkaSignal
{

class OpenCLContext;
template <class T>
class FilterProcessor;

template<class T>
class SpikedetDataLoader
{
public:
	virtual ~SpikedetDataLoader() {}

	virtual void readSignal(T* data, int64_t firstSample, int64_t lastSample) = 0;
	virtual int64_t sampleCount() = 0;
	virtual int channelCount() = 0;
};

typedef struct bandwidth
{
	/**
	 * A constructor
	 * @param bl lower limit filtering
	 * @param bh upper limit filtering
	 */
	bandwidth(const int& bl, const int& bh) : m_bandLow(bl), m_bandHigh(bh) {}

	/// Lower limit of filtering.
	int m_bandLow;

	/// Upper limit filtering.
	int m_bandHigh;
} BANDWIDTH;

/**
 * Output class containing output data from the detector.
 */
class CDetectorOutput
{
public:
	/**
	 * A constructor.
	 */
	CDetectorOutput();

	/**
	 * A virtual desctructor.
	 */
	virtual ~CDetectorOutput();

	/**
	 * Add data to the vectors.
	 * @param pos spike position (second).
	 * @param dur spike duration (second) - fix value 5 ms.
	 * @param chan channel.
	 * @param con spike condition (1-obvious 0.5-ambiguous).
	 * @param weight statistical significance "CDF".
	 * @param pdf sstatistical significance "PDF".
	 */
	void Add(const double& pos, const double& dur, const int& chan, const double& con, const double& weight, const double& pdf);

	/**
	 * Erase records at positions.
	 * @param pos position of records to erase.
	 */
	void Remove(const wxVector<int>& pos);

	/// spike position (second)
	wxVector<double>  m_pos;
	/// channel
	wxVector<int>     m_chan;
	/// spike duration (second) - fix value 5 ms
	wxVector<double>  m_dur;
	/// spike condition (1-obvious 0.5-ambiguous)
	wxVector<double>  m_con;
	/// statistical significance "CDF"
	wxVector<double>  m_weight;
	/// statistical significance "PDF"
	wxVector<double>  m_pdf;
};

/**
 * Discharges
 */
class CDischarges
{
public:
	/**
	 * A constructor.
	 * @param countChannels count channels.
	 */
	CDischarges(const int& countChannels);

	/**
	 * A virtual desctructor.
	 */
	virtual ~CDischarges();

	/**
	 * Erase records at positions.
	 * @param pos positions of record to erase.
	 */
	void Remove(const wxVector<int>& pos);

	/**
	 *	Return count channels.
	 * @return count channels.s
	 */
	inline unsigned GetCountChannels() const
	{
		return m_countChannels;
	}

	/// spike type 1-obvious, 0.5- ambiguous
	std::vector<double>* m_MV;
	/// max. amplitude of envelope above backround
	std::vector<double>* m_MA;
	/// event start position
	std::vector<double>* m_MP;
	/// duration of event
	std::vector<double>* m_MD;
	/// statistical significance "CDF"
	std::vector<double>* m_MW;
	/// probability of occurence
	std::vector<double>* m_MPDF;

private:
	/// count channels
	unsigned 			 m_countChannels;
};

typedef struct detectorSettings
{
	int    m_band_low = 10;                // -fl
	int    m_band_high = 60;               // -fh
	double m_k1 = 3.65;                    // -k1
	double m_k2 = 3.65;                    // -k2
	double m_k3 = 0;                       // -k3
	int    m_winsize = 5;                  // -w
	double m_noverlap = 4;                 // -n
	int    m_buffering = 300;              // -buf
	int    m_main_hum_freq = 50;           // -h
	double m_discharge_tol = 0.005;        // -dt
	double m_polyspike_union_time = 0.12;  // -pt
	int    m_decimation = 200;             // -dec

	/// A constructor
	detectorSettings() {}
	detectorSettings(int band_low, int band_high, double k1, double k2, double k3, int winsize, double noverlap, int buffering, int main_hum_freq,
		double discharge_tol, double polyspike_union_time, int decimation)
		: m_band_low(band_low), m_band_high(band_high), m_k1(k1), m_k2(k2), m_k3(k3), m_winsize(winsize), m_noverlap(noverlap), m_buffering(buffering),
		m_main_hum_freq(main_hum_freq), m_discharge_tol(discharge_tol), m_polyspike_union_time(polyspike_union_time), m_decimation(decimation) {}
} DETECTOR_SETTINGS;

template<class T>
class Spikedet
{
	int fs;
	int channelCount;
	FilterProcessor<T>* filterProcessor = nullptr;
	int decimationF;
	cl_command_queue queue = nullptr;
	cl_mem inBuffer = nullptr, outBuffer = nullptr;
	std::vector<T> segmentBuffer, stepBuffer;
	int progressComplete;
	std::atomic<int> progressCurrent;
	std::atomic<bool> cancelComputation;

	DETECTOR_SETTINGS settings;
	OpenCLContext* context;
	DETECTOR_SETTINGS* m_settings = &settings;
	/// output object - structure out
	CDetectorOutput*   m_out;
	/// output object - strucutre discharges
	CDischarges*       m_discharges;

public:
	Spikedet(int fs, int channelCount, DETECTOR_SETTINGS settings, OpenCLContext* context);
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

private:
	void getIndexStartStop(wxVector<int64_t>& indexStart, wxVector<int64_t>& indexStop, int64_t cntElemInCh, int64_t T_seg, int fs, int winsize);

	void spikeDetector(SpikedetDataLoader<T>* loader, int startSample, int stopSample, const int& countChannels, const int& inputFS, const BANDWIDTH& bandwidth,
	                   CDetectorOutput*& out, CDischarges*& discharges);

	std::vector<T>* prepareSegment(SpikedetDataLoader<T>* loader, int start, int stop);
};

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_SPIKEDET_H
