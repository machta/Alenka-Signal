#ifndef SPIKEDET_H
#define SPIKEDET_H

#include "openclcontext.h"

#include <cstdio>
#include <cmath>
#include <vector>

#include <CL/cl_gl.h>
#include <clFFT.h>

#define wxVector std::vector

#define SIGNALTYPE float

template<class T>
class SpikedetDataLoader
{
public:
	void readSignal(T* data, int64_t firstSample, int64_t lastSample) = 0;
	int64_t sampleCount() = 0;
	int channelCount() = 0;
};

typedef struct bandwidth
{
public:
	/**
	 * A constructor
	 * @param bl lower limit filtering
	 * @param bh upper limit filtering
	 */
	bandwidth(const int& bl, const int& bh)
	    : m_bandLow(bl), m_bandHigh(bh)
	{
		/* empty */
	}

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
// methods
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
private:

// variables
public:
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
private:
	/* none */
};

/**
 * Discharges
 */
class CDischarges
{
// methods
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
private:

// variables
public:
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

/*struct DETECTOR_SETTINGS
{
	int    band_low = 10;          // -fl
	int    band_high = 60;         // -fh
	double k1 = 3.65;              // -k1
	double k2 = -1000;             // -k2
	double k3 = 0;                 // -k3
	int    winsize = -1000;        // -w
	double noverlap = -1000;       // -n
	int    buffering = 300;        // -buf
	int    main_hum_freq = 50;     // -h
	double discharge_tol;          // -dt
	double polyspike_union_time;   // -pt
	int    decimation = 200;       // -dec
};*/

typedef struct detectorSettings {
public:
	int    m_band_low;                // (-fl)
	int    m_band_high;               // (-fh)
	double m_k1;               		  // (-k1)
	double m_k2;               		  // (-k2)
	double m_k3;					  // (-k3)
	int    m_winsize;     			  // (-w)
	double m_noverlap;     			  // (-n)
	int    m_buffering;               // (-buf)
	int    m_main_hum_freq;           // (-h)
	double m_discharge_tol;           // (-dt)
	double m_polyspike_union_time ;	  // (-pt)
	int    m_decimation;

	/// A constructor
	detectorSettings(int band_low, int band_high, double k1, double k2, double k3, int winsize, double noverlap, int buffering, int main_hum_freq,
	    double discharge_tol, double polyspike_union_time, int decimation)
	    : m_band_low(band_low), m_band_high(band_high), m_k1(k1), m_k2(k2), m_k3(k3), m_winsize(winsize), m_noverlap(noverlap), m_buffering(buffering),
	    m_main_hum_freq(main_hum_freq), m_discharge_tol(discharge_tol), m_polyspike_union_time(polyspike_union_time), m_decimation(decimation)
	{
		/* empty */
	}

} DETECTOR_SETTINGS;

template<class T>
class Spikedet
{
public:
	Spikedet(int fs, DETECTOR_SETTINGS settings, OpenCLContext* context) : fs(fs), settings(settings), context(context)
	{
		/*if (settings.k2 < 0)
			this->settings.k2 = settings.k1;

		if (settings.winsize < 0)
			this->settings.winsize = 5*fs;

		if (settings.noverlap < 0)
			this->settings.noverlap = 4*fs;*/
	}

	~Spikedet()
	{}

	void runAnalysis(SpikedetDataLoader<T>* loader, CDetectorOutput* out, CDischarges* discharges);

private:
	int fs;
	DETECTOR_SETTINGS settings;
	OpenCLContext* context;
	DETECTOR_SETTINGS* m_settings = &settings;
	/// output object - structure out
	CDetectorOutput*   m_out;
	/// output object - strucutre discharges
	CDischarges*       m_discharges;

	void getIndexStartStop(wxVector<int>& indexStart, wxVector<int>& indexStop, const int& cntElemInCh, const double& T_seg, const int& fs, const int& winsize);

	void spikeDetector(wxVector<SIGNALTYPE>& data, int countRecords, const int& countChannels, const int& inputFS, const BANDWIDTH& bandwidth,
	                   CDetectorOutput* out, CDischarges* discharges);
};



#endif // SPIKEDET_H
