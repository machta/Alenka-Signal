#ifndef CSpikeDetector_H
#define	CSpikeDetector_H

#include <wx/wxprec.h>
#ifndef WX_PRECOMP
#include <wx/wx.h>
#endif
#include <wx/thread.h>

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <numeric>
#include <vector>

#include "CException.h"
#include "CSettingsModel.h"
#include "Definitions.h"
#include "CInputModel.h"
#include "CDSP.h"

#include "lib/Alglib/interpolation.h"

class CDetectorOutput;
class CMarker;
class CDischarges;

/**
 * Implementation of the spike detector.
 * @see http://isarg.feld.cvut.cz/downloads.html
 * @author Jakub Drabek
 */
class CSpikeDetector : public wxThread
{
	// methods
public:
	/**
	 * A constructor
	 * @param frame pointer to object to which are sents events
	 * @param model pointer to model
	 * @param settings pointer to struct \ref DETECTOR_SETTINGS
	 * @param out pointer to output class \ref CDetectorOutput
	 * @param disch pointer to output class \ref CDischarges
	 */
	CSpikeDetector(wxEvtHandler* frame, CInputModel* model, DETECTOR_SETTINGS* settings, CDetectorOutput*& out, CDischarges*& disch);
	
	/**
	 * A virtual destructor.
	 */
	virtual ~CSpikeDetector();
	
	/**
	 * Start the analysis process and store the results.
	 * @throws CException
	 */
	void RunAnalysis();
	
	/**
	 * This is the entry point of the thread.
	 */
	virtual ExitCode Entry();

private:
	/**
	 * Calculate the starts and ends of indexes for CSpikeDetector::spikeDetector
	 * @param indexStart output vector with starts
	 * @param indexStart output vector with ends
	 * @param cntElemInCh count elements in channel
	 * @param T_seg it's equal round(countRecords / N_seg / fs)
	 * @param fs sample rate
	 * @param winsize size of the window
	 */
	void getIndexStartStop(wxVector<int>& indexStart, wxVector<int>& indexStop, const int& cntElemInCh, const double& T_seg, const int& fs, const int& winsize);

	/**
	 * Run analysis for a segment of data.
	 * @param data inpud data - iEEG
	 * @param counChanles count channles of input data
	 * @param inpuFS sample rate of input data
	 * @param bandwidth bandwifth
	 * @param out a pointer to output object of \ref CDetectorOutput
	 * @param discharges a pointer to output object of \ref CDischarges
	 */
	void spikeDetector(wxVector<SIGNALTYPE>*& data, const int& countChannels, const int& inputFS, const BANDWIDTH& bandwidth,
					   CDetectorOutput*& out, CDischarges*& discharges);

	// variables
public:
private:
	/// object to which are sents events
	wxEvtHandler* 	   m_frame;
	/// object for access to data
	CInputModel* 	   m_model;
	/// settings of the detector
	DETECTOR_SETTINGS* m_settings;

	/// output object - structure out
	CDetectorOutput*   m_out;
	/// output object - strucutre discharges
	CDischarges*       m_discharges;

	/// progress of the analysis
	int 			   m_progress;
};


/**
 * Structure containing output data from \ref COneChannelDetect
 */
typedef struct oneChannelDetectRet
{
public:
	wxVector<bool>*   	 m_markersHigh;
	wxVector<bool>*   	 m_markersLow;
	wxVector<double>  	 m_prahInt[2];
	wxVector<double>  	 m_envelopeCdf;
	wxVector<double>  	 m_envelopePdf;
	wxVector<SIGNALTYPE> m_envelope;

	/// A constructor
	oneChannelDetectRet(wxVector<bool>*& markersHigh, wxVector<bool>*& markersLow, const wxVector<double> prahInt[2],
	const wxVector<double> envelopeCdf, const wxVector<double> envelopePdf, const wxVector<SIGNALTYPE>& envelope)
		: m_markersHigh(markersHigh), m_markersLow(markersLow)
	{
		m_prahInt[0].assign(prahInt[0].begin(), prahInt[0].end());
		m_prahInt[1].assign(prahInt[1].begin(), prahInt[1].end());
		m_envelopeCdf.assign(envelopeCdf.begin(), envelopeCdf.end());
		m_envelopePdf.assign(envelopePdf.begin(), envelopePdf.end());
		m_envelope.assign(envelope.begin(), envelope.end());
	}

	/// A destructor
	~oneChannelDetectRet()
	{
		/* empty */
	}

} ONECHANNELDETECTRET;

/**
 * Implementation "one_channel_detect" function from reference implementation of spike detector.
 */
class COneChannelDetect : public wxThread
{
	// methods
public:
	/**
	 * A constructor.
	 * @param data input data
	 * @param settings settings od the detector
	 * @param fs sample rate
	 * @param index indexs
	 * @param channel number of channel
	 */
	COneChannelDetect(const wxVector<SIGNALTYPE>* data, const DETECTOR_SETTINGS* settings, const int& fs, const wxVector<int>* index, const int& channel);

	/**
	 * A virtual desctructor.
	 */
	virtual ~COneChannelDetect();

	/**
	 * This is the entry point of the thread.
	 */
	virtual ExitCode Entry();

private:
	/**
	 * Calculating a mean from data in vector.
	 * @param data a vector of input data
	 * @return a mean
	 */
	double mean(wxVector<double>& data);
	
	/**
	 * Calculating a variance from data in vector.
	 * @param data input vector with data
	 * @param mean mean of data in vector
	 * @return a variance
	 */
	double variance(wxVector<double>& data, const double & mean);

	/**
	 * Detection of local maxima in envelope.
	 * @param envelope envelope of input channel
	 * @param prah_int threeshold curve
	 * @param polyspike_union_time polyspike union time
	 * @return vector cintaining markers of local maxima
	 */
	wxVector<bool>* localMaximaDetection(wxVector<SIGNALTYPE>& envelope, const wxVector<double>& prah_int, const double& polyspike_union_time);
	
	/**
	 * Detecting of union and their merging.
	 * @param marker1 markers of spikes
	 * @param envelope envelope of input channel
	 * @param union_samples union samples time
	 */
	void detectionUnion(wxVector<bool>* marker1, wxVector<SIGNALTYPE>& envelope, const double& union_samples);

	/**
	 * Finding of the highes maxima of the section with local maxima.
	 * implement:
	 * point(:,1)=find(diff([0;marker1])>0); % start
	 * point(:,2)=find(diff([marker1;0])<0); % end
	 * @param point
	 * @param marker1
	 */
	void findStartEndCrossing(wxVector<int> point[2], const wxVector<bool>* marker1);

	// variables
public:
	/* none */

private:
	/// input data
	const wxVector<SIGNALTYPE>* m_data;
	/// settinggs of the detector
	const DETECTOR_SETTINGS*    m_settings;
	/// sample rate
	const int 			  		m_fs;
	/// indexs of suspects areas
	const wxVector<int>* 		m_index;
	/// channel number
	const int 					m_channel;
};

// ------------------------------------------------------------------------------------------------
// OUTPUT CLASSES
// ------------------------------------------------------------------------------------------------

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

#endif
