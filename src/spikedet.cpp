#include <AlenkaSignal/spikedet.h>

#include <AlenkaSignal/openclcontext.h>
#include <AlenkaSignal/filter.h>
#include <AlenkaSignal/filterprocessor.h>

#include "interpolation.h"
#include "fasttransforms.h"

#include <climits>
#include <algorithm>
#include <complex>
#include <numeric>

#include <Eigen/Dense>
using namespace std;
using namespace AlenkaSignal;

namespace
{

namespace CDSP
{

/// Calculation of the absolute values of the Hilbert transform
void AbsHilbert(wxVector<SIGNALTYPE>& data)
{
	int i, sizeInput;
	alglib::complex_1d_array in;

	// create complex array
	sizeInput = data.size();
	in.setlength(sizeInput);

	for (i = 0; i < sizeInput; i++)
	{
		in[i].x = data.at(i);
		in[i].y = 0;
	}

	// FFT
	alglib::fftc1d(in);

	// H
	wxVector<int> h(data.size(), 0);
	if (2*floor(sizeInput/2) == sizeInput)
	{
		// even
		h.at(0) = 1;
		h.at(sizeInput/2) = 1;
		for (i = 1; i < sizeInput/2; i++)
			h.at(i) = 2;
	}
	else
	{
		// odd
		h.at(0) = 1;
		for (i = 1; i < (sizeInput+1)/2; i++)
			h.at(i) = 2;
	}

	// IFFT
	for (i = 0; i < sizeInput; i++)
	{
		in[i].x *= h[i];
		in[i].y *= h[i];
	}

	// IFFT
	alglib::fftc1dinv(in);

	// Absolute value
	for (i = 0; i < sizeInput; i++)
	{
		complex<SIGNALTYPE> tmp(in[i].x, in[i].y);
		data.at(i) = abs(tmp);
	}
}

// Code from here up to filtfilt was taken from: http://stackoverflow.com/a/27270420/287933
void add_index_range(vector<int> &indices, int beg, int end, int inc = 1)
{
	for (int i = beg; i <= end; i += inc)
	{
	   indices.push_back(i);
	}
}

void add_index_const(vector<int> &indices, int value, size_t numel)
{
	while (numel--)
	{
		indices.push_back(value);
	}
}

template<class T>
void append_vector(vector<T> &vec, const vector<T> &tail)
{
	vec.insert(vec.end(), tail.begin(), tail.end());
}

template<class T>
vector<T> subvector_reverse(const vector<T> &vec, int idx_end, int idx_start)
{
	vector<T> result(&vec[idx_start], &vec[idx_end+1]);
	std::reverse(result.begin(), result.end());
	return result;
}

inline int max_val(const vector<int>& vec)
{
	return std::max_element(vec.begin(), vec.end())[0];
}

template<class T>
void filter(vector<T> B, vector<T> A, const vector<T> &X, vector<T> &Y, vector<T> &Zi)
{
	if (A.empty())
	{
		throw std::domain_error("The feedback filter coefficients are empty.");
	}
	if (std::all_of(A.begin(), A.end(), [](double coef){ return coef == 0; }))
	{
		throw std::domain_error("At least one of the feedback filter coefficients has to be non-zero.");
	}
	if (A[0] == 0)
	{
		throw std::domain_error("First feedback coefficient has to be non-zero.");
	}

	// Normalize feedback coefficients if a[0] != 1;
	auto a0 = A[0];
	if (a0 != 1.0)
	{
		std::transform(A.begin(), A.end(), A.begin(), [a0](double v) { return v / a0; });
		std::transform(B.begin(), B.end(), B.begin(), [a0](double v) { return v / a0; });
	}

	size_t input_size = X.size();
	size_t filter_order = std::max(A.size(), B.size());
	B.resize(filter_order, 0);
	A.resize(filter_order, 0);
	Zi.resize(filter_order, 0);
	Y.resize(input_size);

	const T *x = &X[0];
	const T *b = &B[0];
	const T *a = &A[0];
	T *z = &Zi[0];
	T *y = &Y[0];

	for (size_t i = 0; i < input_size; ++i)
	{
		size_t order = filter_order - 1;
		while (order)
		{
			if (i >= order)
			{
				z[order - 1] = b[order] * x[i - order] - a[order] * y[i - order] + z[order];
			}
			--order;
		}
		y[i] = b[0] * x[i] + z[0];
	}
	Zi.resize(filter_order - 1);
}

// This is a function overload hack (because templates didn't work for some reason).
#define T float
#define M MatrixXf
#define V VectorXf
#include "filtfilt.h"
#undef T
#undef M
#undef V

#define T double
#define M MatrixXd
#define V VectorXd
#include "filtfilt.h"
#undef T
#undef M
#undef V

const double PId = M_PIl;

enum FILTERTYPE
{
	/// Flag of LOWPASS filter type
	LOWPASS = 0,
	/// Flag of HIGHPASS filter type
	HIGHPASS = 1
};

/**
 * Coefficients of filter design
 */
struct coeff {
	/// A numerator
	double* num;
	/// A denumerator
	double* den;
};

/// Data type specifications - FLOAT
#define RTF_FLOAT	0
/// Data type specifications - DOUBLE
#define RTF_DOUBLE	1

/// Butterworth filter order selection
void Buttord(const double& wp, double& ws, const double& rp, const double& rs, int& order)
{
	int ftype = 0;
	double WP, WS, WA, WN, W0;
	double ord;
	double x, y, f1, f2, log;

	//wxString str;

	if ( wp <= 0 || wp >= 1 || ws <= 0 || ws >= 1)
	{
		//str.Printf(wxT("Bad input parameters for Buttord! \n wp = %f, ws = %f"), wp, ws);
		//throw new CException(str, wxT("CDSP::Buttord"));
		throw runtime_error("Bad input parameters for Buttord!");
	}

	if (wp < ws)
		ftype = ftype + 1;  // low (1)
	else
		ftype = ftype + 2;  // high (2)

	WP = tan(PId*wp/2);
	WS = tan(PId*ws/2);

	if (ftype == 1)         // low (1)
		WA = WS/WP;
	else if (ftype == 2)    // high (2)
		WA = WP/WS;

	WA = abs(WA);

	x = 0.1 * abs(rs);
	y = 0.1 * abs(rp);
	f1 = pow(10, x) - 1;
	f2 = pow(10, y) - 1;
	log = log10(f1/f2);
	ord = ceil( log / (2*log10(WA)) );

	x = pow(10, 0.1*abs(rs)) - 1;
	y = 1/(2*(abs(ord)));
	W0 = WA/pow(x,y);

	if (ftype == 1)         // low
		WN = W0 * WP;
	else if (ftype == 2)    // high
		WN = WP / W0;

	order = ord;
	ws = ((2/PId)*atan(WN));
}

void getPoleCoefs(double p, double np, double fc, double r, int highpass, double a[3], double b[3])
{
	double rp, ip, es, vx, kx, t, w, m, d, x0, x1, x2, y1, y2, k;

	// calculate pole locate on the unit circle
	rp = -cos(PId / (np * 2.0) + (p - 1.0) * PId / np);
	ip = sin(PId / (np * 2.0) + (p - 1.0) * PId / np);

	// Warp from a circle to an ellipse
	if (r != 0.0) {
		es = sqrt(pow(1.0 / (1.0 - r), 2) - 1.0);
		vx = asinh(1.0/es) / np;
		kx = acosh(1.0/es) / np;
		kx = cosh( kx );
		rp *= sinh(vx) / kx;
		ip *= cosh(vx) / kx;
	}

	// s to z domains conversion
	t = 2.0*tan(0.5);
	w = 2.0*PId*fc;
	m = rp*rp + ip*ip;
	d = 4.0 - 4.0*rp*t + m*t*t;
	x0 = t*t/d;
	x1 = 2.0*t*t/d;
	x2 = t*t/d;
	y1 = (8.0 - 2.0*m*t*t)/d;
	y2 = (-4.0 - 4.0*rp*t - m*t*t)/d;

	// LP(s) to LP(z) or LP(s) to HP(z)
	if (highpass)
		k = -cos(w/2.0 + 0.5)/cos(w/2.0 - 0.5);
	else
		k = sin(0.5 - w/2.0)/sin(0.5 + w/2.0);
	d = 1.0 + y1*k - y2*k*k;
	a[0] = (x0 - x1*k + x2*k*k)/d;
	a[1] = (-2.0*x0*k + x1 + x1*k*k - 2.0*x2*k)/d;
	a[2] = (x0*k*k - x1*k + x2)/d;
	b[1] = (2.0*k + y1 + y1*k*k - 2.0*y2*k)/d;
	b[2] = (-k*k - y1*k + y2)/d;
	if (highpass) {
		a[1] *= -1.0;
		b[1] *= -1.0;
	}
}

int computeChebyIir(double *num, double *den, unsigned int num_pole,
			  int highpass, double ripple, double cutoff_freq)
{
	double *a, *b, *ta, *tb;
	double ap[3], bp[3];
	double sa, sb, gain;
	unsigned int i, p;
	int retval = 1;

	// Allocate temporary arrays
	a = new double[num_pole + 3];
	b = new double[num_pole + 3];
	ta = new double[num_pole + 3];
	tb = new double[num_pole + 3];

	if (!a || !b || !ta || !tb) {
		retval = 0;
		goto exit;
	}
	memset(a, 0, (num_pole + 3) * sizeof(*a));
	memset(b, 0, (num_pole + 3) * sizeof(*b));

	a[2] = 1.0;
	b[2] = 1.0;

	for (p = 1; p <= num_pole / 2; p++) {
		// Compute the coefficients for this pole
		getPoleCoefs(p, num_pole, cutoff_freq, ripple, highpass, ap, bp);

		// Add coefficients to the cascade
		memcpy(ta, a, (num_pole + 3) * sizeof(*a));
		memcpy(tb, b, (num_pole + 3) * sizeof(*b));
		for (i = 2; i <= num_pole + 2; i++) {
			a[i] = ap[0]*ta[i] + ap[1]*ta[i-1] + ap[2]*ta[i-2];
			b[i] = tb[i] - bp[1]*tb[i-1] - bp[2]*tb[i-2];
		}
	}

	// Finish combining coefficients
	b[2] = 0.0;
	for (i = 0; i <= num_pole; i++) {
		a[i] = a[i + 2];
		b[i] = -b[i + 2];
	}

	// Normalize the gain
	sa = sb = 0.0;
	for (i = 0; i <= num_pole; i++) {
		sa += a[i] * ((highpass && i % 2) ? -1.0 : 1.0);
		sb += b[i] * ((highpass && i % 2) ? -1.0 : 1.0);
	}
	gain = sa / (1.0 - sb);
	for (i = 0; i <= num_pole; i++)
		a[i] /= gain;

	// Copy the results to the num and den
	for (i = 0; i <= num_pole; i++) {
		num[i] = a[i];
		den[i] = -b[i];
	}
	// den[0] must be 1.0
	den[0] = 1.0;

exit:
	delete [] a;
	delete [] b;
	delete [] ta;
	delete [] tb;

	return retval;
}

/// Calculate coefficients for Butterworth filter.
int calcButterCoeff(unsigned int nchann, int proctype, double fc,
					  unsigned int num_pole, int highpass, struct coeff *coeff)
{
	double* num = NULL,* den = NULL;
	double ripple = 0.0;
	int res = 0;

	if (num_pole % 2 != 0)
		return -1;

	num = new double[num_pole+1];
	if (num == NULL)
		return -2;
	den = new double[num_pole+1];
	if (den == NULL) {
		res = -3;
		goto err1;
	}

	/* Prepare the z-transform of the filter */
	if (!computeChebyIir(num, den, num_pole, highpass, ripple, fc)) {
		res = -4;
		goto err2;
	}

	coeff->num = num;
	coeff->den = den;

	return 0;

err2:
	delete [] den;
err1:
	delete [] num;
	return res;
}

/// Designs an Nth order lowpass digital Butterworth filter and returns the filter coefficients in length N+1 vectors B (numerator) and A (denominator)
bool Butter(wxVector<double>& b, wxVector<double>& a, const int& order, const double& Wn, FILTERTYPE ftype)
{
	struct   coeff coeff;
	int      res;
	unsigned i;

	unsigned int nchann     = 1;                /* channels number */
	int proctype            = RTF_DOUBLE;       /* samples have float type */
	double fc               = Wn/2;             /* normalized cut-off frequency, Hz */
	unsigned int num_pole   = order;            /* filter order */
	int highpass            = ftype;            /* lowpass filter */

	res = calcButterCoeff(nchann, proctype, fc, num_pole, highpass, &coeff);
	if (res != 0) {
		//cout << "Error: unable to calculate coefficients: " << res << endl;
		//return false;
		throw runtime_error("Error: unable to calculate coefficients");
	}

	b.clear();
	a.clear();

	for (i = 0; i < num_pole+1; i++)
	{
		a.push_back(coeff.den[i]);
		b.push_back(coeff.num[i]);
	}

	delete [] coeff.num;
	delete [] coeff.den;
	return true;
}

/// Digital signal filtering 10-60Hz
void Filtering(wxVector<SIGNALTYPE>* data, const int& countChannels, const int& fs, const BANDWIDTH& bandwidth)
{
	int 			 i;
	wxVector<double> a, b;
	//CFiltFilt** 	 threads = new CFiltFilt*[countChannels];

	double           Wp, Ws, Rp = 6, Rs = 60;
	int              order;

	//CException*      e = NULL;

	// Precomputation values - Chebyshev II, for band low = 10, band high = 60
	// high
	const double bh[7] = {0.609690631014657, -3.63559011319046, 9.05535356429358, -12.0589078771451, 9.05535356429359, -3.63559011319046, 0.609690631014657};
	const double ah[7] = {1, -4.98146241953972, 10.4387837885612, -11.7670343111137, 7.51972533817281, -2.58144797118803, 0.371722665567175};

	// low
	const double bl[9] = {0.0756582074206937, 0.473325643095456, 1.40026011018072, 2.53803955902203, 3.07132128851322, 2.53803955902203, 1.40026011018072, 0.473325643095457, 0.0756582074206939};
	const double al[9] = {1, 2.04764667081498, 3.02818564408939, 2.79774379822106, 1.90643217355482, 0.903943414012386, 0.296359133636092, 0.0598505403817098, 0.00572695324055361};


	// HIGH pass filtering
	if (bandwidth.m_bandLow != 10 || bandwidth.m_bandHigh != 60)
	{
		Wp = 2*bandwidth.m_bandLow/(double)fs;
		Ws = (2*bandwidth.m_bandLow/(double)fs)-0.05;
		if (Ws < 0) Ws = 0.05;

		// calc filter order
		Buttord(Wp, Ws, Rp, Rs, order);

		// if order is odd + 1
		if (order % 2 != 0)
			order ++;

		// filter design
		if (!Butter(b, a, order, Ws, HIGHPASS))
		{
			//throw new CException(wxT("Error calculating filter design for HIGHPASS filter!"), wxT("CDSP::Filtering"));
			throw runtime_error("Error calculating filter design for HIGHPASS filter!");
		}
	}
	else
	{
		a.assign(ah, ah+7);
		b.assign(bh, bh+7);
	}

	// run filtering
	vector<float> out, af, bf;

	for (double e : a)
		af.push_back(e);

	for (double e : b)
		bf.push_back(e);

	for (i = 0; i < countChannels; i++)
	{
		filtfilt(bf, af, data[i], out);

		assert(data[i].size() == out.size());
		data[i] = out;

		//threads[i] = new CFiltFilt(b, a, &data[i]);
		//threads[i]->Run();
	}
	// wait until they work
	/*for (i = 0; i < countChannels; i++)
	{
		threads[i]->Wait();
		delete threads[i];
	}*/

	if (bandwidth.m_bandHigh == fs/2)
	{
		//delete [] threads;
		return;
	}

	// LOW pass filtering
	if (bandwidth.m_bandLow != 10 || bandwidth.m_bandHigh != 60)
	{
		Wp = 2*bandwidth.m_bandHigh/(double)fs;
		Ws = (2*bandwidth.m_bandHigh/(double)fs)+0.1;
		if (Ws > 1) Ws = 1;

		// calc filter order
		Buttord(Wp, Ws, Rp, Rs, order);

		// if order is odd + 1
		if (order % 2 != 0)
			order ++;

		// filter design
		if (!Butter(b, a, order, Ws, LOWPASS))
		{
			//throw new CException(wxT("Error calculating filter design for HIGHPASS filter!"), wxT("CDSP::Filtering"));
			throw runtime_error("Error calculating filter design for HIGHPASS filter!");
		}
	}
	else
	{
		a.assign(al, al+9);
		b.assign(bl, bl+9);
	}

	af.clear();
	for (double e : a)
		af.push_back(e);

	bf.clear();
	for (double e : b)
		bf.push_back(e);

	for (i = 0; i < countChannels; i++)
	{
		filtfilt(bf, af, data[i], out);

		assert(data[i].size() == out.size());
		data[i] = out;

		//threads[i] = new CFiltFilt(b, a, &data[i]);
		//threads[i]->Run();
	}
	// wait until they work
	/*for (i = 0; i < countChannels; i++)
	{
		e = (CException*)threads[i]->Wait();
		delete threads[i];

		if (e != NULL)
			throw e;
	}*/

	//delete [] threads;
}

/// Digital signal filtering Nx50hz
void Filt50Hz(wxVector<SIGNALTYPE>* data, const int& countChannels, const int& fs, const int& hum_fs, const BANDWIDTH& bandwidth)
{
	double 			 R = 1, r = 0.985, tmp, i;
	wxVector<int> 	 f0;
	wxVector<float> b, a;

	for (i = hum_fs; i <= fs/2 && i <= bandwidth.m_bandHigh; i+= hum_fs)
		f0.push_back(i);

	for (i = 0; i < f0.size(); i++)
	{
		b.clear();
		a.clear();

		b.push_back(1);
		tmp = -2 * R * cos( 2*M_PI*f0[i] / fs);
		b.push_back(tmp);
		b.push_back(R*R);

		a.push_back(1);
		tmp = -2 * r * cos( 2*M_PI*f0[i] / fs);
		a.push_back(tmp);
		a.push_back(r*r);

		vector<float> out;

		//filt50Hz(data, countChannels, b, a);
		for (int j = 0; j < countChannels; j++)
		{
			filtfilt(b, a, data[j], out);

			assert(data[j].size() == out.size());
			data[j] = out;
		}
	}
}

/*/// Digital signal filtering Nx50hz
void filt50Hz(wxVector<SIGNALTYPE>* data, const int& countChannels, const wxVector<double>& B, const wxVector<double>& A)
{
	int i;
	CFiltFilt** threads = new CFiltFilt*[countChannels];

	for (i = 0; i < countChannels; i++)
	{
		threads[i] = new CFiltFilt(B, A, &data[i]);
		threads[i]->Run();
	}

	for (i = 0; i < countChannels; i++)
	{
		threads[i]->Wait();
		delete threads[i];
	}

	delete [] threads;
}*/

} // namespace CDSP

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
		// These thwo vectors were allocated in localMaximaDetection(). Sometimes they point to the same vector.
		delete m_markersHigh;
		if (m_markersHigh != m_markersLow)
			delete m_markersLow;
	}

} ONECHANNELDETECTRET;

/**
 * Implementation "one_channel_detect" function from reference implementation of spike detector.
 */
class COneChannelDetect
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
	virtual ONECHANNELDETECTRET* Entry();

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
// COneChannelDetect
// ------------------------------------------------------------------------------------------------

/// A constructor
COneChannelDetect::COneChannelDetect(const wxVector<SIGNALTYPE>* data, const DETECTOR_SETTINGS* settings, const int& fs, const wxVector<int>* index,
									 const int& channel)
	: m_data(data), m_settings(settings), m_fs(fs), m_index(index), m_channel(channel)
{
	/* empty */
}

/// A destructor
COneChannelDetect::~COneChannelDetect()
{
	/* empty */
}

/// This is the entry point of the thread.
ONECHANNELDETECTRET* COneChannelDetect::Entry()
{

	wxVector<SIGNALTYPE>  envelope(m_data->begin(), m_data->end());
	int 				  start, stop, tmp, i, j;
	int 				  indexSize = m_index->size();
	int     			  envelopeSize = envelope.size();
	double 				  std, l, m;

	wxVector<double>      logs;
	alglib::real_1d_array r1a_logs;

	//wxVector<SIGNALTYPE>  phatMedian, phatStd;
	vector<double>  phatMedian, phatStd;

	ONECHANNELDETECTRET*  ret = NULL;

	// Hilbert's envelope (intense envelope)
	CDSP::AbsHilbert(envelope);
	//OpenCLContext::printBuffer("envelope.txt", envelope.data(), envelope.size());

	for (i = 0; i < indexSize; i++)
	{
		start = m_index->at(i);
		stop = start + m_settings->m_winsize*m_fs - 1;

		for (j = start; j <= stop; j++)
		{
			if (envelope[j] > 0)
			{
				l = log(envelope[j]);
				logs.push_back(l);
			}
		}

		m = mean(logs);
		std = sqrt(variance(logs, m));

		phatMedian.push_back(m);
		phatStd.push_back(std);
		logs.clear();
	}

	double r = (double)envelope.size() / (double)indexSize;
	double n_average = m_settings->m_winsize;

	wxVector<double> b, a;
	tmp = round(n_average*m_fs/r);

	if (tmp > 1)
	{
		a.push_back(1);
		for (i = 0; i < tmp; i++)
		{
			l = 1.0/tmp;
			b.push_back(l);
		}

		// I don't know why this is here, so better keep these two filters.
		vector<double> out;
		CDSP::filtfilt(b, a, phatMedian, out);
		//OpenCLContext::printBufferDouble("phat.txt", out.data(), out.size());
		//OpenCLContext::printBufferDouble("old_phat.txt", phatMedian.data(), phatMedian.size());
		phatMedian = out;

		out.clear();
		CDSP::filtfilt(b, a, phatStd, out);
		phatStd = out;
	}

	// interpolation of thresholds value to threshold curve (like backround)
	wxVector<double> phat_int[2];
	wxVector<double> x, y;
	alglib::real_1d_array xreal, yrealMedian, yrealStd, x2real, retMedian, retStd;

	if (phatMedian.size() > 1)
	{
		for (i = 0; i < indexSize; i++)
			x.push_back( m_index->at(i) + round( (m_settings->m_winsize * m_fs)/2 ) );

		start = *m_index->begin();
		stop = m_index->back();
		y.reserve(m_index->back());
		for (i = start; i <= stop; i++)
			y.push_back(i + round( (m_settings->m_winsize * m_fs)/2 ));

		try
		{
			// interpolation Median and Std
			wxVector<double> phatMedVecDoub(phatMedian.begin(), phatMedian.end());
			wxVector<double> phatStdVecDoub(phatStd.begin(), phatStd.end());

			xreal.setcontent(x.size(), &x[0]);
			yrealMedian.setcontent(phatMedVecDoub.size(), &phatMedVecDoub[0]);
			yrealStd.setcontent(phatStdVecDoub.size(), &phatStdVecDoub[0]);
			x2real.setcontent(y.size(), &y[0]);

			alglib::spline1dconvcubic(xreal, yrealMedian, x2real, retMedian);
			alglib::spline1dconvcubic(xreal, yrealStd, x2real, retStd);
		}
		catch(alglib::ap_error e)
		{
			//cout << "Aglib msg: " << e.msg.c_str() << endl;
			return NULL;
		}

		for (i = 0; i < retMedian.length(); i++)
		{
			phat_int[0].push_back(retMedian[i]);
			phat_int[1].push_back(retStd[i]);
		}

		// DOPLNENI
		double temp_elem0 = phat_int[0].front();
		double temp_elem1 = phat_int[1].front();

//		for (i = 0; i < floor(m_settings->m_winsize * m_fs / 2); i++)
//		{
//			phat_int[0].insert(phat_int[0].begin(), temp_elem0);
//			phat_int[1].insert(phat_int[1].begin(), temp_elem1);
//		}
		// My optimization: Insert all elements at once to prevent repeated copying in the vector. Makes is about 2.5x faster for an 1 kHz sample with no decimation.
		int n = floor(m_settings->m_winsize * m_fs / 2);
		phat_int[0].insert(phat_int[0].begin(), n, temp_elem0);
		phat_int[1].insert(phat_int[1].begin(), n, temp_elem1);

		temp_elem0 = phat_int[0].back();
		temp_elem1 = phat_int[1].back();
		for (i = phat_int[0].size(); i < envelopeSize; i++)
		{
			phat_int[0].push_back(temp_elem0);
			phat_int[1].push_back(temp_elem1);
		}
	}
	else
	{
		for (i = 0; i < envelopeSize; i++)
		{
			phat_int[0].push_back( phatMedian.at(1) * 1 );
			phat_int[1].push_back( phatStd.at(1) * 1 );
		}
	}

	// LOGNORMAL distr.
	double tmp_diff, tmp_exp, tmp_square, lognormal_mode, lognormal_median, tmp_prah_int, lognormal_mean, tmp_sum;
	wxVector<double> prah_int[2];

	double tmp_sqrt_one, tmp_to_erf, tmp_log, tmp_erf, tmp_pdf, tmp_x, tmp_x2;
	double tmp_sqrt = sqrt(2*M_PI);
	wxVector<double> envelope_cdf, envelope_pdf;
	int phatIntSize = phat_int[0].size();

	for (i = 0; i < phatIntSize; i++)
	{
		tmp_square = phat_int[1].at(i) * phat_int[1].at(i);
		tmp_diff   = phat_int[0].at(i) - tmp_square;
		tmp_sum    = phat_int[0].at(i) + tmp_square/2;

		lognormal_mode   = exp(tmp_diff);
		lognormal_median = exp(phat_int[0].at(i));
		lognormal_mean   = exp(tmp_sum);

		tmp_prah_int = (m_settings->m_k1 * (lognormal_mode + lognormal_median)) - (m_settings->m_k3 * (lognormal_mean-lognormal_mode));
		prah_int[0].push_back(tmp_prah_int);

		if (m_settings->m_k2 != m_settings->m_k1)
		{
			tmp_prah_int = (m_settings->m_k2 * (lognormal_mode + lognormal_median)) - (m_settings->m_k3 * (lognormal_mean-lognormal_mode));
			prah_int[1].push_back(tmp_prah_int);
		}

		// CDF of lognormal distribution, PDF of lognormal distribution
		// CDF
		tmp_sqrt_one = sqrt( 2.0f * phat_int[1].at(i) * phat_int[1].at(i));
		tmp_log = log(envelope[i]);
		tmp_to_erf = (tmp_log - phat_int[0].at(i)) / tmp_sqrt_one;
		tmp_erf = erf(tmp_to_erf);
		envelope_cdf.push_back(0.5 + 0.5 * tmp_erf);

		// PDF
		tmp_x = (tmp_log - phat_int[0].at(i)) / phat_int[1].at(i);
		tmp_x *= tmp_x;
		tmp_exp = exp( -0.5 * tmp_x );
		tmp_x2 = (envelope[i] * phat_int[1].at(i) * tmp_sqrt);
		tmp_pdf = tmp_exp / tmp_x2;
		envelope_pdf.push_back(tmp_pdf);
	}

	wxVector<bool>* markers_high = NULL,* markers_low = NULL;
	try {
		markers_high = localMaximaDetection(envelope, prah_int[0], m_settings->m_polyspike_union_time);
		detectionUnion(markers_high, envelope, m_settings->m_polyspike_union_time * m_fs);

		if (m_settings->m_k2 != m_settings->m_k1 && prah_int[1].size() != 0)
		{
			markers_low = localMaximaDetection(envelope, prah_int[1], m_settings->m_polyspike_union_time);
			detectionUnion(markers_low, envelope, m_settings->m_polyspike_union_time * m_fs);
		} else markers_low = markers_high;
	}
	//catch (CException* e)
	catch (exception e)
	{
		fprintf(stderr, "Exception: %s", e.what());
		return NULL;
	}

	ret = new ONECHANNELDETECTRET(markers_high, markers_low, prah_int, envelope_cdf, envelope_pdf, envelope);
	return ret;
}

/// Calculating a mean from data in vector
double COneChannelDetect::mean(wxVector<double>& data)
{
	float sum = 0;
	wxVector<double>::iterator b = data.begin();
	wxVector<double>::iterator e = data.end();

	while (b != e)
	{
		sum = sum + *b;
		++b;
	}

	return sum / data.size();
}

/// Calculating a variance from data in vector
double COneChannelDetect::variance(wxVector<double>& data, const double & mean)
{
	wxVector<double> v(data.begin(), data.end());

	wxVector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(),
			std::bind2nd(std::minus<double>(), mean));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

	return sq_sum / (v.size()-1);
}

/// Detection of local maxima in envelope
wxVector<bool>* COneChannelDetect::localMaximaDetection(wxVector<SIGNALTYPE>& envelope, const wxVector<double>& prah_int, const double& polyspike_union_time)
{
	unsigned int         size = envelope.size();
	wxVector<bool>* 	 marker1 = new wxVector<bool>(size, 0); // This leak was fixed in ~oneChannelDetectRet().
	wxVector<int>   	 point[2];
	wxVector<SIGNALTYPE> seg, seg_s;
	wxVector<int>        tmp_diff_vector;
	int        			 pointer_max;
	SIGNALTYPE 			 tmp_max;
	int 				 tmp_pointer;
	unsigned int 		 i, j, tmp;

	wxVector<int> 		 pointer;
	bool 				 state_previous = false;
	int 				 tmp_ceil, tmp_stop, start;
	float 				 tmp_sum_elems;

	wxVector<int> 		 lokal_max, lokal_max_poz;
	wxVector<SIGNALTYPE> lokal_max_val;
	float 				 tmp_diff;
	int 				 tmp_sign, tmp_sign_next;
	wxVector<float> 	 tmp_diff_vector2;

	for (i = 0; i < size; i++)
		if (envelope[i] > prah_int[i])
			marker1->at(i) = 1;

	// start + end crossing
	findStartEndCrossing(point, marker1);
	if (point[0].size() != point[1].size())
	{
		//throw new CException(wxT("local_maxima_detection: point sizes are different"), wxT("COneChannelDetect::localMaximaDetection"));
		throw runtime_error("local_maxima_detection: point sizes are different");
	}

	marker1->assign(size, 0);

	for (i = 0; i < point[0].size(); i++)
	{
		seg_s.clear();
		seg.clear();
		tmp_diff_vector.clear();

		// detection of local maxima in section which crossed threshold curve
		if (point[1].at(i) - point[0].at(i) > 2)
		{
			seg.assign(&envelope[0] + point[0].at(i), &envelope[0] + point[1].at(i)+1);
			for (j = 0; j < seg.size()-1; j++)
				seg_s.push_back(seg.at(j+1) - seg.at(j));

			for (j = 0; j < seg_s.size(); j++)
			{
				if (seg_s.at(j) > 0) seg_s.at(j) = 1;
				else if(seg_s.at(j) < 0) seg_s.at(j) = -1;
				else seg_s.at(j) = 0;
			}

			// positions of local maxima in the section
			seg_s.insert(seg_s.begin(), 0);
			for (j = 0; j < seg_s.size()-1; j++)
				tmp_diff_vector.push_back(seg_s.at(j+1) - seg_s.at(j));

			seg_s.clear();
			for (j = 0; j < tmp_diff_vector.size(); j++)
				if (tmp_diff_vector.at(j) < 0)
					seg_s.push_back(j);

			for (j = 0; j < seg_s.size(); j++)
			{
				tmp_pointer = point[0].at(i) + seg_s.at(j);
				if (tmp_pointer < (int)size)
					marker1->at(tmp_pointer) = 1;
			}
		}
		else if (point[1].at(i) - point[0].at(i) <= 2)
		{
			pointer_max = 1;
			tmp_max = 0;
			seg.assign(&envelope[0] + point[0].at(i), &envelope[0] + point[1].at(i)+1);
			for (j = 0; j < seg.size(); j++)
			{
				if (seg.at(j) > tmp_max)
				{
					pointer_max = j;
					tmp_max = seg.at(j);
				}
			}

			tmp = point[0].at(i) + pointer_max;
			if (tmp < size)
				marker1->at(tmp) = 1;
		}
	}

	// union of section, where local maxima are close together <(1/f_low + 0.02 sec.)~ 120 ms
	for (i = 0; i < marker1->size(); i++)
		if (marker1->at(i))
			pointer.push_back(i);

	state_previous = false;
	for (i = 0; i < pointer.size(); i++)
	{
		tmp_sum_elems = 0;
		seg.clear();

		tmp_ceil = ceil(pointer.at(i) + polyspike_union_time * m_fs);
		if (tmp_ceil >= (int)size)
			tmp_stop = size;
		else
			tmp_stop = tmp_ceil + 1;

		for(j = pointer.at(i)+1; j < (unsigned)tmp_stop; j++)
				seg.push_back(marker1->at(j));

		for (wxVector<SIGNALTYPE>::iterator j = seg.begin() ; j != seg.end(); ++j)
			tmp_sum_elems += *j;

		if (state_previous)
		{
			if (tmp_sum_elems > 0)
				state_previous = true;
			else
			{
				state_previous = false;
				for (j = start; j <= (unsigned)pointer.at(i) && j < size; j++)
					marker1->at(j) = true;
			}
		}
		else
		{
			if (tmp_sum_elems > 0)
			{
				state_previous = true;
				start = pointer.at(i);
			}
		}
	}

	// finding of the highes maxima of the section with local maxima
	findStartEndCrossing(point, marker1);

	if (point[0].size() != point[1].size())
	{
		//throw new CException(wxT("local_maxima_detection: point sizes are different 2"), wxT("COneChannelDetect::localMaximaDetection"));
		throw runtime_error("local_maxima_detection: point sizes are different 2");
	}

	// local maxima with gradient in souroundings
	for (i = 0; i < point[0].size(); i++)
	{
		lokal_max.clear();
		lokal_max_val.clear();
		lokal_max_poz.clear();

		tmp_diff = point[1].at(i) - point[0].at(i);
		if (tmp_diff > 1)
		{
			for (j = 0; j < pointer.size(); j++)
				if (pointer.at(j) >= point[0][i] && pointer.at(j) <= point[1][i])
					lokal_max.push_back(pointer.at(j));

			for (j = point[0].at(i); j <= (unsigned)point[1].at(i) && j < size; j++)
				marker1->at(j) = false;

			// envelope magnitude in local maxima
			for (j = 0; j < lokal_max.size(); j++)
				lokal_max_val.push_back(envelope[lokal_max.at(j)]);

			// lokal_max_poz=(diff(sign(diff([0;lokal_max_val;0]))<0)>0);
			// diff([0;lokal_max_val;0])
			tmp_diff_vector2.clear();
			tmp_diff = 0;
			for (j = 0; j < lokal_max_val.size(); j++)
			{
				tmp_diff_vector2.push_back( lokal_max_val.at(j) - tmp_diff );
				tmp_diff = lokal_max_val.at(j);
			}
			tmp_diff_vector2.push_back( 0 - lokal_max_val.back() );

			// sign(diff([0;lokal_max_val;0]))<0
			lokal_max_poz.clear();
			for (j = 0; j < tmp_diff_vector2.size()-1; j++)
			{
				if (tmp_diff_vector2.at(j) > 0) tmp_sign = 1;
				else tmp_sign = 0;

				if (tmp_diff_vector2.at(j+1) > 0) tmp_sign_next = 1;
				else tmp_sign_next = 0;

				tmp_diff = tmp_sign - tmp_sign_next;
				lokal_max_poz.push_back(tmp_diff);
			}

			for (j = 0; j < lokal_max.size(); j++)
			{
				if (lokal_max_poz.at(j) == 1 && lokal_max.at(j) < (int)size)
				{
					marker1->at( lokal_max.at(j) ) = true;
				}
			}
		}
	}

	return marker1;
}

/// Detecting of union and their merging.
void COneChannelDetect::detectionUnion(wxVector<bool>* marker1, wxVector<SIGNALTYPE>& envelope, const double& union_samples)
{
	int 				  i, j, start, stop, sum = round(union_samples);
	float 			 	  max; // maximum value in segment of envelope
	int 			 	  max_pos; // position of maximum
	wxVector<double> 	  vec_MASK(sum, 1.0);
	wxVector<double> 	  vec_marker2(marker1->begin(), marker1->end());
	wxVector<int>    	  point[2];
	alglib::real_1d_array r1a_marker2, r1a_MASK, r1a_ret;

	// dilatation
	r1a_marker2.setcontent(vec_marker2.size(), &vec_marker2[0]);
	r1a_MASK.setcontent(vec_MASK.size(), &vec_MASK[0]);
	vec_marker2.clear();
	marker1->assign(marker1->size(), false);
	alglib::convr1d(r1a_marker2, r1a_marker2.length(), r1a_MASK, r1a_MASK.length(), r1a_ret);

	start = (vec_MASK.size() / 2);
	stop  = r1a_ret.length() - (vec_MASK.size() - start);

	for (i = start; i <= stop; i++)
		if (r1a_ret[i] > 0)
			marker1->at(i-start) = true;

	// erosion
	vec_marker2.assign(marker1->begin(), marker1->end());
	marker1->assign(marker1->size(), false);
	r1a_marker2.setcontent(vec_marker2.size(), &vec_marker2[0]);
	alglib::convr1d(r1a_marker2, r1a_marker2.length(), r1a_MASK, r1a_MASK.length(), r1a_ret);

	stop  = r1a_ret.length() - (vec_MASK.size() - start);
	for (i = start; i <= stop; i++)
		if (r1a_ret[i] == sum)
			marker1->at(i-start) = true;

	// start + end crossing
	findStartEndCrossing(point, marker1);

	marker1->assign(marker1->size(), false);
	for (i = 0; i < (int)point[0].size(); i++)
	{
		max = -1;
		for (j = point[0].at(i); j <= point[1].at(i); j++)
		{
			if (envelope[j] > max)
			{
				max = envelope[j];
				max_pos = j;
			}
		}

		marker1->at(max_pos) = true;
	}
}

// Finding of the highes maxima of the section with local maxima
void COneChannelDetect::findStartEndCrossing(wxVector<int> point[2], const wxVector<bool>* marker1)
{
	bool one = false;
	bool tmp_marker;
	int  size = marker1->size();
	int i;

	point[0].clear();
	point[1].clear();

	for (i = 0; i < size; i++)
	{
		tmp_marker = marker1->at(i);
		if (!one) // start crossing
		{
			if (tmp_marker == 1)
			{
				one = true;
				point[0].push_back(i);

				if (i == size-1)
					point[1].push_back(i);
			}
		} else { // end crossing
			if (tmp_marker == 0 || i == size-1)
			{
				one = false;
				point[1].push_back(i-1);
			}
		}
	}
}

const int BLOCK_SIZE = 1024*16/*512*/;

int nearestGreaterDivisor(int a, int b)
{
	for (int i = a; a < b; i++)
	{
		if (b%i == 0)
			return i;
	}
	return b;
}

} // namespace

namespace AlenkaSignal
{

/// A constructor.
CDetectorOutput::CDetectorOutput()
{
	/* empty */
}

/// A virtual destructor.
CDetectorOutput::~CDetectorOutput()
{
	/* empty */
}


/// Add data to the vectors.
void CDetectorOutput::Add(const double& pos, const double& dur, const int& chan, const double& con, const double& weight, const double& pdf)
{
	m_pos.push_back(pos);
	m_dur.push_back(dur);
	m_chan.push_back(chan);
	m_con.push_back(con);
	m_weight.push_back(weight);
	m_pdf.push_back(pdf);
}

///Erase records at positions.
void CDetectorOutput::Remove(const wxVector<int>& pos)
{
	unsigned i, counter = 0;

	for (i = 0; i < pos.size(); i++)
	{
		if (m_pos.size() < pos.at(i)-counter || pos.at(i)-counter < 0)
			continue;

		m_pos.erase(m_pos.begin()+pos.at(i)-counter);
		m_dur.erase(m_dur.begin()+pos.at(i)-counter);
		m_chan.erase(m_chan.begin()+pos.at(i)-counter);
		m_con.erase(m_con.begin()+pos.at(i)-counter);
		m_weight.erase(m_weight.begin()+pos.at(i)-counter);
		m_pdf.erase(m_pdf.begin()+pos.at(i)-counter);

		counter++;
	}
}

// ------------------------------------------------------------------------------------------------
// CDischarges
// ------------------------------------------------------------------------------------------------

/// A constructor.
CDischarges::CDischarges(const int& countChannels)
{
	m_countChannels = countChannels;

	m_MV   = new std::vector<double>[countChannels];
	m_MA   = new std::vector<double>[countChannels];
	m_MP   = new std::vector<double>[countChannels];
	m_MD   = new std::vector<double>[countChannels];
	m_MW   = new std::vector<double>[countChannels];
	m_MPDF = new std::vector<double>[countChannels];
}

/// A virual destructor.
CDischarges::~CDischarges()
{
	delete [] m_MV;
	delete [] m_MA;
	delete [] m_MP;
	delete [] m_MD;
	delete [] m_MW;
	delete [] m_MPDF;
}

/**
 * Erase records.
 * @param pos positions of records.
 */
void CDischarges::Remove(const wxVector<int>& pos)
{
	unsigned i, channel, counter = 0;

	for (i = 0; i < pos.size(); i++)
	{

		for (channel = 0; channel < m_countChannels; channel++)
		{
			if (m_MV[channel].size() < pos.at(i)-counter || pos.at(i)-counter < 0)
				continue;

			m_MV[channel].erase(m_MV[channel].begin()+pos.at(i)-counter);
			m_MA[channel].erase(m_MA[channel].begin()+pos.at(i)-counter);
			m_MP[channel].erase(m_MP[channel].begin()+pos.at(i)-counter);
			m_MW[channel].erase(m_MW[channel].begin()+pos.at(i)-counter);
			m_MPDF[channel].erase(m_MPDF[channel].begin()+pos.at(i)-counter);
			m_MD[channel].erase(m_MD[channel].begin()+pos.at(i)-counter);
		}
		counter++;
	}

}

template<class T>
Spikedet<T>::Spikedet(int fs, int channelCount, DETECTOR_SETTINGS settings, OpenCLContext* context) :
	fs(fs), channelCount(channelCount), settings(settings), context(context)
{
	decimationF = settings.m_decimation;
	decimationF = nearestGreaterDivisor(min(fs, decimationF), fs);

	if (decimationF >= fs)
		return;

	int M = fs + 1;
	Filter<T> filter(M, fs);
	//filter.notch(true);
	//filter.highpass(true);
	filter.lowpass(true);

	//filter.setNotch(settings.m_main_hum_freq);
	//filter.setHighpass(settings.m_band_low);
	//filter.setLowpass(min(decimationF, settings.m_band_high));
	filter.setLowpass(decimationF);

	filterProcessor = new FilterProcessor<T>(BLOCK_SIZE, channelCount, context, WindowFunction::Hamming);
	filterProcessor->changeSampleFilter(M, filter.computeSamples());

#ifndef NDEBUG
	//FILE* file = fopen("spikedet_coefficients.txt", "w");
	//filter.printCoefficients(file, filterProcessor->getCoefficients());
	//fclose(file);
#endif

	cl_int err;
	cl_mem_flags flags = CL_MEM_READ_WRITE;

	queue = clCreateCommandQueue(context->getCLContext(), context->getCLDevice(), 0, &err);
	checkClErrorCode(err, "clCreateCommandQueue");

	inBuffer = clCreateBuffer(context->getCLContext(), flags, (BLOCK_SIZE + 4)*channelCount*sizeof(T), nullptr, &err);
	checkClErrorCode(err, "clCreateBuffer");

	outBuffer = clCreateBuffer(context->getCLContext(), flags, (BLOCK_SIZE + 4)*channelCount*sizeof(T), nullptr, &err);
	checkClErrorCode(err, "clCreateBuffer");
}

template<class T>
Spikedet<T>::~Spikedet()
{
	delete filterProcessor;

	cl_int err;

	if (queue != nullptr)
	{
		err = clReleaseCommandQueue(queue);
		checkClErrorCode(err, "clReleaseCommandQueue()");
	}

	if (inBuffer != nullptr)
	{
		err = clReleaseMemObject(inBuffer);
		checkClErrorCode(err, "clReleaseMemObject()");
	}
	if (outBuffer != nullptr)
	{
		err = clReleaseMemObject(outBuffer);
		checkClErrorCode(err, "clReleaseMemObject()");
	}
}

template<class T>
void Spikedet<T>::runAnalysis(SpikedetDataLoader<T>* loader, CDetectorOutput*& out, CDischarges*& discharges)
{
	m_out = out;
	m_discharges = discharges;
	int 				  i, j, k, indexSize;
	int64_t				  start, stop, tmp;
	//wxVector<SIGNALTYPE>* segments = NULL;
	wxVector<int64_t>         indexStart, indexStop; // TODO: Make sure all variables holding the sample index have 8 bytes like these ones.

	BANDWIDTH 			  bandwidth(m_settings->m_band_low, m_settings->m_band_high);

	CDetectorOutput*      subOut 		= NULL;
	CDischarges*          subDischarges = NULL;

	int 				  posSize, disSize, tmpFirst, tmpLast;
	double 				  minMP, tmpShift;

	wxVector<int>         removeOut;
	wxVector<int>         removeDish;

//	int 				  countSamples = m_model->GetCountSamples();
//	int 				  countChannels = m_model->GetCountChannels();
//	int 				  fs = m_model->GetFS();
	int64_t				  countSamples = loader->sampleCount();
	int 				  countChannels = loader->channelCount();

	int    	  			  winsize  = m_settings->m_winsize * fs;

	//wxThreadEvent         event(wxEVT_THREAD, DETECTOR_EVENT);
	int 				  progressBy;

	// verify buffering
	tmp = countSamples / fs;
	if (m_settings->m_buffering > tmp)
		m_settings->m_buffering = tmp;

	// Signal buffering
	int64_t N_seg = floor(countSamples/(m_settings->m_buffering * fs));
	if (N_seg < 1)
		N_seg = 1;
	int64_t T_seg = round((double)countSamples/(double)N_seg/fs);
	// Indexs of segments with two-side overlap
	getIndexStartStop(indexStart, indexStop, countSamples, T_seg, fs, winsize);

	progressBy = round(100/(float)indexStart.size());

	// starting analysis on the segmented data
	indexSize = indexStop.size();
	assert(indexStart.size() == indexStop.size());
	for (i = 0; i < indexSize; i ++)
	{
		start = indexStart.at(i);
		stop = indexStop.at(i);

//		if (TestDestroy())
//		{
//			std::cout << "cancel" << std::endl;
//			break;
//		}

//		segments = m_model->GetSegment(start, stop);
//		if (segments == NULL)
//		{
//			// error - end of file?
//			break;
//		}
		spikeDetector(loader, start, stop, countChannels, fs, bandwidth, subOut, subDischarges);
		//continue;
		//delete [] segments;
		//segments = NULL;

		// send progress to main frame
		//m_progress += progressBy;
		//event.SetInt(m_progress);
		//wxQueueEvent((wxEvtHandler*)m_frame, event.Clone());

		// removing of two side overlap detections
		posSize = subOut->m_pos.size();
		disSize = subDischarges->m_MP[0].size();

		if (i > 0)
			tmpFirst = 1;
		else tmpFirst = 0;

		if (i < (int)indexStop.size()-1)
			tmpLast = 1;
		else tmpLast = 0;

		if (posSize > 0)
		{
			if (indexStop.size() > 1)
			{
				for (j = 0; j < posSize; j++)
				{
					if (subOut->m_pos.at(j) < tmpFirst*3*m_settings->m_winsize ||
						subOut->m_pos.at(j) > ((stop - start) - tmpLast*3*m_settings->m_winsize*fs)/fs )
							removeOut.push_back(j);
				}
				subOut->Remove(removeOut);

				for (j = 0; j < disSize; j++)
				{
					minMP = INT_MAX;
					for (k = 0; k < countChannels; k++)
						if (subDischarges->m_MP[k].at(j) < minMP)
								minMP = subDischarges->m_MP[k].at(j);

					if (minMP < tmpFirst*3*m_settings->m_winsize ||
						minMP > ((stop-start) - tmpLast*3*m_settings->m_winsize*fs)/fs )
								removeDish.push_back(j);
				}
				subDischarges->Remove(removeDish);
			}
		}

		posSize = subOut->m_pos.size();
		disSize = subDischarges->m_MP[0].size();
		tmpShift = (indexStart.at(i)+1)/(double)fs - 1/(double)fs;

		// connect out
		for (j = 0; j < posSize; j++)
		{
			m_out->Add(
					subOut->m_pos.at(j) + tmpShift,
					subOut->m_dur.at(j),
					subOut->m_chan.at(j),
					subOut->m_con.at(j),
					subOut->m_weight.at(j),
					subOut->m_pdf.at(j)
				);
		}

		// connect discharges
		for (j = 0; j < countChannels; j++)
		{
			for (k = 0; k < (int)subDischarges->m_MP[j].size(); k++)
			{
				subDischarges->m_MP[j].at(k) += tmpShift;
			}
		}

		for (j = 0; j < countChannels; j++)
		{
			m_discharges->m_MV[j].insert(m_discharges->m_MV[j].end(), subDischarges->m_MV[j].begin(), subDischarges->m_MV[j].end());
			m_discharges->m_MA[j].insert(m_discharges->m_MA[j].end(), subDischarges->m_MA[j].begin(), subDischarges->m_MA[j].end());
			m_discharges->m_MP[j].insert(m_discharges->m_MP[j].end(), subDischarges->m_MP[j].begin(), subDischarges->m_MP[j].end());
			m_discharges->m_MD[j].insert(m_discharges->m_MD[j].end(), subDischarges->m_MD[j].begin(), subDischarges->m_MD[j].end());
			m_discharges->m_MW[j].insert(m_discharges->m_MW[j].end(), subDischarges->m_MW[j].begin(), subDischarges->m_MW[j].end());
			m_discharges->m_MPDF[j].insert(m_discharges->m_MPDF[j].end(), subDischarges->m_MPDF[j].begin(), subDischarges->m_MPDF[j].end());
		}

		// clear
		if (subOut)
			delete subOut;
		if (subDischarges)
			delete subDischarges;
		removeOut.clear();
		removeDish.clear();
	}
}

/// Calculate the starts and ends of indexes for @see #spikeDetector
template<class T>
void Spikedet<T>::getIndexStartStop(wxVector<int64_t>& indexStart, wxVector<int64_t>& indexStop, int64_t cntElemInCh, int64_t T_seg, int fs, int winsize)
{
	int64_t start = 0, end;
	int i, startSize;

	while (start < cntElemInCh)
	{
		indexStart.push_back(start);
		start += T_seg * fs;
	}

	startSize = indexStart.size();
	if (indexStart.size() > 1)
	{
		for (i = 1; i < startSize; i++)
			indexStart.at(i) -= 3*winsize;

		for (i = 0; i < startSize; i++)
		{
			end = indexStart.at(i) + T_seg*fs + 2*(3*winsize);
			indexStop.push_back(end);
		}

		indexStop.front() -= 3*winsize;
		indexStop.back() = cntElemInCh;

		if (indexStop.back() - indexStart.back() < T_seg * fs)
		{
			indexStart.pop_back();
			//indexStart.pop_back(); // No idea what this shit was supposed to mean...
			//indexStop.back() = cntElemInCh;
			auto tmp = indexStop.back();
			indexStop.pop_back();
			indexStop.back() = tmp;
		}
	}
	else
	{
		indexStop.push_back(cntElemInCh);
	}
}

/// Spike detector
template<class T>
void Spikedet<T>::spikeDetector(SpikedetDataLoader<T>* loader, int startSample, int stopSample, const int& countChannels, const int& inputFS, const BANDWIDTH& bandwidth,
									CDetectorOutput*& out, CDischarges*& discharges)
{
	double 				  k1 = m_settings->m_k1;
	double 				  k2 = m_settings->m_k2;
	double 				  discharge_tol = m_settings->m_discharge_tol;
	int 				  decimation = /*m_settings->m_decimation*/decimationF;
	int                   fs = inputFS;

	int    		  		  countRecords = stopSample - startSample;
	wxVector<int> 		  index;
	int 		  		  stop, step;
	int	  	        	  i, j;
	float 				  k;
	int 				  tmp_start;
	COneChannelDetect**   threads;
	ONECHANNELDETECTRET** ret;

	int    	  			  winsize  = m_settings->m_winsize * fs;
	double 	  			  noverlap = m_settings->m_noverlap * fs;

	// OUT
	double 				  t_dur = 0.005;
	wxVector<bool> 		  ovious_M(countRecords, false);
	double 				  position;
	bool 				  tmp_sum = false;

	wxVector<double>** 	  m;
	int 				  tmp_round;
	float 				  tmp_start2, tmp_stop;

	// definition of multichannel events vectors
	wxVector<int>* 		  point = /*new wxVector<int>[2]*/nullptr; // Why the fuck is this here?
	int 				  tmp_old = 0;
	int 				  tmp_act = 0;
	int 				  channel;

	// MV && MA && MW && MPDF && MD && MP
	double 				  tmp_seg;
	double   			  tmp_mv;
	double 				  tmp_max_ma;
	double 				  tmp_max_mw;
	double 				  tmp_max_mpdf;
	double 				  tmp_md;
	double 				  tmp_mp;
	int    				  tmp_row;

	wxVector<T>* data = prepareSegment(loader, startSample, stopSample);

	// If sample rate is > "decimation" the signal is decimated => 200Hz default.
	if (fs > decimation)
	{
		//CDSP::Resample(data, countChannels, fs, decimation);

		fs = decimation;
		winsize  = m_settings->m_winsize * fs;
		noverlap = m_settings->m_noverlap * fs;
		countRecords = data[0].size();
	}

	// Segmentation index
	stop = countRecords - winsize + 1;

	if (noverlap < 1)
		step = round(winsize * (1 - noverlap));
	else
		step = winsize - noverlap;

	for (i = 0; i < stop; i += step)
		index.push_back(i);

	// FILTERING
	// filtering Nx50Hz
	CDSP::Filt50Hz(data, countChannels, fs, m_settings->m_main_hum_freq, bandwidth);

	// filtering 10-60Hz
	CDSP::Filtering(data, countChannels, fs, bandwidth);

	// local maxima detection
	ret = new ONECHANNELDETECTRET*[countChannels];
	threads = new COneChannelDetect*[countChannels];
	for (i = 0; i < countChannels; i++)
	{
		threads[i] = new COneChannelDetect(&data[i], m_settings, fs, &index, i);
		//threads[i]->Run();
		ret[i] = threads[i]->Entry();
	}

	for (i = 0; i < countChannels; i++)
	{
		//ret[i] = (ONECHANNELDETECTRET*)threads[i]->Wait();
		delete threads[i];
	}

	delete[] data;
	//return;

	delete [] threads;
	// processing detection results
	for (i = 0; i < countChannels; i++)
	{
		if (ret[i] == NULL)
			continue;

		//% first and last second is not analyzed (filter time response etc.)
		// first section
		for (j = 0; j < fs; j++)
		{
			ret[i]->m_markersHigh->at(j) = false;
			ret[i]->m_markersLow->at(j) = false;
		}
		// last section
		tmp_start = ret[i]->m_markersHigh->size() - fs - 1;
		for (j = tmp_start; j < (int)ret[i]->m_markersHigh->size(); j++)
		{
			ret[i]->m_markersHigh->at(j) = false;
			ret[i]->m_markersLow->at(j) = false;
		}
	}

	// OUT
	out = new CDetectorOutput();
	for (channel = 0; channel < countChannels; channel++)
	{
		for (j = 0; j < countRecords; j++)
		{
			if (ret[channel] == NULL)
				continue;

			if (ret[channel]->m_markersHigh->at(j) == true)
			{
				ovious_M.at(j) = true;
				position = (j+1)/(double)fs;

				out->Add(position, t_dur, channel + 1, 1, ret[channel]->m_envelopeCdf.at(j), ret[channel]->m_envelopePdf.at(j));
			}
		}
	}

	// ambiguous spike events output
	if (k1 != k2)
	{
		for (channel = 0; channel < countChannels; channel++)
		{
			if (ret[channel] == NULL)
				continue;

			for (j = 0; j < countRecords; j++)
			{
				if (ret[channel]->m_markersLow->at(j) == true)
				{
					if (ret[channel]->m_markersHigh->at(j) == true)
						continue;

					tmp_sum = false;
					for (k = round(j - 0.01*fs); k <= (j - 0.01*fs); k++)
						if(ovious_M.at(k))
							tmp_sum = true;

					if(tmp_sum)
					{
						position = (j+1)/(double)fs;
						out->Add(position, t_dur, channel + 1, 0.5, ret[channel]->m_envelopeCdf.at(j), ret[channel]->m_envelopePdf.at(j));
					}
				}
			}
		}
	}

	// making M stack pointer of events
	m = new wxVector<double>*[countChannels];
	for (i = 0; i < countChannels; i++)
		m[i] = new wxVector<double>(countRecords, 0.0);

	for (i = 0; i < (int)out->m_pos.size(); i++)
	{
		tmp_start2 = out->m_pos.at(i) * fs;
		tmp_stop = out->m_pos.at(i) * fs + discharge_tol * fs;
		for (k = tmp_start2; k <= tmp_stop; k += 1)
		{
			tmp_round = round(k) - 1;
			m[out->m_chan.at(i)-1]->at(tmp_round) = out->m_con.at(i);
		}
	}

	// definition of multichannel events vectors
	delete [] point; // TODO: remove this
	point = new wxVector<int>[2];
	for (i = 0; i < countRecords; i++)
	{
		tmp_act = 0;
		for (j = 0; j < countChannels; j++)
			if (m[j]->at(i) > 0)
				tmp_act = 1;

		if (tmp_old != tmp_act && tmp_act - tmp_old > 0)
			point[0].push_back(i);
		else if (tmp_old != tmp_act && tmp_act - tmp_old < 0)
			point[1].push_back(i-1);

		tmp_old = tmp_act;
	}

	// MV && MA && MW && MPDF && MD && MP
	discharges = new CDischarges(countChannels);
	for (i = 0; i < (int)point[0].size(); i++)
	{
		for (channel = 0; channel < countChannels; channel++)
		{
			if (ret[channel] == NULL) continue;

			tmp_mv 	     = 0;
			tmp_max_ma   = 0;
			tmp_max_mw   = 0;
			tmp_max_mpdf = 0;
			tmp_mp  	 = NAN;
			tmp_row 	 = 0;

			for (j = point[0].at(i) - 1; j < point[1].at(i); j++)
			{
				// MV
				if (m[channel]->at(j) > tmp_mv)
				{
					tmp_mv = m[channel]->at(j);
					// MP
					if (std::isnan(tmp_mp))
						tmp_mp = ((double)tmp_row + point[0].at(i)+1) / (double)fs;
				}

				// MA
				tmp_seg = ret[channel]->m_envelope.at(j) - (ret[channel]->m_prahInt[0].at(j) / (double)k1);
				tmp_seg = std::fabs(tmp_seg);
				if (tmp_seg > tmp_max_ma)
					tmp_max_ma = tmp_seg;

				// MW
				tmp_seg = ret[channel]->m_envelopeCdf.at(j);
				if (tmp_seg > tmp_max_mw)
					tmp_max_mw = tmp_seg;

				// MPDF
				tmp_seg = ret[channel]->m_envelopePdf.at(j) * m[channel]->at(j);
				if (tmp_seg > tmp_max_mpdf)
				   tmp_max_mpdf = tmp_seg;

			   tmp_row++;
			}

			discharges->m_MV[channel].push_back(tmp_mv);
			discharges->m_MA[channel].push_back(tmp_max_ma);
			discharges->m_MPDF[channel].push_back(tmp_max_mpdf);
			discharges->m_MW[channel].push_back(tmp_max_mw);
			discharges->m_MP[channel].push_back(tmp_mp);

			// MD
			tmp_md = (point[1].at(i) - point[0].at(i)) / (double)fs;
			discharges->m_MD[channel].push_back(tmp_md);
		}
	}

	for (i = 0; i < countChannels; i++)
	{
		if (ret[i] == NULL) continue;

		delete ret[i];
		delete m[i];
	}

	delete [] point;
	delete [] m;
	delete [] ret;
}

template<class T>
vector<T>* Spikedet<T>::prepareSegment(SpikedetDataLoader<T>* loader, int start, int stop)
{
	if (decimationF >= fs)
	{
		vector<T>* output = new vector<T>[channelCount];

		stepBuffer.resize(BLOCK_SIZE*channelCount);
		for (int i = 0; i < channelCount; i++)
			output[i].reserve(stop - start);

		for (int i = start; i < stop; i += BLOCK_SIZE)
		{
			int end = min(stop, i + BLOCK_SIZE);
			int len = end - i;

			loader->readSignal(stepBuffer.data(), i, end - 1);

			for (int j = 0; j < channelCount; j++)
			{
				for (int k = 0; k < len; k++)
					output[j].push_back(stepBuffer[j*len + k]);
			}
		}

		return output;
	}

	cl_int err;

	int discard = filterProcessor->discardSamples();
	int delay = filterProcessor->delaySamples();
	int step = BLOCK_SIZE - discard;
	int len = stop - start;
	len = (len + step - 1)/step*step;

	stepBuffer.resize((BLOCK_SIZE + 4)*channelCount);

	segmentBuffer.resize(len*channelCount);
	vector<T*> channelPointers(channelCount);
	for (int i = 0; i < channelCount; i++)
		channelPointers[i] = segmentBuffer.data() + i*len;

	for (int i = 0; i < len; i += step)
	{
		loader->readSignal(stepBuffer.data(), start + i - discard + delay, start + i + delay + step - 1 + 4);

		err = clEnqueueWriteBuffer(queue, inBuffer, CL_TRUE, 0, (BLOCK_SIZE + 4)*channelCount*sizeof(T), stepBuffer.data(), 0, nullptr, nullptr);
		checkClErrorCode(err, "clEnqueueWriteBuffer()");

		filterProcessor->process(inBuffer, outBuffer, queue);

		//OpenCLContext::printBuffer("spikedet_after_filter.txt", outBuffer, queue);

		for (int j = 0; j < channelCount; j++)
		{
			err = clEnqueueReadBuffer(queue, outBuffer, /*CL_FALSE*/CL_TRUE, (j*(BLOCK_SIZE + 4) + discard)*sizeof(T), step*sizeof(T), channelPointers[j], 0, nullptr, nullptr);
			checkClErrorCode(err, "clEnqueueReadBuffer()");

			channelPointers[j] += step;
		}
	}

	//OpenCLContext::printBuffer("spikedet_segment.txt", segmentBuffer.data(), len*channelCount);

	err = clFinish(queue);
	checkClErrorCode(err, "clFinish()");

	// Create the output segment;
	int D = fs/decimationF;
	assert(D >= 1);
	assert(D*decimationF == fs);

	vector<T>* output = new vector<T>[channelCount];
	for (int i = 0; i < channelCount; i++)
	{
		int size = (stop - start)/D;
		output[i].resize(size);
		T* channelPointer = segmentBuffer.data() + i*len;

		for (int j = 0; j < size; j++)
		{
			output[i][j] = *channelPointer;
			channelPointer += D;
		}
	}

	return output;
}

template class Spikedet<float>;
//template class Spikedet<double>; // TODO: try it with doubles

} // namespace AlenkaSignal
