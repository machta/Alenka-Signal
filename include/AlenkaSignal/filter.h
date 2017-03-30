#ifndef ALENKASIGNAL_FILTER_H
#define ALENKASIGNAL_FILTER_H

#include <cmath>
#include <vector>
#include <cstdio>

namespace AlenkaSignal
{

/**
 * @brief A class for computing FIR filter coefficients.
 *
 * After construction the filter has no effect, i.e. it is an all-pass filter.
 *
 * The filter can be configured to work as any combination of the following
 * basic filter types:
 * * low-pass,
 * * high-pass and
 * * notch filter.
 *
 * The filter can be configured by the appropriate set functions.
 *
 * The coefficients are computed using the frequency-sampling method.
 */
template<class T>
class Filter
{
public:
	/**
	 * @brief Filter constructor.
	 * @param M The length of the filter and the number of coefficients that
	 * will be returned.
	 * @param Fs The sampling frequency.
	 */
	Filter(unsigned int M, double Fs, double notchWidth = 3) : M(M), Fs(Fs), notchWidth(notchWidth) {}

	/**
	 * @brief Returns a vector with the coefficients.
	 */
	std::vector<T> computeSamples();

	void lowpass(bool on)
	{
		lowpassOn = on;
	}
	bool lowpass()
	{
		return lowpassOn;
	}
	double getLowpass() const
	{
		return lowpassF*Fs/2;
	}
	void setLowpass(double value)
	{
		lowpassF = value/Fs*2;
	}

	void highpass(bool on)
	{
		highpassOn = on;
	}
	bool highpass()
	{
		return highpassOn;
	}
	double getHighpass() const
	{
		return highpassF*Fs/2;
	}
	void setHighpass(double value)
	{
		highpassF = value/Fs*2;
	}

	void notch(bool on)
	{
		notchOn = on;
	}
	bool notch()
	{
		return notchOn;
	}
	double getNotch() const
	{
		return notchF;
	}
	void setNotch(double value)
	{
		notchF = value/Fs*2;
	}

private:
	unsigned int M;
	double Fs, lowpassF, highpassF, notchF, notchWidth;
	bool notchOn = false, lowpassOn = false, highpassOn = false;
};

} // namespace AlenkaSignal

#endif // ALENKASIGNAL_FILTER_H
