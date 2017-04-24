#!/bin/bash

# Usage: ./download-libraries.sh
#
# This script downloads all dependant libraries.
# Use Git Bash or a similar tool to run this on Windows.

if [ -d clFFT ]
then
	clFFT=skipped
else
	git clone --depth 1 https://github.com/clMathLibraries/clFFT.git -b develop && # TODO: Switch to master for the new release that includes the locale fix.
	clFFT=OK || clFFT=fail
fi

if [ -d unit-test/googletest ]
then
	googletest=skipped
else
	git clone --depth 1 https://github.com/google/googletest.git unit-test/googletest &&
	googletest=OK || googletest=fail
fi

if [ -d alglib ]
then
	alglib=skipped
else
	curl http://www.alglib.net/translator/re/alglib-3.10.0.cpp.gpl.zip > alglib.zip &&
	unzip -q alglib.zip -dalglib &&
	rm alglib.zip &&
	mv alglib/cpp/* alglib &&
	rm -r alglib/cpp &&
	alglib=OK || alglib=fail
fi

if [ -d eigen ]
then
	eigen=skipped
else
	git clone --depth 1 https://github.com/RLovelett/eigen.git &&
	eigen=OK || eigen=fail
fi

if [ -d libsamplerate ]
then
	libsamplerate=skipped
else
	git clone --depth 1 https://github.com/erikd/libsamplerate.git &&
	cd libsamplerate &&
	./autogen.sh &&
	cd - &&
	libsamplerate=OK || libsamplerate=fail
fi

echo
echo ========== Download summary ==========
echo "Library path            Status"
echo ======================================
echo "clFFT                   $clFFT"
echo "unit-test/googletest    $googletest"
echo "alglib                  $alglib"
echo "eigen                   $eigen"
echo "libsamplerate           $libsamplerate"

