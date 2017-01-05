#!/bin/bash

if [ -d clFFT ]
then
	clFFT=skipped
else
	#git clone --depth 1 https://github.com/clMathLibraries/clFFT.git &&
	git clone --depth 1 https://github.com/machta/clFFT.git -b locale-bug &&
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

echo
echo ========== Download summary ==========
echo "Library                 Status"
echo ======================================
echo "clFFT                   $clFFT"
echo "unit-test/googletest    $googletest"
echo "alglib                  $alglib"
echo "eigen                   $eigen"
