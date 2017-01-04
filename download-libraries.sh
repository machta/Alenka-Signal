#!/bin/bash

#git clone --depth 1 https://github.com/clMathLibraries/clFFT.git
git clone --depth 1 https://github.com/machta/clFFT.git -b locale-bug

git clone --depth 1 https://github.com/google/googletest.git unit-test/googletest

curl http://www.alglib.net/translator/re/alglib-3.10.0.cpp.gpl.zip > alglib.zip
unzip alglib.zip -dalglib
rm alglib.zip
mv alglib/cpp/* alglib
rm -r alglib/cpp

git clone --depth 1 https://github.com/RLovelett/eigen.git
