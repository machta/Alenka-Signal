#!/bin/bash

# Downloads the test files from Dropbox.
# This is a workaround for GitHub's screwed up LFS system.

function download
{
	curl -L 'https://www.dropbox.com/s/'$1/$2'.zip?dl=0' > $2.zip &&
	unzip -o $2.zip &&
	rm $2.zip
}

data0=skipped
data1=skipped
data2=skipped

md5sum -c --quiet data.md5 ||
{
	{ download hr3l0c2f79ka8kx alenka-signal-data0 && data0=OK || data0=fail; }
	{ download e5louj39081tnju alenka-signal-data1 && data1=OK || data1=fail; }
	{ download rlroatzft75m4d9 alenka-signal-data2 && data2=OK || data2=fail; }
} &&
md5sum -c --quiet data.md5 &&
md5=OK || md5=fail

echo
echo ======= Download summary =======
echo "File                Status"
echo ================================
echo "data0               $data0"
echo "data1               $data1"
echo "data2               $data2"
echo "md5sum              $md5"

