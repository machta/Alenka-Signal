#!/bin/bash

function download
{
	curl -L 'https://www.dropbox.com/s/'$1/$2'.zip?dl=0' > $2.zip &&
	unzip -o $2.zip &&
	rm $2.zip
}

data0=skipped

md5sum -c --quiet data.md5 ||
{
	{ download hr3l0c2f79ka8kx alenka-signal-data0 && data0=OK || data0=fail; }
} &&
md5sum -c --quiet data.md5 &&
md5=OK || md5=fail

echo
echo ======= Download summary =======
echo "File                Status"
echo ================================
echo "data0               $data0"
echo "md5sum              $md5"
