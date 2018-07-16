#!/bin/bash
set -x 
set -e

# small script to compare different spades versions
# Datasets and spades versions can be adjusted
data=Test

for version in 3.5.0 3.6.2 3.7.0 3.8.0 3.8.1 3.9.0 3.10.0 3.10.1  
do
echo "$version" 
if [ ! -d bactofidia_"$version" ]; then
git clone https://github.com/aschuerch/bactofidia bactofidia_"$version"
fi
cd bactofidia_"$version"
git checkout newestversions
ln -s ../*fastq.gz .
for configfiles in config_miseq.yaml.versions  config.yaml.versions  package-list.txt.versions
do 
echo "$configfiles"

filename="${configfiles%.*}"
echo "$filename"
sed s/"spadesversion"/$version/g  "$configfiles"  > "$filename" 
done

./bactofidia.sh "$data"
cd ..
done
