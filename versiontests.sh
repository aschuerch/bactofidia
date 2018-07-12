#!/bin/bash
set -x 
set -e

# small script to compare different spades versions
# Datasets and spades versions can be adjusted
data=Test

for version in 3.5.0 3.6.2 3.7.0 3.8.0 3.8.1 3.9.0 3.10.0 3.10.1 3.11.0 3.11.1 3.12.0  
do
echo "$version" 
for configfiles in config_miseq.yaml.versions  config.yaml.versions  package-list.txt.versions
do 
echo "$configfiles"

filename="${configfiles%.*}"
echo "$filename"
sed s/"spadesversion"/$version/g  "$configfiles"  > "$filename" 
./bactofidia.sh "$data" && mv results finished_"$version"
done
done
