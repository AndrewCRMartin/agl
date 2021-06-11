#!/bin/bash

dest=${HOME}/bin

if [ "X$1" != "X" ]; then
    dest=$1
fi

datadest=${dest}/share/agl/data

mkdir -p $dest
if [ ! -d $dest ]; then
    echo "Unable to create destination directory: $dest"
    exit 1
fi

mkdir -p $datadest
if [ ! -d $datadest ]; then
    echo "Unable to create data directory: $datadest"
    exit 1
fi


(cd src; make)

echo -n "Building data..."
./makedb.pl
echo "done"

echo -n "Copying files to ${dest}..."
cp src/agl $dest
cp -p share/agl/data/* $datadest
echo "done"
