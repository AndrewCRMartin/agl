#!/bin/bash

dest=${HOME}/bin

if [ "X$1" == "X-h" ]; then
    echo ""
    echo "Usage: ./install.sh [destination]"
    echo ""
    echo "If destination is not provided then software will be installed"
    echo "in $HOME/bin"
    echo ""
    exit 0
fi

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

echo "Grabbing the latest D-SEGMENT data"
(cd util; ./grab_imgt.sh)
echo "done"

echo -n "Building data..."
./util/makedb.pl
echo "done"

echo -n "Copying files to ${dest}..."
cp agl $dest
cp -p share/agl/data/* $datadest
echo "done"
