#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker build --progress=plain -t trinityctat/lr_fusioninspector_to_orf:${VERSION} .
docker build -t trinityctat/lr_fusioninspector_to_orf .

