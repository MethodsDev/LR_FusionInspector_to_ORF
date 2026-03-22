#!/bin/bash

set -ev

VERSION=`cat VERSION.txt`

docker push trinityctat/lr_fusioninspector_to_orf:${VERSION} 
docker push trinityctat/lr_fusioninspector_to_orf:latest


