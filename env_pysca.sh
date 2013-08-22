#!/bin/bash

export PYSCA_DIR=$(cd `dirname "${BASH_SOURCE[0]}"` && pwd -P)

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH=${PYSCA_DIR}
else
    export PYTHONPATH=${PYSCA_DIR}:${PYTHONPATH}
fi

export PATH=${PYSCA_DIR}/bin:$PATH
