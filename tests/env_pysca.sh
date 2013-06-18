#!/bin/bash

__base_dir=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd -P)

if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH=${__base_dir}
else
    export PYTHONPATH=${__base_dir}:${PYTHONPATH}
fi

export PATH=$PATH:${__base_dir}/bin

unset __base_dir
