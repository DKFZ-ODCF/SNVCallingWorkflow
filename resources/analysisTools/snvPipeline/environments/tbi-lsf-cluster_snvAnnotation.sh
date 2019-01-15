#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#


source $(dirname "${BASH_SOURCE[0]}")/tbi-lsf-cluster.sh

export HTSLIB_INCLUDE_PATH="/tbi/software/x86_64/htslib/htslib-${HTSLIB_VERSION_FOR_HTS_PYTHON}/el7/include"
export PYPY_LOCAL_LIBPATH="$HOME/.local/lib/pypy"
export HTS_PYTHON_COMMIT="cd35d6f15221bb09c218d0f978c0b92c560fcdd6"

usePypy() {
    if [[ $("$PYPY_OR_PYTHON_BINARY" --version 2>&1 | grep -P "PyPy") != "" ]]; then
        echo true
    else
        echo false
    fi
}

lock() {
    if [[ -z `which lockfile` ]]; then
       echo "lockfile-{create,remove} is not well tested. Please install 'lockfile' from procmail package."
       lockfile-create "$1"
    else
       lockfile -s 1 -r 50 "$1"
    fi
}

unlock() {
    if [[ -z `which lockfile` ]]; then
       echo "lockfile-{create,remove} is not well tested. Please install 'lockfile' from procmail package."
       lockfile-remove "$1"
    else
       rm -f "$1"
    fi
}


## This is not provided as a
cloneAndBuildHtsPython() {
    local sitePackageDir="$PYPY_LOCAL_LIBPATH/site-packages"
    local LOCK="$sitePackageDir/hts-egg~"
    lock "$LOCK"
    if [[ ! -d "$sitePackageDir" ]]; then
        mkdir -p "$sitePackagedir"

        module load git/"${GIT_VERSION:?GIT_VERSION undefined}"
        export GIT_BINARY=git

        # Needed for CFFI
        export C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH"
        module switch "htslib/$HTSLIB_VERSION_FOR_HTS_PYTHON"

        if [[ ! -d `echo "$sitePackageDir/hts-*.egg"` ]]; then
            echo "Installing hts-python on your local directory..."
            "$GIT_BINARY" clone "https://github.com/pjb7687/hts-python" "$RODDY_SCRATCH/hts-python"
            "$GIT_BINARY" -C "$RODDY_SCRATCH/hts-python" checkout "$HTS_PYTHON_COMMIT"
            # Below assume that 'nose' is already installed with PyPy
            pushd ${RODDY_SCRATCH}/hts-python
            "$PYPY_LOCAL_LIBPATH" "$PYPY_OR_PYTHON_BINARY" setup.py build
            "$PYPY_LOCAL_LIBPATH" "$PYPY_OR_PYTHON_BINARY" setup.py install --prefix="$PYPY_LOCAL_LIBPATH"
            popd
        fi
        unset C_INCLUDE_PATH
        module switch "htslib/$HTSLIB_VERSION"
    fi
    unlock "$LOCK"
}

pypyCopySam() {
   C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH"
   module switch "htslib/$HTSLIB_VERSION_FOR_HTS_PYTHON"
   PYTHONPATH="$PYPY_LOCAL_LIBPATH/site-packages" "$PYPY_OR_PYTHON_BINARY" "$@"
   unset C_INCLUDE_PATH
   module switch "htslib/$HTSLIB_VERSION"
}

if [[ $(usePypy) == "true" ]]; then
    cloneAndBuildHtsPython
    export -f pypyCopySam
    export PYPY_WITH_PYSAM=pypyCopySam
else
    export PYPY_WITH_PYSAM="$PYPY_OR_PYTHON_BINARY"
fi

filterPeOverlap() {
    "$PYPY_WITH_PYSAM" -u "$TOOL_FILTER_PE_OVERLAP" "$@"
}
export -f filterPeOverlap
