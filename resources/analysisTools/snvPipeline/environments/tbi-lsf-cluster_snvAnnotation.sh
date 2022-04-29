#!/usr/bin/env bash
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#


source $(dirname "${BASH_SOURCE[0]}")/tbi-lsf-cluster.sh

export HTSLIB_INCLUDE_PATH="/tbi/software/x86_64/htslib/htslib-$HTSLIB_VERSION_FOR_HTS_PYTHON/el7/include"

# Have a central installation for the user. This is to avoid DoS on Github!
# It includes the commit-hash to ensure multiple versions don't interfere with each other.
export PYPY_LOCAL_LIBPATH="$HOME/.local/lib/pypy/hts-python/$HTS_PYTHON_COMMIT"

usePypy() {
    # shellcheck disable=SC1073
    if [[ "$(which "$PYPY_OR_PYTHON_BINARY" 2> /dev/null)" == "" ]]; then
        echo "Couldn't find binary '$PYPY_OR_PYTHON_BINARY'" >> /dev/stderr
        exit 1
    elif [[ $("$PYPY_OR_PYTHON_BINARY" --version 2>&1 | grep -P "PyPy") != "" ]]; then
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
    # Note that both the site-packages/ and the repo-dir are specific for the commit.
    # Thus, by selecting the commit one selects the local PyPy libdir!
    local sitePackageDir="$PYPY_LOCAL_LIBPATH/site-packages"
    local HTS_PYTHON_REPODIR="$PYPY_LOCAL_LIBPATH/repository"
    if [[ ! -d "$sitePackageDir" ]]; then
        mkdir -p "$sitePackageDir"
    fi
    local LOCK="$sitePackageDir/hts-egg~"
    lock "$LOCK"

    if [[ ! -d "$(echo "$sitePackageDir/hts-*.egg")" ]]; then
        echo "Installing hts for pypy" >> /dev/stderr
        module load git/"${GIT_VERSION:?GIT_VERSION undefined}"
        export GIT_BINARY=git

        "$GIT_BINARY" clone "$HTS_PYTHON_GIT_REPOSITORY" "$HTS_PYTHON_REPODIR"
        "$GIT_BINARY" -C "$HTS_PYTHON_REPODIR" checkout "$HTS_PYTHON_COMMIT"
        pushd "$HTS_PYTHON_REPODIR"
        ls .eggs/nose* > /dev/null 2>&1
        if [[ $? -ne 0 ]]; then
            pip3 install -t .eggs/ nose==1.3.7
        fi

        # Needed for CFFI
        module switch "htslib/$HTSLIB_VERSION_FOR_HTS_PYTHON"

        C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH" PYTHONPATH="$sitePackageDir" \
          "$PYPY_OR_PYTHON_BINARY" setup.py build
        C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH" PYTHONPATH="$sitePackageDir" \
          "$PYPY_OR_PYTHON_BINARY" setup.py install --prefix="$PYPY_LOCAL_LIBPATH"

        popd
        unset C_INCLUDE_PATH
        module switch "htslib/$HTSLIB_VERSION"
    else
        echo "hts for pypy is already installed" >> /dev/stderr
    fi
    unlock "$LOCK"
}

pypyCopySam() {
   module switch "htslib/$HTSLIB_VERSION_FOR_HTS_PYTHON"
   local sitePackageDir="$PYPY_LOCAL_LIBPATH/site-packages"
   C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH" PYTHONPATH="$sitePackageDir" "$PYPY_OR_PYTHON_BINARY" "$@"
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
