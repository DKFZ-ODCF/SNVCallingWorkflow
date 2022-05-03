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
    # Registering a trap to ensure that the lock is removed also in case of an error. To not
    # pollute the outer trap namespace, the block is put into a subshell by parentheses.
    (
        trap "echo 'Forced unlock: $LOCK' >> /dev/stderr; unlock $LOCK" ERR
        # The trap must be set after the lock was successfully acquired, or it may delete a lock
        # set by another process.

        if [[ $(ls -d $sitePackageDir/hts-*.egg 2>> /dev/null | wc -l) -eq 0 ]]; then
            module load git/"${GIT_VERSION:?GIT_VERSION undefined}"
            export GIT_BINARY=git

            if [[ ! -d "$(echo "$HTS_PYTHON_REPODIR")" ]]; then
                echo "Cloning hts-python repository" >> /dev/stderr
                "$GIT_BINARY" clone "$HTS_PYTHON_GIT_REPOSITORY" "$HTS_PYTHON_REPODIR"
                "$GIT_BINARY" -C "$HTS_PYTHON_REPODIR" checkout "$HTS_PYTHON_COMMIT"
            fi
            pushd "$HTS_PYTHON_REPODIR"

            if [[ $(ls .eggs/nose* 2> /dev/null | wc -l) -eq 0 ]]; then
                echo "Installing nose (needed for installation)" >> /dev/stderr
                pip3 install -t .eggs/ nose==1.3.7
            fi

            # Needed for CFFI
            module switch "htslib/$HTSLIB_VERSION_FOR_HTS_PYTHON"

            echo "Building and installing hts-python for pypy" >> /dev/stderr
            C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH" PYTHONPATH="$sitePackageDir" \
              "$PYPY_OR_PYTHON_BINARY" setup.py build
            C_INCLUDE_PATH="$HTSLIB_INCLUDE_PATH" PYTHONPATH="$sitePackageDir" \
              "$PYPY_OR_PYTHON_BINARY" setup.py install --prefix="$PYPY_LOCAL_LIBPATH"

            popd
            unset C_INCLUDE_PATH
            module switch "htslib/$HTSLIB_VERSION"
        else
            echo "hts-python for pypy is already installed" >> /dev/stderr
        fi
        unlock "$LOCK"
        # Delete the trap
        trap - ERR
    )
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
