#!/usr/bin/env python
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#


import sys
import re
import argparse
from Bio import bgzf

class BGZFType(argparse.FileType):
    def __call__(self, string):
        # the special argument "-" means sys.std{in,out}
        if string == '-':
            if 'r' in self._mode:
                return sys.stdin
            elif 'w' in self._mode:
                return sys.stdout
            else:
                raise ValueError('argument "-" with mode %r' % self._mode)

        # all other arguments are used as file names
        try:
            if string[-3:] == ".gz":
                if 'r' in self._mode:
                    return bgzf.BgzfReader(string, self._mode)
                elif 'w' in self._mode or 'a' in self._mode:
                    return bgzf.BgzfWriter(string, self._mode)
            else:
                return open(string, self._mode, self._bufsize, self._encoding, self._errors)
        except OSError as e:
            raise ArgumentTypeError("can't open '%s': %s" % (string, e))

class LineParser(object):
    def __init__(self, entries, header_indices):
        self.data = {}
        for header in header_indices:
            if header_indices[header] != -1:
                self.data[header] = entries[header_indices[header]]

    def __getitem__(self, header):
        if header[-6:] == "_VALID":
            entry = self.data[header[:-6]]
            return True if entry != "0" and entry != "." else False
        else:
            return self.data[header]

    def __setitem__(self, header, value):
        self.data[header] = value

def load_headername_from_conf_file(configfile):
    if configfile is None:
        return {}

    p = re.compile(r"^.+_COL=.+$", re.MULTILINE)
    headers = p.findall(configfile.read())

    return dict(header.split("=") for header in headers)

def get_header_indices(headers, configfile, fixed_headers, variable_headers={}):
    header_indices = {}
    header_names_from_conf = load_headername_from_conf_file(configfile)

    for header_key, variable_header in variable_headers.items():
        conf_header = header_names_from_conf.get(header_key, "")
        if conf_header == "":
            selected_header = variable_header
        else:
            selected_header = "^" + conf_header + "$"

        found = False
        for idx, header in enumerate(headers):
            m = re.search(selected_header, header)
            if m is not None:
                header_indices[header_key] = idx
                found = True
        if not found:
            header_indices[header_key] = -1

    for fixed_header in fixed_headers:
        found = False
        for idx, header in enumerate(headers):
            m = re.search(fixed_header, header)
            if m is not None:
                header_indices[m.group(0)] = idx
                found = True
        if not found:
            header_indices[fixed_header.strip("$").strip("^")] = -1

    return header_indices