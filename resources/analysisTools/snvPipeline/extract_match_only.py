#!/usr/bin/env python
# Author: Jeongbin Park
# Purpose: Reporting actual values of 'MATCH=exact'
# Usage: perl annotate_vcf.pl ... --reportMatchType | python -u extract_match_only.py > output.vcf
# Explanation: This script checks whether the value of the last column contains 'MATCH=exact', which is reported
#              by 'annotate_vcf.pl'. If so, the column value will be replaced with its actual value, otherwise
#              with '.' (dot). This script should be used with '-u' option of Python so that it does not use
#              stdin/stdout buffering.

import sys
for line in sys.stdin:
    if line[0] == '#':
        print line.rstrip()
        continue
    es = line.strip().split('\t')
    ees = es[-1].split('&')
    val = '.'
    for ee in ees:
        if 'MATCH=exact;' in ee:
            val = ee.split(';')[-1]
            break
        else:
            continue
    es[-1] = val
    print('\t'.join(es))
