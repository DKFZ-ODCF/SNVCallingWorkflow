#!/usr/bin/env python
# Author: Jeongbin Park
# Purpose: Filter out rows from VCF files based on column values
# Usage: vcf_filter_by_crit.py input.vcf[.gz] output.vcf HEADER1 KEY1 CRIT1 HEADER2 KEY2 CRIT2 ...
# Explanation: This script is used to filter out rows which have values that does not meet the criteria.
#              Multiple headers can be specified. KEY is used to get a value from the values in a
#              column (e.g. KEY1=1;KEY2=2). If CRIT is numeric, then it should end with '+' or '-', which specifies
#              'more than' or 'less than' respectively. The row with the value will be filtered out when it meets
#              with the criteria. If CRIT is 'nonexist' then the rows that don't contain KEY will be filtered out.

import sys
from gzip import GzipFile as gzopen

def compare_numeric_value(value, crit_type, crit_value):
    filtered = False
    found_value = float(value) if '.' in value else int(value)
    if (crit_type == '+' and crit_value < found_value) or (crit_type == '-' and crit_value > found_value):
        filtered = True
    return filtered

def main(argv=sys.argv):
    vcffn = argv[1]
    outfn = argv[2]

    headers = argv[3::3]

    header_idxs = []
    with gzopen(vcffn) if vcffn[-3:] == '.gz' else open(vcffn) as f, open(outfn, "w") as fo:
        for line in f:
            if line[:2] == "##":
                fo.write(line)
                continue

            entries = line.strip().split('\t')
            if line[0] == '#': # Header
                fo.write(line)
                entries[0] = entries[0][1:] # remove hash
                for header in headers:
                    for idx, entry in enumerate(entries):
                        if header == entry:
                            header_idxs.append(idx)
                            break
                    else:
                        print("Error: header column %s is not found in vcf file."%header)
                        sys.exit(1)
                filters = zip(header_idxs, argv[4::3], argv[5::3])
            else:
                filtered = False
                for header_idx, filter_keys, crit in filters:
                    entry = entries[header_idx]
                    if entry == '.': continue

                    if crit == "nonexist":
                        for filter_key in filter_keys.split(','):
                            keys = [splitted_value.split('=')[0] for splitted_value in entry.split(';')]
                            if not filter_key in keys:
                                filtered = True
                                break
                    else: # Numeric
                        crit_type = crit[-1]
                        if crit_type == '+' or crit_type == '-':
                            crit_value = crit[:-1]
                            crit_value = float(crit_value) if '.' in crit_value else int(crit_value)
                        else:
                            print("Error: not supported criterion %s."%crit)
                            sys.exit(1)

                        if filter_keys == '.':
                            value = entries[header_idx]
                            if compare_numeric_value(value, crit_type, crit_value):
                                filtered = True
                                break
                        else:
                            filter_key_list= filter_keys.split(',')
                            for key, values in filter(lambda e: len(e) == 2, [entry.split('=') for entry in entries[header_idx].split(';')]):
                                if key in filter_key_list:
                                    for value in values.split(','): # There could be multiple values assigned for one key
                                        if compare_numeric_value(value, crit_type, crit_value):
                                            filtered = True
                                            break
                                    if filtered:
                                        break
                    if filtered:
                        break
                if not filtered:
                    fo.write(line)

if __name__ == '__main__':
    main()