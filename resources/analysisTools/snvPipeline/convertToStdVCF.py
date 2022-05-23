#!/usr/bin/env python
#
# Copyright (c) 2022 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#
# Authors: Jules Kerssemakers, Sophia Stahl, Philip Kensche
#
from __future__ import print_function

import gzip
import json
import os
import sys
import textwrap
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# Import tool for converting DKFZ VCF files to vcf files conforming to the specifications of the
# standard VCF version:
vcf_version = 'VCFv4.2'

StdVCF_header_cols = [
    # These first 8 are required
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
    # Format is required when we specify genotype information, followed by patient-ID
    # headernames added in convert()
    'FORMAT']


def read_meta_information(file):

    def dict_without_comments(d):
        return dict(filter(lambda (k, v): not k.startswith("__"), d.items()))

    with open(file, "r") as f:
        raw_infos = json.load(f)
        return dict(map(lambda (k, v): (k, dict_without_comments(v)),
                        dict_without_comments(raw_infos).items()))


def open_maybe_compressed_file(filename, mode="r"):
    """
    Checks if file is compressed or not and opens it accordingly.

    Checks if the last three letters equal ".gz" and uses gzip.open() or open() depending on the
    outcome.

    :param filename: Path to the file that is going to be opened.
    :param mode: Uses read-mode as default.
    :return: Returns the opened fileobject.
    """
    if filename[-3:] == ".gz":
        #print("detected " + filename + " as gzip file", file=sys.stderr)
        return gzip.open(filename, mode)
    else:
        #print("detected " + filename + " as non-zip file", file=sys.stderr)
        return open(filename, mode)


def convert_str_to_dict(str_to_convert, pair_separator, key_val_separator):
    """
    Converts a flat string with key/value pairs into a dictionary.

    :param str_to_convert: String to be converted to a dict
    :param pair_separator: String or character that separates the key-value-pairs from each other
    :param key_val_separator: String or character that separates keys from values
    :return: dictionary that was created from the string
    """
    return dict(
        pair.split(key_val_separator)
        for pair
        in str_to_convert.split(pair_separator))


def convert_dict_to_str(dict_to_convert, pair_separator, key_val_separator):
    """
    Converts a flat dictionary with just a single level of keys and values into a string.
    Pair- and key/value-separators are configurable.

    :param dict_to_convert: Dictionary that is to be converted
    :param key_val_separator: String or character that should separate keys from values in the
           string
    :param pair_separator: String or character that should separate the key-value-pairs from each
           other
    :return: string that was created from the dictionary
    """

    return pair_separator.join(
        key + key_val_separator + val
        for (key,val)
        in dict_to_convert.iteritems())


def infer_vcf_column_indices(file_object):
    """
    Maps the header of the VCF file to indices to make them comfortably accessible.

    Returns a dictionary indices["HEADER"]=index of each column-heading, indices are zero-based,
    so indices[CHROM]===0. This reads the first couple of lines of the file. To do so it `seek()`s
    to the beginning of the file. The filepointer is returned to its starting location afterwards.
    If no header line was found before the first data line (first line without '#'), an empty
    dictionary is returned.

    :param file_object: The vcf file to read the header from.
    :return: Dictionary with the mapped header and indices.
    """

    indices = {}

    # store where we were in the file, and look from the beginning (where headers usually hang out :-) )
    # TODO: seek and tell don't work on compressed files!
    #before_index = file_object.tell()
    #file_object.seek(0)

    for lineNumber, line in enumerate(file_object):
        if line.startswith("##"):
            # "##INFO" or "##FORMAT" line, metadata we don't need; ignore
            continue

        elif line.startswith("#CHROM") or line.startswith("CHROM"):
            # Title line, with all the column headings.
            titleFields = line.lstrip("#").split()

            for column, header in enumerate(titleFields):
                # unfortunately, the INFO_CONTROL header has some parameters in it:
                #   "INFO_control(VAF=variant_allele_fraction;TSR=total_variant_supporting_reads_incl_lowqual)"  # noqa
                # that (..) shouldn't be there in the first place
                # for easier access in the rest of the program: trim the excess stuff before
                # it infects everything.
                if header.startswith("INFO_control"):
                    header = "INFO_control"

                indices[header] = column
            break
        else:
            # normal data line, header not found...
            break

    # erase traces: go back to where we started
    #file_object.seek(before_index)
    return indices


def write_metadata_definitions(keys_to_write, category, output_vcf_object, meta_information):
    """
    Write the annotation for used-and-desired metadata keys to the VCF file meta-information block.
    See also: VCF 4.2 spec 1.2

    :param keys_to_write: The keys whose definition should be written to the new VCF file's
            meta-information section. Every key included MUST have a complete definition
            in meta_information.py.
    :param category: Category of the keys we're writing (one of FORMAT, INFO, FILTER, ALT)
    :param output_vcf_object: opened vcf file object to which the new meta-information is written
    """
    key = None
    try:
        for key in sorted(keys_to_write):
            key_info = meta_information[category][key]

            # If we have defined a 'better' name for this field, use that instead
            key = key_info.get("new_info_id", key)

            if category == "INFO" or category == "FORMAT":
                output_vcf_object.write(
                    "##{category}=<ID={id},Number={number},Type={type},Description=\"{description}\">\n".format(
                        category=category, id=key, **key_info)
                ),
            elif category == "FILTER":
                output_vcf_object.write(
                    "##{category}=<ID={id},Description=\"{description}\">\n".format(
                        category=category, id=key, **key_info)
                )
    except KeyError:
        raise KeyError("Meta-information incomplete for '%s' in meta_information.py" % key)


def update_keys_used(line, col_indices, keys_used):
    """
    Sets of keys used in the chosen columns are updated for each data_line until they contain
    all the meta-information keys used in the entire file contents.

    :param line: Data line of vcf file that is to be converted into the new format
    :param col_indices: Dictionary of column headers as keys and respective indices as values
    :param keys_used: Set of already observed keys, updated in place
    """
    line_as_list      = line.rstrip().split("\t")
    INFO_old_keys     = convert_str_to_dict(line_as_list[col_indices["INFO"]], ';', '=').keys()

    INFO_control_keys = convert_str_to_dict(line_as_list[col_indices["INFO_control"]], ';', '=').keys()
    INFO_control_keys = [ key + "_ctrl" for key in INFO_control_keys ]

    freestanding_INFO_columns = (column_name for (column_name, index) in col_indices.iteritems() if index > 9)

    keys_used["INFO"].update(INFO_old_keys, INFO_control_keys, freestanding_INFO_columns)

    FORMAT_keys = line_as_list[col_indices["FORMAT"]].split(":")
    keys_used["FORMAT"].update(FORMAT_keys)

    FILTER_keys = line_as_list[col_indices["FILTER"]].split(";")
    keys_used["FILTER"].update(FILTER_keys)


def get_all_used_metadata_keys(col_indices, input_file):
    """
    Iterates over the input_file to gather all occurring metadata keys.

    Used keys can come from the sources:
        column ALT    -> tags, only for structural variants, not yet considered by this script
        column FORMAT -> tags, :-separated
        column FILTER -> tags, ,-separated
        column INFO   -> as dict
        column INFO_control -> as dict, mapped into INFO
        freestanding extra columns, mapped into INFO

    :param col_indices: Dictionary of column headers as keys and respective indices as values
    :param input_file: file object opened for reading
    :return: Dictionary with the column categories as keys and a set of meta-info keys used in
             each as values
    """

    keys_used = {
        "INFO": set(),
        "FORMAT": set(),
        "FILTER": set()
    }

    for line in input_file:
        # skip header
        if line.startswith('#'):
            continue

        # unfortunately, we must go through all lines, because keys can be different per line.
        # even the last line may be the first to use a particular FILTER, for example.
        update_keys_used(line, col_indices, keys_used)

    return keys_used


def extract_freestanding_columns_into_dict(line_original, col_indices, meta_information):
    """
    Convert column names and contents to a dictionary that can be inserted into the new INFO field.

    :param line_original: Data line in vcf file that is to be converted
    :param col_indices: Dictionary of column headers as keys and respective indices as values
    :return: Dictionary for the new INFO field with column names as keys and their contents for
             this line as values
    """

    columns_to_extract = [
        column_name for (column_name, details)
        in meta_information["INFO"].items()
        if "new_info_id" in details
    ]

    cols_as_INFO_dict = {}
    for column_name in columns_to_extract:
        try:
            old_column_contents = line_original[col_indices[column_name]]
        except KeyError as e:
            raise RuntimeError("Is this a DKFZ VCF? Avoid the `_raw.vcf.gz`")

        new_INFO_id = meta_information["INFO"][column_name]["new_info_id"]
        cols_as_INFO_dict[new_INFO_id] = old_column_contents
    return cols_as_INFO_dict


def extract_desired_INFO_fields_from_all(raw_INFO_field, used_and_desired_keys, rename_control):
    """
    Convert entries of INFO or INFO_control field to a dict so that keys can be used to retrieve
    only relevant meta-information.

    :param raw_INFO_field: The textual cell contents of an INFO or INFO_control field
    :param used_and_desired_keys: Dictionary with the intersections of available and used
           meta-information keys per category
    :param rename_control: Boolean, if true, all keys are appended with "_ctrl", to differentiate
           between control and tumor values
    :return: Dictionary with contents from the INFO or INFO_control field
    """

    old_INFO_dict = convert_str_to_dict(raw_INFO_field, ";", "=")

    if rename_control:
        for key in old_INFO_dict.keys():
            key_control = key + "_ctrl"
            # remove old 'key' from dict, and replace it under new name 'key_ctrl'
            old_INFO_dict[key_control] = old_INFO_dict.pop(key)

    return dict( filter(
        lambda (k,v): k in used_and_desired_keys,
        old_INFO_dict.iteritems()
    ))


def merge_and_format_INFO_sources(col_indices, line_original, used_and_desired_keys,
                                  meta_information):
    """
    Combines the three sources for the new INFO field (old, control and other columns) into
    first a dictionary, that is then converted into a string.

    :param col_indices: Dictionary of column headers as keys and respective indices as values
    :param line_original: Data line in vcf file that is to be converted
    :param used_and_desired_keys: Dictionary with the intersections of available and used
           meta-information keys per category
    :param meta_information: Meta information configuration dictionary.
    :return: String that replaces the old INFO field with filtered and new information
    """
    new_info = {}
    new_info.update(extract_desired_INFO_fields_from_all(line_original[col_indices["INFO"]],
                                                         used_and_desired_keys["INFO"],
                                                         rename_control=False))
    new_info.update(extract_desired_INFO_fields_from_all(line_original[col_indices["INFO_control"]],
                                                         used_and_desired_keys["INFO"],
                                                         rename_control=True))
    new_info.update(extract_freestanding_columns_into_dict(line_original,
                                                           col_indices,
                                                           meta_information))

    return convert_dict_to_str(new_info, ";", "=")


def convert(input_filename, output_filename, sample_id, meta_information):
    """
    Does the actual conversion of one VCF file format into another (here DKFZ -> standard).

    :param input_filename: VCF file that is to be converted
    :param output_filename: VCF file to which the converted content of the input is written
    :param sample_id: label to use for the sample column (instead of the useless filename in the
           input file)
    :return: None, the output file is the outcome of this function
    """

    # Input file consists of multiple lines of tab-separated strings
    with open_maybe_compressed_file(input_filename, 'r') as input_file:
        # since the column ordering of the in-house columns is never the same, dynamically
        # determine their index in this file.
        col_indices = infer_vcf_column_indices(input_file)
        keys_used_dict = get_all_used_metadata_keys(col_indices, input_file)

    # pre-define some column indices we'll use a lot below
    FORMAT_col_index = col_indices["FORMAT"]
    INFO_col_index   = col_indices["INFO"]
    # By the official VCF spec, all columns following the FORMAT column are sample-columns.
    # In the DKFZ case, there is only one sample column and an unknown number of free-form
    # annotation columns. You shouldn't be using this script if there is more than one sample
    # column.
    SAMPLE_col_index = FORMAT_col_index + 1

    with open_maybe_compressed_file(output_filename, 'w') as output_file:
        # Write Meta-information lines
        output_file.write('##fileformat=%s\n' % vcf_version)  # Required as first line in the file

        # First pass through input file, check which tags/keys occur so that we may write our
        # meta-information block
        used_and_desired_keys = {}
        for category in ["INFO", "FILTER", "FORMAT"]: # ALT not (yet?) considered here
            # Figure out which of all the keys we have documented/available are actually
            # present/used in the current file
            category_keys_available = set(meta_information[category].keys())
            category_keys_used = keys_used_dict[category]
            # only the overlap/intersection is relevant, ignore everything else
            used_and_desired_category_keys = \
                category_keys_available.intersection(category_keys_used)

            write_metadata_definitions(used_and_desired_category_keys,
                                       category,
                                       output_file,
                                       meta_information)

            used_and_desired_keys[category] = used_and_desired_category_keys

        # Second pass: extract directly-copyable meta-information lines, and subsequently
        # convert all data lines
        with open_maybe_compressed_file(input_filename, 'r') as input_file:
            for line in input_file:

                # copy reference info
                if line.startswith("##referenc"):
                    output_file.write(line)

                # copy contig info
                if line.startswith("##contig"):
                    output_file.write(line)

                # CHROM marks the header-line, by definition the last metadata line
                elif line.startswith("CHROM") or line.startswith("#CHROM"):
                    # close off the metadata section with our newer, sexier, shorter, header line
                    output_file.write( "#" +
                                       "\t".join(StdVCF_header_cols) +
                                       "\t" + sample_id +
                                       "\n")

                # other, unhandled, header or meta-information line, ignore
                elif line.startswith("#"):
                    pass

                # Data line: convert
                else:
                    line_original = line.rstrip().split("\t")

                    # First 7 columns left as is
                    new_line = (line_original[:8])

                    # gather all the different tidbits that end up in the final INFO field
                    new_info_str = merge_and_format_INFO_sources(col_indices,
                                                                 line_original,
                                                                 used_and_desired_keys,
                                                                 meta_information)
                    new_line[INFO_col_index] = new_info_str

                    # FORMAT Column copied as is
                    new_line.append(line_original[FORMAT_col_index])

                    # sample information copied as is
                    new_line.append(line_original[SAMPLE_col_index])

                    # Write new data line to output file
                    output_file.write("\t".join(new_line)+"\n")


def parse_options(argv):
    parser = ArgumentParser(description=textwrap.dedent("""
    The VCFs produced by the SNVCallingWorkflow are not standard conform in that some values are not added as additional columns after a single variant column. By contrast, in the standard format, additional columns should only be used to show variants occurring in additional samples.

    This script can be used to convert the DKFZ VCFs to standard VCFs (version 4.2).
    """),
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_file", required=True, type=str,
                        help="Input DKFZ-formatted VCF file")
    parser.add_argument("-s", "--sample-id", dest="sample_id", required=True, type=str,
                        help="Name to use for the sample column in the output VCF")
    parser.add_argument("-o", "--output", dest="output_file", default="/dev/stdout",
                        required=False, type=str,
                        help="Output standard 4.2 VCF file")
    parser.add_argument("-c", "--config", dest="config_file", required=False, type=str,
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "convertToStdVCF.json"),
                        help="Configuration JSON file")
    if len(argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)
    else:
        return parser.parse_args(argv[1:])


if __name__ == '__main__':
    args = parse_options(sys.argv)
    print(" ".join(["Converting", args.input_file, "into", args.output_file,
                    "for", args.sample_id, "using", args.config_file]),
          file=sys.stderr)
    meta_information = read_meta_information(args.config_file)
    convert(args.input_file, args.output_file, args.sample_id, meta_information)
