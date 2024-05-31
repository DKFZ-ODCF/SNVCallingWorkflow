"""
parse_VEP_annotations.py

This script parses VEP annotations from the input and writes the parsed output to STDOUT.
It contains functions to parse VEP format, gene consequences and HGVSc annotations, and format transcript information.

Usage:
    cat VEP_annotated.vcf | python parse_VEP_annotations.py > VEP_annotated_parsed.vcf

"""

import sys

def parse_vep_annotations():
    """
    Parses VEP annotations from the input and writes the parsed output to STDOUT.

    This function reads input from `sys.stdin` line by line and processes each line.
    If a line starts with "#" and contains "##INFO=<ID=CSQ", it extracts the VEP format.
    If a line starts with "#CHROM", it appends the header for the VEP gene consequence and HGVSc columns.
    For all other lines, it splits the line by tab and extracts the info field.
    The info field is then split by semicolon and parsed to extract gene consequence and HGVSc information.
    Finally, the transcript information is formatted and written to STDOUT.

    Args:
        None

    Returns:
        None
    """
    vep_format = {}

    for line in sys.stdin:
        line = line.strip()
        gene_consequence_hgvsc_ds = {}

        if line.startswith("#"):
            if line.startswith("##INFO=<ID=CSQ"):
                vep_format = parse_vep_format(line)
            if line.startswith("#CHROM"):
                line = line + "\tVEP_TRANSCRIPTS"
            # write to STDOUT
            sys.stdout.write("{0}\n".format(line))
        else:
            variant_infos = line.split("\t")
            info = variant_infos[7]
            info = info.split(";")
            gene_consequence_hgvsc_ds = parse_gene_consequence_hgvsc(info, vep_format, gene_consequence_hgvsc_ds)
            line = format_transcript_info(line, gene_consequence_hgvsc_ds)
            
            # write to STDOUT
            sys.stdout.write(line)

def parse_vep_format(line):
    """
    Parses the VEP format line and returns a dictionary mapping each field to its index.

    Args:
        line (str): The VEP format line.

    Returns:
        dict: A dictionary mapping each field to its index.
    """
    vep_format = {}
    vep_format_line = line.split("Format: ")[1]
    vep_format_line = vep_format_line.split("|")
    for i, field in enumerate(vep_format_line):
        vep_format[field] = i
    return vep_format

def parse_gene_consequence_hgvsc(info, vep_format, gene_consequence_hgvsc_ds):
    """
    Parses gene consequences and HGVSc annotations from VEP annotations.

    Args:
        info (list): List of VEP annotations.
        vep_format (dict): Dictionary mapping VEP annotation fields to their indices.
        gene_consequence_hgvsc_ds (dict): Dictionary to store gene consequences and HGVSc annotations.

    Returns:
        dict: Updated gene_consequence_hgvsc_ds dictionary with gene consequences and HGVSc annotations.

    """
    for i in info:
        if i.startswith("CSQ="):
            i = i.replace("CSQ=", "")
            i = i.split(",")
            for j in i:
                j = j.split("|")
                gene_name = j[vep_format["SYMBOL"]]
                consequence = j[vep_format["Consequence"]]
                transcript = j[vep_format["Feature"]]
                hgvsc = j[vep_format["HGVSc"]].split(":")[1] if j[vep_format["HGVSc"]] != "" else ""
                hgvsp = j[vep_format["HGVSp"]].split(":")[1] if j[vep_format["HGVSp"]] != "" else ""
                exon = j[vep_format["EXON"]].split("/")[0] if j[vep_format["EXON"]] != "" else ""
                intron = j[vep_format["INTRON"]].split("/")[0] if j[vep_format["INTRON"]] != "" else ""

                if exon != "" and intron == "":
                    gene_part = "exon" + exon
                elif intron != "" and exon == "":
                    gene_part = "intron" + intron
                else:
                    gene_part = "exon" + exon + "-intron" + intron

                biotype = j[vep_format["BIOTYPE"]]
                impact = j[vep_format["IMPACT"]]

                # Format transcript information with exon number, HGVSc, and HGVSp
                tran_exon_hgvsc_p = "{0}:{1}:{2}:{3}".format(transcript, gene_part, hgvsc, hgvsp)

                if biotype == "protein_coding" and impact in ["HIGH", "MODERATE"]:
                    if gene_name in gene_consequence_hgvsc_ds:
                        if consequence in gene_consequence_hgvsc_ds[gene_name]:
                            gene_consequence_hgvsc_ds[gene_name][consequence].append(tran_exon_hgvsc_p)
                        else:
                            gene_consequence_hgvsc_ds[gene_name][consequence] = [tran_exon_hgvsc_p]
                    else:
                        gene_consequence_hgvsc_ds[gene_name] = {consequence: [tran_exon_hgvsc_p]}
    return gene_consequence_hgvsc_ds

def format_transcript_info(line, gene_consequence_hgvsc_ds):
    """
    Formats the transcript information based on the gene, consequence, and HGVSc data.

    Args:
        line (str): The input line to be formatted.
        gene_consequence_hgvsc_ds (dict): A dictionary containing gene, consequence, and HGVSc data.

    Returns:
        str: The formatted line with transcript information.

    """
    if len(gene_consequence_hgvsc_ds) > 0:
        transcript_info = ""
        for gene in gene_consequence_hgvsc_ds:
            for consequence in gene_consequence_hgvsc_ds[gene]:
                transcript_info += "{0}|{1}({2});".format(gene, consequence, ','.join(gene_consequence_hgvsc_ds[gene][consequence]))
        line = line + "\t" + transcript_info + "\n"
    else:
        line = line + "\t.\n"
    return line


if __name__ == "__main__":
    
    # Parse VEP annotations and write to STDOUT
    parse_vep_annotations()
