"""
Check datasets assemblies (i.e. fna files) downloaded from NCBI.
Compares the expected number of contigs and total sequence length vs. what's contained in the file.

Outputs either:
    1. If no issues found, no output (apart from optional verbose mode)
    2. If mismatches are found: Mismatched accession numbers, expected and found values, and filenames
        a. If the associated file is not found, there will be an exception thrown (#TODO handle more gracefully)
        b. n.b. This does not check if there are fna files in the extracted set which do not appear in the catalog or
           assembly data report

By default:
    1. runs in the toplevel directory of an unzipped ncbi_dataset
        a. (i.e a directory containing a 'ncbi_datset/data/<assembly>/<seqfile.fna>' hierarchy
    2. gets the expected values for each assembly from the assembly data report.
    3. gets the filename and location from the dataset catalog
    4. outputs results to STDOUT (pretty or machine-friendly)

Potential improvements include:
    1. Specifying alternate locations for  assembly report and datset catalog
    2. Specifying a custom catalog associating accession numbers with files
       (e.g. for work after flattening, renaming, or otherwise altering the unzipped download)
    3. Specifying a custom set of expected values
    4. Report any accession number that is not in both the catalog and assembly report
    5. Handling compressed dataset
    6. Better automated testing

contact: joe.e.weaver@gmail.com
"""

import os
import sys
import json
from argparse import ArgumentParser
from Bio import SeqIO
import pandas as pd

VERSION = "0.0.1"


def parse_args():
    """
    Standard argument parsing
    :return:
    """
    parser = ArgumentParser(description="Check datasets assemblies (i.e. fna files) downloaded from NCBI. "
                                        f'version {VERSION}')
    parser.add_argument("-d", "--directory", type=str, default="", help="Directory containing unzipped dataset")
    parser.add_argument("-V", "--version", action='store_true', help="Show version")
    parser.add_argument("-t", "--tsv_out", required=False, action='store_true',
                        help="Output results in a format suitable for ingestion as a TSV "
                        "Otherwise output is formatted for easier reading.")
    parser.add_argument("-v", "--verbose", required=False, action='store_true',
                        help="Display parsing progress")
    return parser.parse_args()


def main():
    """
    Determine expected FNA stats, read what was actually downloaded, and report on mismatches
    :return:
    """

    # just show the version and return
    args = parse_args()
    if args.version:
        print(VERSION)
        exit(0)

    # confirm we have the assembly data report and catalog
    expected_values_fname = os.path.join(args.directory, "ncbi_dataset", "data", "assembly_data_report.jsonl")
    expected_catalog_fname = os.path.join(args.directory, "ncbi_dataset", "data", "dataset_catalog.json")
    confirm_file_present(expected_values_fname)
    confirm_file_present(expected_catalog_fname)

    # determine expected FNA file stats from assembly data report
    expectations = read_expected_values(expected_values_fname, args.verbose)

    # associate accession numbers with filenames and paths
    accession_files = associate_accessions_with_files(expected_catalog_fname, args.verbose)

    # determine observed FNA file stats from downloaded FNA files
    observations = read_observed_values(accession_files, args.directory, args.verbose)

    # find any mismatches between observed and expected total sequence lengths or number of contigs
    mismatches = determine_mismatches(expectations, observations)

    # report on mismatches
    output_results(mismatches, args.verbose, args.tsv_out)


def confirm_file_present(file_name: str):
    """
    Checks if the file is present. Compalins to STDERR if not and exits with an error code
    :param file_name:
    :return:
    """
    if not os.path.isfile(file_name):
        print(f'Could not locate expected values file: {file_name}', file=sys.stderr)
        exit(1)


def read_expected_values(expected_values_fname: str, verbose: bool) -> pd.DataFrame:
    """
    Given the FNA information in the assembly data record, set our expectations for FNA file statistics
    :param expected_values_fname: Location of assembly data record
    :param verbose: Provide extra output while processing
    :return: Dataframe associating accessions with total number of seqs and total seq length
    """
    if verbose:
        print('Reading assembly data report')
    # get our expected accessions, sequence counts, and sequence lengths from the assembly report
    parsed_expectations = []
    with open(expected_values_fname, 'r') as f:
        for line in f.readlines():
            entry = json.loads(line)
            parsed_expectations.append((entry['accession'],
                                        int(entry['assemblyStats']['numberOfComponentSequences']),
                                        int(entry['assemblyStats']['totalSequenceLength'])))
    expectations = pd.DataFrame(parsed_expectations, columns=['Accession', 'NumSeqsExp', 'TotalSeqLenExp'])
    if verbose:
        print(f'\tRead {expectations.shape[0]} assembly expectations')
    return expectations


def associate_accessions_with_files(expected_catalog_fname: str, verbose: bool) -> dict:
    """
    Given a dataset catalog, list the FNA files associated with each accession
    :param expected_catalog_fname: Filename of dataset catalog to use
    :param verbose: Provide extra output while processing
    :return: Dictionary where each key is an accession number and each file is the location of the FNA file
    """
    if verbose:
        print('Getting FNA files associated with accessions')
    accession_files = {}
    with open(expected_catalog_fname, 'r') as f:
        cat_item = json.load(f)
        for assembly in cat_item['assemblies']:
            # not every element in assemblies relates to an accession, c.f. first etnry which points to the data report
            if 'accession' in assembly:
                for file_nfo in assembly['files']:
                    if file_nfo['fileType'] == "GENOMIC_NUCLEOTIDE_FASTA":
                        # shouldn't have more than one file per assembly
                        if assembly['accession'] in accession_files:
                            print(f'More than one file listed is associated with accession: {assembly["accession"]}',
                                  file=sys.stderr)
                            exit(1)
                        accession_files[assembly['accession']] = file_nfo['filePath']
    if verbose:
        print(f'\tAssociated {len(accession_files.keys())} accessions with files.')
    return accession_files


def read_observed_values(accession_files: dict, directory: str, verbose: bool) -> pd.DataFrame:
    """
    Determine the stats for the FNA files as downloaded (e.g. sequences per accession)
    :param accession_files: Dictionary with accession number as key and FNA file as value
    :param directory: Location of toplevel directory containing unzipped ncbi_datasets
    :param verbose: Provide extra output while processing
    :return: Dataframe listing accession, observed number of sequences, total seq length, and filename of FNA
    """
    if verbose:
        print(f'Processing received FNA files.')
    seqcount = 0
    fna_stats = []
    for k, v in accession_files.items():
        seqname = os.path.join(directory, "ncbi_dataset", "data", v)
        lens = 0
        for i, record in enumerate(SeqIO.parse(seqname, "fasta")):
            lens = lens + len(record.seq)
        nseq = i + 1
        fna_stats.append((k, nseq, lens, v))
        seqcount = seqcount + 1
        if verbose and (not seqcount % 100):
            print(f'\tProcessed {seqcount} FNA files so far.')
    if verbose:
        print(f'\tProcessed {seqcount} total FNA.')
    return pd.DataFrame(fna_stats, columns=['Accession', 'NumSeqsObs', 'TotalSeqLenObs', 'File'])


def determine_mismatches(expectations: pd.DataFrame, observations: pd.DataFrame) -> pd.DataFrame:
    """
    Determine if there's any mismatches between expectations and observations
    :param expectations: Dataframe with expected # seqs and total seq lengths per accession
    :param observations: Dataframe with observed # seqs and total seq lengths per accession
    :return: Dataframe of all accession where expectations and observations do not match
    """
    combined = expectations.merge(observations, on='Accession', how='left')
    mismatches = combined[(combined['NumSeqsObs'] != combined['NumSeqsExp']) |
                          (combined['TotalSeqLenObs'] != combined['TotalSeqLenExp'])]
    return mismatches


def output_results(mismatches: pd.DataFrame, verbose: bool, tsv_format: bool):
    """
    Report any mismatches
    :param mismatches: Dataframe listing accessions with mismatched expectations and observations
    :param verbose: Provide additional output during processing
    :param tsv_format: Should output be in TSV format or more human-readable
    :return:
    """
    if verbose:
        print(f'Found {mismatches.shape[0]} mismatched accessions.')

    for _, mismatch in mismatches.iterrows():
        # creating these mainly for shorter fstrings
        accno = mismatch["Accession"]
        exp_nseq = mismatch["NumSeqsExp"]
        exp_seqlen = mismatch["TotalSeqLenExp"]
        obs_nseq = mismatch["NumSeqsObs"]
        obs_seqlen = mismatch["TotalSeqLenObs"]
        dif_nseq = mismatch["NumSeqsExp"] - mismatch["NumSeqsObs"]
        dif_seqlen = mismatch["TotalSeqLenExp"] - mismatch["TotalSeqLenObs"]
        fname = mismatch["File"]

        if tsv_format:
            print(f'{accno}\t{exp_nseq}\t{obs_nseq}\t{dif_nseq}\t{exp_seqlen}\t{obs_seqlen}\t{dif_seqlen}\t{fname}')
        else:
            print(f'Accession {accno} mismatch. Expected/Observed Seqs: {exp_nseq}/{obs_nseq}(delta:{dif_nseq})  '
                  f'Expected/Observed Total Sequence Lengths: {exp_seqlen}/{obs_seqlen} (delta: {dif_seqlen})'
                  f' File: {fname}')


if __name__ == "__main__":
    main()
