# Check DataSets Genome Assemblies

The NCBI [```datasets```](https://github.com/ncbi/datasets) tool is very handy, but I've noticed that for large datasets some of the FNA files are truncated.  This may be more common when doing a direct download rather than using the ```--dehydrate``` and ```--rehydrate``` workflows.

This script was written to be used immediately after unarchiving a ```datasets``` download.  It:
    1. Scans the directory tree for the associated ```assembly_data_report.jsonl``` and ```dataset_catalog.json``` files
    2. Sets up expectations for the total number of sequences and total sequence length per assembly
    3. Compares those expecations to what's in the acutal downloaded files
    4. Complains about any mismatches

## Usage
```
check-dsga.py [-h] [-d DIRECTORY] [-V] [-t] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Directory containing unzipped dataset
  -V, --version         Show version
  -t, --tsv_out         Output results in a format suitable for ingestion as a TSV Otherwise output is formatted for easier reading.
  -v, --verbose         Display parsing progress
```

## Current Behaviour

Outputs either:
1. If no issues found, no output (apart from optional verbose mode)
2. If mismatches are found: Mismatched accession numbers, expected and found values, and filenames
   - If the associated file is not found, there will be an exception thrown (#TODO handle more gracefully)
   - n.b. This does not check if there are fna files in the extracted set which do not appear in the catalog or
           assembly data report

By default:
1. runs in the toplevel directory of an unzipped ncbi_dataset
   - (i.e a directory containing a 'ncbi_datset/data/<assembly>/<seqfile.fna>' hierarchy
2. gets the expected values for each assembly from the assembly data report.
3. gets the filename and location from the dataset catalog
4. outputs results to STDOUT (pretty or machine-friendly)
    

## Potential improvements
1. Specifying alternate locations for  assembly report and datset catalog
2. Specifying a custom catalog associating accession numbers with files   
   - (e.g. for work after flattening, renaming, or otherwise altering the unzipped download)
3. Specifying a custom set of expected values
4. Report any accession number that is not in both the catalog and assembly report
5. Handling compressed dataset
6. Better automated testing

contact: joe.e.weaver@gmail.com
