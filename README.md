# Variant annotation tool

Program to annotate each variant in a given VCF file

Each variant is annotated with each of the following:
1. Depth of sequence coverage at the site of variation
1. Number of reads supporting the variant
1. Percentage of reads supporting the variant versus those supporting reference reads
1. The gene of the variant, type of variation (substitution, insertion, etc.)
and their effect (missense, silent, intergenic, etc.)
1. The minor allele frequency of the variant

## Getting started

Install dependencies and run the tool as follows, replacing ``input_vcf`` and ``output_tsv`` with
the file paths to the input VCF file and output TSV file, respectively:

```shell
pip install -r requirements.txt
python annotate_variants.py <input_vcf> <output_tsv>
```

For example:

```shell
pip install -r requirements.txt
python annotate_variants.py test_vcf_data.txt output.tsv
```

## Files

* ``annotate_variants.py``: variant annotation script
* ``vcf_reader.py``: read each variant from a VCF file
* ``ensembl_api_client.py``: asynchronous client for Ensembl API

## Output file format

The annotated file output is a TSV file containing one variant per row with the following columns:
1. ``chrom``: the chromosome containing the variant
1. ``pos``: the position of the variant in the chromosome
1. ``ref``: the reference base(s) at the given position
1. ``alt``: the variant base(s)
1. ``hgvs``: the HGVS nomenclature
1. ``num_reads``: the total number of reads containining variant location in the sample
1. ``variant_reads``: number of reads containing variant in this sample
1. ``variant_percentage``: percent of reads containing variant in this sample
1. ``gene_id``: ID of affected gene
1. ``gene_symbol``: symbol of affected gene
1. ``biotype``: biotype of transcript or regulatory feature
1. ``impact``: subjective impact classification of consequence type
1. ``most_severe_consequence``: most secure consequence according to Ensembl API
1. ``minor_allele_frequency``: minor allele frequency of the variant

Note: some information was not found from the API. These values appear as blanks in the output file

## Compatibility

This code was developed with Ubuntu 22.04 and Python 3.11.9
