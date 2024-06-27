"""
Script to read VCF file, annotate variants, and output results to TSV file

Usage:
    python annotate_variants.py input_vcf output_tsv

Author: Irwin Deng
"""

import argparse
import csv
import logging
import sys
import asyncio
import aiohttp
from typing import Any
from vcf_reader import VcfReader
from ensembl_api_client import EnsemblApiClient
from tqdm import tqdm

logger = logging.getLogger(__name__)


def write_to_tsv(annotated_variants: list[dict[str, Any]], output_tsv: str) -> None:
    """
    Write annotated variants to TSV file

    Args:
        annotated_variants: list of annotated variants
        output_tsv: filepath to output TSV file
    """

    if not annotated_variants:
        logger.info("No annotated variants to write")
        return

    fieldnames = annotated_variants[0].keys()
    try:
        with open(output_tsv, "w", newline="", encoding="utf-8") as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(annotated_variants)
        logger.info("Annotated variants written to %s", output_tsv)
    except IOError as e:
        logger.error("Error writing to output TSV file: %s", str(e))
        raise


async def process_variants(vcf_reader: VcfReader) -> list[dict[str, Any]]:
    """
    Annotate each variant using the VEP API and write the results to a TSV file.

    Args:
        vcf_reader: VCF reader object.

    Returns:
        list of annotated variants
    """

    api_client = EnsemblApiClient()

    async with aiohttp.ClientSession() as session:
        # Create a list of tasks for each variant
        tasks = [api_client.update_variant_info(session, variant) for variant in vcf_reader]
        annotated_variants = []

        # Use asyncio.as_completed to get futures as they complete
        for future in tqdm(asyncio.as_completed(tasks), total=len(tasks), desc="Annotating variants"):
            # Await the future and append the result to the list of annotated variants
            result = await future
            annotated_variants.append(result)

    return annotated_variants


def main():
    """Parse command-line arguments and run variant annotation process"""

    parser = argparse.ArgumentParser(description="Annotate variants from a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_tsv", help="Output TSV file")
    args = parser.parse_args()

    try:
        vcf_reader = VcfReader(args.input_vcf)
        annotated_variants = asyncio.run(process_variants(vcf_reader))
        vcf_reader.close()
        write_to_tsv(annotated_variants, args.output_tsv)
    # I/O errors already logged in VcfReader and write_to_tsv()
    except IOError:
        sys.exit(1)
    # Catch any unexpected errors
    except Exception as e:
        logger.error("An unexpected error occurred: %s", str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
