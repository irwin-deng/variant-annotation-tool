"""
Read each variant from a VCF file, using PyVCF3 to parse each line

Author: Irwin Deng
"""

import logging
import vcf
from typing import Any, Iterator


logger = logging.getLogger(__name__)

class VcfReader:
    """
    Read each variant from a VCF file, using PyVCF3 to parse each line
    """
    vcf_reader: vcf.Reader

    def __init__(self, vcf_file: str):
        """
        Args:
            vcf_file: path to input VCF file.
        """

        try:
            self.vcf_reader = vcf.Reader(filename=vcf_file)
        except FileNotFoundError:
            logger.error("Input VCF file not found: %s", vcf_file)
            raise
        except IOError as e:
            logger.error("Error reading from input VCF file: %s", str(e))
            raise

    @staticmethod
    def get_coverage_info(record: vcf.model._Record, index: int) -> dict[str, Any]:
        """
        Extract depth of sequence coverage and number of reads for each
        variant, and then calculate percentage of reads supporting the
        variant

        Args:
            variant: A VCF record object
            index: the index of the variant within the record

        Returns:
            dictionary containing coverage info
        """

        try:
            num_reads = record.samples[0]["NR"][index]
            variant_reads = record.samples[0]["NV"][index]
            variant_proportion = (variant_reads / num_reads) if num_reads > 0 else 0
        except (IndexError, KeyError):
            logger.warning("Could not find coverage info in record %s", record)
            return {}

        return {
            "num_reads": num_reads,
            "variant_reads": variant_reads,
            "variant_percentage": variant_proportion * 100
        }

    def __iter__(self) -> Iterator[dict[str, Any]]:
        """
        Iterate through each variant of each line of the input VCF file

        Yields:
            dictionary containing variant info
        """

        for record in self.vcf_reader:
            for i, alt in enumerate(record.ALT):
                variant_info = {
                    "chrom": record.CHROM,
                    "pos": record.POS,
                    "ref": record.REF,
                    "alt": str(alt),
                    "hgvs": f"{record.CHROM}:g.{record.POS}{record.REF}>{alt}"
                }

                variant_info.update(self.get_coverage_info(record, i))
                yield variant_info

    def close(self):
        """Close the reader"""
        try:
            self.vcf_reader.close()
        except AttributeError:
            pass
