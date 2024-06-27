"""
Asynchronous client for VEP API, with rate limiting and error handling

Implementation based on Ensembl's example Python client:
https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client

Author: Irwin Deng
"""

import asyncio
import logging
import time
from typing import Any, Optional
from urllib.parse import urlencode
import aiohttp

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


class EnsemblApiClient:
    """
    Asynchronous client for VEP API, with rate limiting and error handling
    """
    server: str             # The base URL for the Ensembl API server
    reqs_per_sec: float     # Maximum number of requests sent per second
    req_count: int          # Current count of requests made
    last_req: float         # Timestamp of the last request
    lock: asyncio.Lock      # Lock to ensure rate limiting

    def __init__(self, server: str = "https://grch37.rest.ensembl.org",
                 reqs_per_sec: int = 10):
        """
        Args:
            server: the base URL for the Ensembl API server
            reqs_per_sec: maximum number of requests allowed per second
        """
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.lock = asyncio.Lock()

    async def perform_rest_action(self, session: aiohttp.ClientSession, endpoint: str,
                                  params: Optional[dict] = None) -> dict:
        """
        Call Ensembl's REST API

        Args:
            session: the aiohttp session for making requests
            endpoint: the API endpoint to query
            params: query parameters for the request

        Returns:
            the JSON response, or empty dict if the response failed
        """
        url = f"{self.server}{endpoint}"
        if params:
            url += "?" + urlencode(params)

        await self._apply_rate_limit()

        try:
            async with session.get(url, headers={"Content-Type": "application/json"}) as response:
                response.raise_for_status()
                return (await response.json())[0]
        except aiohttp.ClientResponseError as e:
            # check if we are being rate limited by the server
            if e.status == 429:
                if "Retry-After" in e.headers:
                    retry = e.headers["Retry-After"]
                    logger.warning("Rate limit exceeded. Retrying after %s seconds.", retry)
                    await asyncio.sleep(float(retry))
                    return await self.perform_rest_action(session, endpoint, params)
            elif e.status == 400:
                logger.info("Endpoint not found: %s", endpoint)
            else:
                logger.error("Request failed for %s: %s - %s", endpoint, e.status, e.message)
        except Exception as e:
            logger.error("Error making request for %s: %s", endpoint, str(e))
        return {}

    async def _apply_rate_limit(self) -> None:
        """
        Apply rate limiting to API requests, ensuring the number of
        requests does not exceed self.reqs_per_sec over the past second
        """
        async with self.lock:
            if self.req_count >= self.reqs_per_sec:
                delta = time.time() - self.last_req
                if delta < 1:
                    await asyncio.sleep(1 - delta)
                self.last_req = time.time()
                self.req_count = 0
            self.req_count += 1

    async def update_variant_info(self, session: aiohttp.ClientSession,
                                  variant: dict[str, Any]) -> dict[str, Any]:
        """
        Update variant using Ensembl REST API response

        Args:
            session: the aiohttp session for making requests
            variant: variant info as dict

        Returns:
            the same variant, updated with data from the API
        """

        # Call REST API
        vep_info = await self.perform_rest_action(
            session = session,
            endpoint = f"/vep/human/hgvs/{variant['hgvs']}",
            params = {"pick": 1}  # Pick one consequence, chosen according to API's criteria
        )

        # Extract consequences from API response
        transcript_consequences = {}
        try:
            transcript_consequences = vep_info["transcript_consequences"][0]
        except (IndexError, KeyError):
            logger.debug("No transcript consequences found for %s: %s", variant["hgvs"], vep_info)

        # Extract minor allele frequency from API response
        maf = None
        try:
            maf = vep_info["colocated_variants"][-1]["frequencies"][variant["alt"]]["af"]
        except (IndexError, KeyError):
            logger.debug("No minor allele frequency found for %s: %s", variant["hgvs"], vep_info)

        # Update variant with data from API response
        variant.update({
            "gene_id": transcript_consequences.get("gene_id"),
            "gene_symbol": transcript_consequences.get("gene_symbol"),
            "biotype": transcript_consequences.get("biotype"),
            "impact": transcript_consequences.get("impact"),
            "most_severe_consequence": vep_info.get("most_severe_consequence"),
            "minor_allele_frequency": maf
        })

        return variant
