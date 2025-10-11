"""
Open ME/CFS Data Pipeline — PubMed Ingest

Fetches ME/CFS-related research papers from PubMed using Biopython’s Entrez API,
extracts key metadata (title, authors, abstract, PMID), and saves the results
to a local JSON file for later AI summarization.

Author: Jonathan Spruance
Project: Open ME/CFS
"""

from datetime import datetime
import json
from typing import List, Dict
from Bio import Entrez

# Required by NCBI Entrez — this identifies you as the requesting user
Entrez.email = "jgspruance@gmail.com"  # type: ignore

# Search parameters
SEARCH_TERM = (
    '"myalgic encephalomyelitis"'
    ' OR "chronic fatigue syndrome"'
    ' OR "ME/CFS"'
    ' OR "MECFS"'
)

MAX_RESULTS = 100


def fetch_pubmed() -> None:
    """
    Fetches ME/CFS papers from PubMed using Entrez API.

    1. Searches PubMed for the specified SEARCH_TERM.
    2. Retrieves up to MAX_RESULTS abstracts in XML format.
    3. Extracts metadata fields and writes to data/raw_papers.json.

    Returns:
        None
    """
    print("🔍 Fetching papers from PubMed...")

    # Step 1: Search PubMed
    handle = Entrez.esearch(db="pubmed", term=SEARCH_TERM, retmax=MAX_RESULTS)
    record = Entrez.read(handle)
    ids: List[str] = record["IdList"]

    print(f"Found {len(ids)} paper IDs. Fetching abstracts...")

    # Step 2: Fetch paper details
    fetch = Entrez.efetch(db="pubmed", id=",".join(ids),
                          rettype="abstract", retmode="xml")
    papers = Entrez.read(fetch)

    results: List[Dict] = []

    # Step 3: Extract data from XML
    for paper in papers["PubmedArticle"]:
        citation = paper["MedlineCitation"]["Article"]
        title = citation.get("ArticleTitle", "")
        abstract = " ".join(citation.get(
            "Abstract", {}).get("AbstractText", []))
        authors = [a.get("LastName", "")
                   for a in citation.get("AuthorList", [])]
        pmid = paper["MedlineCitation"]["PMID"]

        results.append({
            "pmid": str(pmid),
            "title": title,
            "abstract": abstract,
            "authors": authors,
            "fetched_at": datetime.utcnow().isoformat()
        })

    # Step 4: Save results to local JSON
    output_path = "data/raw_papers.json"
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"✅ Saved {len(results)} papers to {output_path}")


if __name__ == "__main__":
    fetch_pubmed()
