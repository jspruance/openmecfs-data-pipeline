"""
Open ME/CFS Data Pipeline — PubMed Ingest

Fetches ME/CFS-related research papers from PubMed using Biopython’s Entrez API,
extracts key metadata (title, authors, abstract, PMID, year), and saves the results
to a local JSON file for later AI summarization.

Author: Jonathan Spruance
Project: Open ME/CFS
"""

from datetime import datetime
import json
import re
from typing import List, Dict
from Bio import Entrez

# Required by NCBI Entrez — identifies you as the requesting user
Entrez.email = "jgspruance@gmail.com"  # type: ignore

# Search parameters
SEARCH_TERM = (
    '"myalgic encephalomyelitis"'
    ' OR "chronic fatigue syndrome"'
    ' OR "ME/CFS"'
    ' OR "MECFS"'
)
MAX_RESULTS = 100


def extract_pub_year(article) -> int | None:
    """Extracts publication year from PubMed article metadata."""
    try:
        journal_info = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]

        # Prefer explicit <Year> tag
        if "Year" in journal_info:
            return int(journal_info["Year"])

        # Fallback: look for 4-digit year in MedlineDate
        elif "MedlineDate" in journal_info:
            match = re.search(r"(19|20)\d{2}", str(
                journal_info["MedlineDate"]))
            if match:
                return int(match.group(0))
    except Exception:
        pass
    return None


def fetch_pubmed() -> None:
    """Fetches ME/CFS papers from PubMed using Entrez API."""
    print("🔍 Fetching papers from PubMed...")

    # Step 1: Search PubMed
    handle = Entrez.esearch(db="pubmed", term=SEARCH_TERM, retmax=MAX_RESULTS)
    record = Entrez.read(handle)
    ids: List[str] = record["IdList"]

    print(f"Found {len(ids)} paper IDs. Fetching abstracts...")

    # Step 2: Fetch paper details
    fetch = Entrez.efetch(db="pubmed", id=",".join(
        ids), rettype="abstract", retmode="xml")
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
        year = extract_pub_year(paper)  # ✅ Added publication year extraction

        results.append({
            "pmid": str(pmid),
            "title": title,
            "abstract": abstract,
            "authors": authors,
            "year": year,
            "fetched_at": datetime.utcnow().isoformat(),
        })

    # Step 4: Wrap results with metadata header
    data_out = {
        "metadata": {
            "search_term": SEARCH_TERM,
            "fetched_at": datetime.utcnow().isoformat(),
            "total_results": len(results),
            "source": "PubMed (NCBI Entrez)",
        },
        "papers": results,
    }

    # Step 5: Save to JSON
    output_path = "data/raw_papers.json"
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data_out, f, indent=2, ensure_ascii=False)

    print(f"✅ Saved {len(results)} papers to {output_path}")


if __name__ == "__main__":
    fetch_pubmed()
