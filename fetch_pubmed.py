"""
Open ME/CFS Data Pipeline — PubMed Ingest (Paginated)
Fetches *all* ME/CFS-related papers from PubMed in batches using Biopython’s Entrez API.
"""

from datetime import datetime
import json
import re
from typing import List, Dict
from Bio import Entrez
import time

# Identify yourself to NCBI
Entrez.email = "jgspruance@gmail.com"

SEARCH_TERM = (
    '"myalgic encephalomyelitis"'
    ' OR "chronic fatigue syndrome"'
    ' OR "ME/CFS"'
    ' OR "MECFS"'
)
BATCH_SIZE = 500  # number of records to fetch per batch
SLEEP_BETWEEN_CALLS = 0.5  # seconds (respect NCBI rate limits)


def extract_pub_year(article) -> int | None:
    try:
        journal_info = article["MedlineCitation"]["Article"]["Journal"]["JournalIssue"]["PubDate"]
        if "Year" in journal_info:
            return int(journal_info["Year"])
        elif "MedlineDate" in journal_info:
            match = re.search(r"(19|20)\d{2}", str(
                journal_info["MedlineDate"]))
            if match:
                return int(match.group(0))
    except Exception:
        pass
    return None


def fetch_pubmed_all() -> None:
    """Fetch all ME/CFS papers from PubMed using batched Entrez queries."""
    print("🔍 Searching PubMed for ME/CFS papers...")
    search_handle = Entrez.esearch(db="pubmed", term=SEARCH_TERM, retmax=1)
    search_record = Entrez.read(search_handle)
    total_count = int(search_record["Count"])
    print(f"📚 Total papers found: {total_count}")

    all_results: List[Dict] = []
    fetched = 0

    while fetched < total_count:
        print(
            f"➡️  Fetching records {fetched + 1}–{min(fetched + BATCH_SIZE, total_count)}...")
        handle = Entrez.esearch(
            db="pubmed",
            term=SEARCH_TERM,
            retstart=fetched,
            retmax=BATCH_SIZE,
        )
        record = Entrez.read(handle)
        ids = record["IdList"]

        if not ids:
            break

        fetch = Entrez.efetch(db="pubmed", id=",".join(
            ids), rettype="abstract", retmode="xml")
        papers = Entrez.read(fetch)

        for paper in papers["PubmedArticle"]:
            citation = paper["MedlineCitation"]["Article"]
            title = citation.get("ArticleTitle", "")
            abstract = " ".join(citation.get(
                "Abstract", {}).get("AbstractText", []))
            authors = [a.get("LastName", "")
                       for a in citation.get("AuthorList", [])]
            pmid = paper["MedlineCitation"]["PMID"]
            year = extract_pub_year(paper)

            all_results.append({
                "pmid": str(pmid),
                "title": title,
                "abstract": abstract,
                "authors": authors,
                "year": year,
                "fetched_at": datetime.utcnow().isoformat(),
            })

        fetched += len(ids)
        time.sleep(SLEEP_BETWEEN_CALLS)  # rate-limit safety

    # Wrap results with metadata
    data_out = {
        "metadata": {
            "search_term": SEARCH_TERM,
            "fetched_at": datetime.utcnow().isoformat(),
            "total_results": len(all_results),
            "source": "PubMed (NCBI Entrez)",
        },
        "papers": all_results,
    }

    output_path = "data/raw_papers_all.json"
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data_out, f, indent=2, ensure_ascii=False)

    print(f"✅ Saved {len(all_results)} papers to {output_path}")


if __name__ == "__main__":
    fetch_pubmed_all()
