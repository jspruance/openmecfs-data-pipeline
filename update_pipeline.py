from fetch_pubmed import fetch_pubmed_all
from summarizer import generate_summaries


def main():
    print("🚀 Running full Open ME/CFS data pipeline")
    fetch_pubmed_all()          # fetch all papers
    generate_summaries()        # summarize + embed
    print("✅ Pipeline complete")
