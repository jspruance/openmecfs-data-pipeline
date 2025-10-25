"""
Uploads summarized papers with embeddings to Supabase.
"""

import os
import json
from supabase import create_client, Client
from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()
url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_ROLE_KEY")
supabase: Client = create_client(url, key)

INPUT_PATH = "data/mecfs_papers_summarized_2025-10-25.json"


def upload_json(path: str):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    print(f"Uploading {len(data)} papers to Supabase...")
    for paper in tqdm(data):
        try:
            payload = {
                "pmid": paper.get("pmid"),
                "title": paper.get("title"),
                "abstract": paper.get("abstract"),
                "authors": paper.get("authors"),
                "year": paper.get("year"),
                "technical_summary": paper.get("technical_summary"),
                "patient_summary": paper.get("patient_summary"),
                "embedding": paper.get("embedding"),
                "fetched_at": paper.get("metadata", {}).get("fetched_at"),
                "summarized_at": paper.get("metadata", {}).get("summarized_at"),
            }
            supabase.table("papers").upsert(
                payload, on_conflict="pmid").execute()
        except Exception as e:
            print("⚠️  Error uploading", paper.get("pmid"), e)
    print("✅ Upload complete.")


if __name__ == "__main__":
    upload_json(INPUT_PATH)
