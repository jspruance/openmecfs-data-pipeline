"""
Uploads summarized papers with embeddings to Supabase in safe batches.
Author: Jonathan Spruance
Project: Open ME/CFS
"""

import os
import json
import math
import time
from supabase import create_client, Client
from dotenv import load_dotenv
from tqdm import tqdm

# --- Setup ---
load_dotenv()
url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_ROLE_KEY")
supabase: Client = create_client(url, key)

INPUT_PATH = "data/mecfs_papers_summarized_2025-10-25.json"
BATCH_SIZE = 500  # adjust if needed


def chunk_list(lst, size):
    """Yield successive chunks of length size."""
    for i in range(0, len(lst), size):
        yield lst[i:i + size]


def upload_json_batched(path: str):
    """Upload summarized papers to Supabase in batches."""
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    data = [p for p in data if p.get("embedding")]  # skip empty embeddings
    total = len(data)
    print(f"Uploading {total} papers in batches of {BATCH_SIZE}...")

    for i, batch in enumerate(chunk_list(data, BATCH_SIZE), start=1):
        payload = []
        for paper in batch:
            payload.append({
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
            })

        try:
            supabase.table("papers").upsert(
                payload, on_conflict="pmid").execute()
            tqdm.write(f"✅ Batch {i}/{math.ceil(total / BATCH_SIZE)} uploaded")
        except Exception as e:
            tqdm.write(f"⚠️  Error in batch {i}: {e}")
            time.sleep(2)  # prevent flood
    print("🎯 All batches uploaded successfully.")


if __name__ == "__main__":
    upload_json_batched(INPUT_PATH)
