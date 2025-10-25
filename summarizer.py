"""
Open ME/CFS Data Pipeline — AI Summarizer + Embedding Generator

Generates:
1. Technical summary (for researchers)
2. Patient-friendly summary (plain language)
3. Embedding vector (numerical representation of meaning)

Author: Jonathan Spruance
Project: Open ME/CFS
"""

# --- Imports ---
import os
import sys
import json
from datetime import date, datetime
from transformers import pipeline
from sentence_transformers import SentenceTransformer
from tqdm import tqdm

# --- Configuration ---
INPUT_PATH = "data/raw_papers_all.json"
OUTPUT_PATH = f"data/mecfs_papers_summarized_{date.today()}.json"

# ⚡ Faster models for dev runs (you can revert later to BART-large)
MODEL_TECHNICAL = "sshleifer/distilbart-cnn-12-6"
MODEL_PATIENT = "sshleifer/distilbart-cnn-12-6"
MODEL_EMBEDDING = "sentence-transformers/all-MiniLM-L6-v2"

# --- Initialize models ---
print("🚀 Loading summarization + embedding models (first run may take a bit)...")
summarizer_technical = pipeline("summarization", model=MODEL_TECHNICAL)
summarizer_patient = pipeline("summarization", model=MODEL_PATIENT)
embedder = SentenceTransformer(MODEL_EMBEDDING)


# --- Utility functions ---
def summarize_text(summarizer, text: str, max_length: int = 180, min_length: int = 60) -> str:
    """Generate a summary from text using the specified summarizer."""
    if not text.strip():
        return ""
    try:
        result = summarizer(
            text,
            max_length=max_length,
            min_length=min_length,
            no_repeat_ngram_size=3,
            clean_up_tokenization_spaces=True,
            do_sample=False,
        )
        return result[0]["summary_text"]
    except Exception as e:
        print(
            f"⚠️  Summarization failed for text starting with: {text[:60]!r}... ({e})")
        return ""


def load_papers(path: str):
    """Load papers from raw_papers.json (handles both old and new schemas)."""
    print(f"📖 Loading abstracts from {path}")
    with open(path, "r", encoding="utf-8") as f:
        raw = json.load(f)

    if isinstance(raw, dict) and "papers" in raw:
        meta = raw.get("metadata", {})
        papers = raw["papers"]
    else:
        meta = {"schema_version": "1.0", "note": "legacy format (no metadata)"}
        papers = raw

    print(
        f"ℹ️  Loaded {len(papers)} papers (schema {meta.get('schema_version', 'unknown')})")
    return papers, meta


# --- Main function ---
def generate_summaries(quick: bool = False):
    """Reads abstracts, generates two summaries and an embedding per paper, saves results."""
    papers, raw_metadata = load_papers(INPUT_PATH)

    # --- Quick mode for development ---
    if quick:
        papers = papers[:100]
        print(
            f"⚡ Quick mode active: processing first {len(papers)} papers only.")

    # --- Resume support ---
    summarized = []
    done_pmids = set()

    if os.path.exists(OUTPUT_PATH):
        print(f"⏸️  Found existing output file: {OUTPUT_PATH}")
        try:
            with open(OUTPUT_PATH, "r", encoding="utf-8") as f:
                summarized = json.load(f)
            done_pmids = {str(p.get("pmid"))
                          for p in summarized if p.get("pmid")}
            print(
                f"🔁 Resuming — {len(done_pmids)} already summarized, {len(papers) - len(done_pmids)} remaining.")
        except Exception as e:
            print(f"⚠️  Could not load existing file, starting fresh: {e}")

    now = datetime.now().isoformat(timespec="seconds")

    print(
        f"🤖 Generating outputs using:\n"
        f"  • Technical: {MODEL_TECHNICAL}\n"
        f"  • Patient:   {MODEL_PATIENT}\n"
        f"  • Embedding: {MODEL_EMBEDDING}"
    )

    for paper in tqdm(papers, desc="Summarizing & embedding"):
        pmid = str(paper.get("pmid"))
        if pmid in done_pmids:
            continue

        abstract = paper.get("abstract", "")
        if not abstract:
            continue

        # --- Generate summaries ---
        tech_summary = summarize_text(summarizer_technical, abstract)
        prompt = f"Explain this research in simple, clear language for patients: {abstract}"
        patient_summary = summarize_text(
            summarizer_patient, prompt, max_length=130, min_length=40)

        # --- Generate embedding (from technical summary if available, else abstract) ---
        embedding = None
        try:
            text_for_embedding = tech_summary or abstract
            embedding = embedder.encode(text_for_embedding).tolist()
        except Exception as e:
            print(f"⚠️  Embedding generation failed for PMID {pmid}: {e}")

        summarized.append({
            "pmid": pmid,
            "title": paper.get("title"),
            "abstract": abstract,
            "authors": paper.get("authors", []),
            "year": paper.get("year"),
            "technical_summary": tech_summary,
            "patient_summary": patient_summary,
            "embedding": embedding,
            "metadata": {
                "technical_summary_model": MODEL_TECHNICAL,
                "patient_summary_model": MODEL_PATIENT,
                "embedding_model": MODEL_EMBEDDING,
                "summarized_at": now,
                "fetched_at": paper.get("fetched_at"),
                "raw_schema_version": raw_metadata.get("schema_version", "unknown"),
                "source": raw_metadata.get("source", "unknown"),
            },
        })

        # --- Periodically save progress every 25 papers ---
        if len(summarized) % 25 == 0:
            with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
                json.dump(summarized, f, indent=2, ensure_ascii=False)

    # --- Final save ---
    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(summarized, f, indent=2, ensure_ascii=False)

    print(
        f"✅ Saved {len(summarized)} summarized papers with embeddings to {OUTPUT_PATH}")


# --- Entrypoint ---
if __name__ == "__main__":
    quick_flag = "--quick" in sys.argv
    generate_summaries(quick=quick_flag)
