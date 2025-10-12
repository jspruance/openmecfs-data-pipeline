"""
Open ME/CFS Data Pipeline — AI Summarizer

Uses pretrained Hugging Face models (BioBART, T5, or BART)
to generate readable summaries of ME/CFS research abstracts.

Outputs two summary types:
1. Technical summary (for researchers)
2. Patient-friendly summary (plain language)

Author: Jonathan Spruance
Project: Open ME/CFS
"""

# --- Imports ---
import json
from datetime import date, datetime
from transformers import pipeline
from tqdm import tqdm

# --- Configuration ---
INPUT_PATH = "data/raw_papers.json"
OUTPUT_PATH = f"data/mecfs_papers_summarized_{date.today()}.json"

# You can experiment with different models:
MODEL_TECHNICAL = "philschmid/bart-large-cnn-samsum"  # for technical summaries
# for patient-friendly summaries
MODEL_PATIENT = "facebook/bart-large-cnn"

# Initialize summarizers
summarizer_technical = pipeline("summarization", model=MODEL_TECHNICAL)
summarizer_patient = pipeline("summarization", model=MODEL_PATIENT)


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

def generate_summaries():
    """Reads abstracts, generates two summaries per paper, and saves results."""
    papers, raw_metadata = load_papers(INPUT_PATH)
    summarized = []
    now = datetime.now().isoformat(timespec="seconds")

    print(
        f"🤖 Generating summaries using:\n  • Technical: {MODEL_TECHNICAL}\n  • Patient:   {MODEL_PATIENT}")

    for paper in tqdm(papers, desc="Summarizing papers"):
        if not isinstance(paper, dict):
            print("⚠️  Skipping malformed entry (not a dict)")
            continue

        abstract = paper.get("abstract", "")
        if not abstract:
            continue

        # --- Technical summary ---
        tech_summary = summarize_text(summarizer_technical, abstract)

        # --- Patient-friendly summary ---
        prompt = f"Explain this research in simple, clear language for patients: {abstract}"
        patient_summary = summarize_text(
            summarizer_patient, prompt, max_length=130, min_length=40
        )

        summarized.append({
            "pmid": paper.get("pmid"),
            "title": paper.get("title"),
            "abstract": abstract,
            "authors": paper.get("authors", []),
            "technical_summary": tech_summary,
            "patient_summary": patient_summary,
            "metadata": {
                "technical_summary_model": MODEL_TECHNICAL,
                "patient_summary_model": MODEL_PATIENT,
                "summarized_at": now,
                "fetched_at": paper.get("fetched_at"),
                "raw_schema_version": raw_metadata.get("schema_version", "unknown"),
                "source": raw_metadata.get("source", "unknown")
            }
        })

    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(summarized, f, indent=2, ensure_ascii=False)

    print(f"✅ Saved {len(summarized)} summarized papers to {OUTPUT_PATH}")


# --- Entrypoint ---
if __name__ == "__main__":
    generate_summaries()
