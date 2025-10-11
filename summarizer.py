"""
Open ME/CFS Data Pipeline — AI Summarizer

Uses a pretrained model from Hugging Face (BioBART, T5, or BART)
to generate readable summaries of ME/CFS research abstracts.

Outputs two summary types:
1. Technical summary (for researchers)
2. Patient-friendly summary (plain language)

Author: Jonathan Spruance
Project: Open ME/CFS
"""

# --- Imports ---
from datetime import date
import json
from datetime import datetime
from transformers import pipeline
from tqdm import tqdm

# --- Configuration ---
INPUT_PATH = "data/raw_papers.json"
OUTPUT_PATH = f"data/mecfs_papers_summarized_{date.today()}.json"


# You can experiment with different models:
MODEL_TECHNICAL = "facebook/bart-large-cnn"   # for technical summary
MODEL_PATIENT = "allenai/biobart-v2"          # for plain-language summary

# Initialize summarizers
summarizer_technical = pipeline("summarization", model=MODEL_TECHNICAL)
summarizer_patient = pipeline("summarization", model=MODEL_PATIENT)


# --- Functions ---
def summarize_text(summarizer, text: str, max_length: int = 180, min_length: int = 60) -> str:
    """Generate a summary from text using the specified summarizer."""
    if not text.strip():
        return ""
    result = summarizer(text, max_length=max_length,
                        min_length=min_length, do_sample=False)
    return result[0]["summary_text"]


def generate_summaries():
    """Reads abstracts, generates two summaries per paper, and saves results."""
    print(f"📖 Loading abstracts from {INPUT_PATH}")
    with open(INPUT_PATH, "r", encoding="utf-8") as f:
        papers = json.load(f)

    summarized = []
    now = datetime.utcnow().isoformat()

    print(
        f"🤖 Generating summaries using:\n  • Technical: {MODEL_TECHNICAL}\n  • Patient: {MODEL_PATIENT}")

    for paper in tqdm(papers, desc="Summarizing papers"):
        abstract = paper.get("abstract", "")
        if not abstract:
            continue

        # --- Technical summary ---
        tech_summary = summarize_text(summarizer_technical, abstract)

        # --- Patient-friendly summary ---
        prompt = f"Explain this research in simple, clear language for patients: {abstract}"
        patient_summary = summarize_text(
            summarizer_patient, prompt, max_length=130, min_length=40)

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
                "fetched_at": paper.get("fetched_at")
            }
        })

    with open(OUTPUT_PATH, "w", encoding="utf-8") as f:
        json.dump(summarized, f, indent=2, ensure_ascii=False)

    print(f"✅ Saved {len(summarized)} summarized papers to {OUTPUT_PATH}")


if __name__ == "__main__":
    generate_summaries()
