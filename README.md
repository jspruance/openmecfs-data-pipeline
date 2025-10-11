# 🧠 Open ME/CFS Data Pipeline

The **Open ME/CFS Data Pipeline** is the first stage of the [Open ME/CFS](https://github.com/yourusername/openmecfs-platform) initiative — an open-source effort to accelerate ME/CFS research through data aggregation and artificial intelligence.

This pipeline automatically:

1. Fetches research papers on **Myalgic Encephalomyelitis / Chronic Fatigue Syndrome (ME/CFS)** from **PubMed** using the NCBI Entrez API.
2. Extracts key metadata (titles, abstracts, authors, PMIDs).
3. Summarizes each abstract using **Hugging Face Transformer models** such as BioBART and T5.
4. Outputs clean, structured JSON files ready for use in the Open ME/CFS API and web platform.

---

## 🚀 Features

- 🔍 Automatic retrieval of ME/CFS studies from PubMed
- 🧠 AI-powered summarization (technical + plain-language summaries)
- 💾 Structured JSON outputs for downstream analysis or database import
- ⚙️ Modular, lightweight Python codebase
- 🧩 Extensible — future integrations planned for ClinicalTrials.gov, NIH RePORTER, and PubMed Central full-text datasets

---

## 🧰 Tech Stack

| Layer            | Tool / Library                                       | Purpose                                 |
| ---------------- | ---------------------------------------------------- | --------------------------------------- |
| Data Access      | [Biopython](https://biopython.org/)                  | Interface to PubMed / NCBI Entrez API   |
| AI Summarization | [Hugging Face Transformers](https://huggingface.co/) | Pretrained models (BioBART, T5)         |
| Compute Engine   | [PyTorch](https://pytorch.org/)                      | Neural network backend for Transformers |
| API-ready        | [FastAPI](https://fastapi.tiangolo.com/)             | Future integration for serving data     |
| Storage          | JSON / (Supabase planned)                            | Structured, open-format outputs         |

---

## 🧩 Repository Purpose

This repo powers the **knowledge ingestion layer** of the Open ME/CFS ecosystem.  
It serves as the foundation for:

- the [`openmecfs-platform`](https://github.com/yourusername/openmecfs-platform) (API and database),
- and the [`openmecfs-ui`](https://github.com/yourusername/openmecfs-ui) (public-facing research explorer).

---

## 🧭 Roadmap

| Phase          | Goal                                                          |
| -------------- | ------------------------------------------------------------- |
| ✅ **Phase 1** | Fetch + summarize PubMed abstracts                            |
| 🧩 **Phase 2** | Integrate with ClinicalTrials.gov + NIH RePORTER              |
| ⚙️ **Phase 3** | Expand to PubMed Central full texts                           |
| 🌐 **Phase 4** | Deploy REST API (FastAPI + Supabase)                          |
| 💬 **Phase 5** | Support a public ME/CFS research dashboard & community portal |

---

## 💙 Mission

The goal of Open ME/CFS is to **democratize biomedical research** — using open data, open source, and AI to connect scientists, clinicians, and patients working toward understanding and treating Myalgic Encephalomyelitis / Chronic Fatigue Syndrome.

---

## 👤 Author

**Jonathan Spruance**  
[GitHub](https://github.com/jspruance) · [openmecfs.org](https://openmecfs.org)

---

## 🪪 License

MIT License © 2025 Jonathan Spruance
