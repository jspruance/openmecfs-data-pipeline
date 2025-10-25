# 🔬 Open ME/CFS — Data Pipeline

**Automated biomedical data ingestion and preprocessing pipeline**  
Part of the [Open ME/CFS](https://github.com/jspruance) ecosystem.

---

## 🧩 Overview

`openmecfs-data-pipeline` automates **collection, cleaning, summarization, and embedding** of ME/CFS-related biomedical research data.

It is designed to continuously ingest data from sources like **PubMed**, **NIH RePORTER**, and **DecodeME**, preparing it for downstream analysis in:

- [`openmecfs-platform`](https://github.com/jspruance/openmecfs-platform) — FastAPI backend & Supabase integration
- [`openmecfs-ai-cure`](https://github.com/jspruance/openmecfs-ai-cure) — clustering & mechanistic hypothesis generation
- [`openmecfs-ui`](https://github.com/jspruance/openmecfs-ui) — research explorer frontend

---

## ⚙️ Pipeline Stages

| Stage                      | Description                                                       | Status         |
| -------------------------- | ----------------------------------------------------------------- | -------------- |
| **1. Fetch papers**        | Retrieve ME/CFS-related abstracts and metadata from PubMed API.   | ✅ Done        |
| **2. Clean + dedupe**      | Normalize author names, titles, and remove duplicates.            | ✅ Done        |
| **3. Summarize**           | Generate technical and lay summaries using Hugging Face models.   | ✅ Done        |
| **4. Embeddings**          | Create semantic embeddings with `sentence-transformers`.          | ✅ Done        |
| **5. JSON export**         | Save enriched data to `/data/processed/` in JSON format.          | ✅ Done        |
| **6. Upload to Supabase**  | Push records and embeddings to Postgres (`papers` table).         | ✅ Done        |
| **7. Incremental updates** | Only fetch and update new or modified papers.                     | ⚙️ In progress |
| **8. Integration**         | Provide data to `openmecfs-ai-cure` for subtyping and clustering. | ✅ Connected   |

---

## 🧰 Tech Stack

| Layer           | Tools                                                                    |
| --------------- | ------------------------------------------------------------------------ |
| Language        | Python 3.12                                                              |
| Package Manager | Poetry                                                                   |
| APIs            | PubMed E-Utilities, NIH RePORTER                                         |
| NLP             | Hugging Face Transformers (`facebook/bart-large-cnn`, `allenai/scibert`) |
| Embeddings      | `sentence-transformers/all-MiniLM-L6-v2`                                 |
| Storage         | Supabase (Postgres + pgvector)                                           |
| Logging         | Rich + tqdm                                                              |

---

## 🛠️ Setup

### 1️⃣ Clone the repository

```bash
git clone https://github.com/jspruance/openmecfs-data-pipeline.git
cd openmecfs-data-pipeline
```

### 2️⃣ Install dependencies

```bash
poetry install
```

### 3️⃣ Create `.env` file

```bash
SUPABASE_URL=https://YOUR_PROJECT.supabase.co
SUPABASE_SERVICE_ROLE_KEY=YOUR_SERVICE_ROLE_KEY
SUPABASE_ANON_KEY=YOUR_ANON_KEY
```

---

## 🚀 Usage

### Fetch new PubMed papers

```bash
poetry run python -m pipeline.fetch_papers --query "myalgic encephalomyelitis chronic fatigue syndrome"
```

### Summarize abstracts

```bash
poetry run python -m pipeline.summarize_abstracts
```

### Generate embeddings

```bash
poetry run python -m pipeline.embed_papers
```

### Upload to Supabase

```bash
poetry run python -m pipeline.json_to_db
```

---

## 📦 Folder Structure

```
openmecfs-data-pipeline/
├── pipeline/
│   ├── fetch_papers.py        # Pulls PubMed results
│   ├── summarize_abstracts.py # BART-based summarization
│   ├── embed_papers.py        # Creates MiniLM embeddings
│   ├── json_to_db.py          # Upload to Supabase
│   └── utils.py               # Shared helper functions
├── data/
│   ├── raw/                   # Raw PubMed JSONs
│   ├── processed/             # Summarized + embedded
├── .env                       # Supabase credentials
├── pyproject.toml
└── README.md
```

---

## 🧠 How It Fits Into the Ecosystem

```
PubMed / NIH / DecodeME
       ↓
openmecfs-data-pipeline (fetch → summarize → embed → upload)
       ↓
openmecfs-platform (API + Supabase)
       ↓
openmecfs-ai-cure (clustering + hypothesis generation)
       ↓
openmecfs-ui (interactive research explorer)
```

---

## 👤 Author

**Jonathan Spruance**  
[GitHub](https://github.com/jspruance) · [openmecfs.org](https://openmecfs.org)

---

## 📚 License

MIT License © 2025 Jonathan Spruance
