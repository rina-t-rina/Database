"""
Build a non-redundant MySQL database for coronaviruses with a strict limit:
→ MAX 30 genes (CDS features) per virus.

Viruses:
- SARS-CoV-2 (NC_045512.2)
- SARS-CoV   (NC_004718.3)
- MERS-CoV   (NC_019843.3)

Pipeline:
NCBI GenBank CDS
 → RefSeq protein
 → UniProt (batch mapping)
 → UniProt PDB cross-references
 → RCSB metadata + chain mapping
"""

from __future__ import annotations
import os
import time
import json
import hashlib
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple

import requests
import mysql.connector
from Bio import Entrez, SeqIO

# ----------------------------
# User configuration
# ----------------------------
MYSQL_CONFIG = {
    "host": os.getenv("MYSQL_HOST", ""), # enter your hostname
    "user": os.getenv("MYSQL_USER", ""), # enter your username
    "password": os.getenv("MYSQL_PASSWORD", ""), # enter your password
    "database": os.getenv("MYSQL_DB", "corona_nrdb"),
    "port": int(os.getenv("MYSQL_PORT", "")), # enter your port number
}

# NCBI requires you to set an email
Entrez.email = os.getenv("NCBI_EMAIL", "rina.tetaj@students.fhnw.ch")
Entrez.api_key = os.getenv("NCBI_API_KEY")  # optional but recommended

# Reference accessions
VIRUSES = [
    ("SARS-CoV-2", "NC_045512.2"),
    ("SARS-CoV", "NC_004718.3"),
    ("MERS-CoV", "NC_019843.3"),
]

MAX_GENES_PER_VIRUS = 30          # cap number of CDS processed per virus
MAX_PDBS_PER_PROTEIN = 100        # cap PDB entries per protein
NCBI_SLEEP = 0.25                 # polite delay between NCBI calls

# ============================
# DATABASE SCHEMA
# ============================

DDL = [
"""
CREATE TABLE IF NOT EXISTS organism (
    organism_id INT AUTO_INCREMENT PRIMARY KEY,
    scientific_name VARCHAR(255) UNIQUE NOT NULL,
    refseq_accession VARCHAR(32) UNIQUE NOT NULL
) ENGINE=InnoDB;
""",
"""
CREATE TABLE IF NOT EXISTS gene (
    gene_id INT AUTO_INCREMENT PRIMARY KEY,
    organism_id INT NOT NULL,
    gene_symbol VARCHAR(64) NOT NULL,
    gene_name VARCHAR(255),
    start_pos INT NOT NULL,
    end_pos INT NOT NULL,
    strand TINYINT NOT NULL,
    ncbi_geneid BIGINT,
    UNIQUE KEY uq_gene (organism_id, gene_symbol, start_pos, end_pos, strand),
    FOREIGN KEY (organism_id) REFERENCES organism(organism_id)
        ON DELETE CASCADE
) ENGINE=InnoDB;
""",
"""
CREATE TABLE IF NOT EXISTS protein (
    protein_pk INT AUTO_INCREMENT PRIMARY KEY,
    gene_id INT NOT NULL,
    refseq_protein_accession VARCHAR(64),
    uniprot_accession VARCHAR(32),
    product_name VARCHAR(255),
    aa_length INT,
    sequence_sha256 CHAR(64) UNIQUE NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES gene(gene_id)
        ON DELETE CASCADE
) ENGINE=InnoDB;
""",
"""
CREATE TABLE IF NOT EXISTS pdb_entry (
    pdb_id CHAR(4) PRIMARY KEY,
    experimental_method VARCHAR(128),
    resolution_angstrom DECIMAL(6,3),
    rcsb_url VARCHAR(255)
) ENGINE=InnoDB;
""",
"""
CREATE TABLE IF NOT EXISTS protein_structure (
    protein_pk INT NOT NULL,
    pdb_id CHAR(4) NOT NULL,
    entity_id VARCHAR(16),
    asym_ids JSON,
    auth_asym_ids JSON,
    mapping_source VARCHAR(32),
    PRIMARY KEY (protein_pk, pdb_id),
    FOREIGN KEY (protein_pk) REFERENCES protein(protein_pk)
        ON DELETE CASCADE,
    FOREIGN KEY (pdb_id) REFERENCES pdb_entry(pdb_id)
        ON DELETE CASCADE
) ENGINE=InnoDB;
"""
]

# ============================
# CACHES
# ============================

REFSEQ2UNIPROT: Dict[str, List[str]] = {}  # RefSeq protein → UniProt IDs
UNIPROT2PDB: Dict[str, List[str]] = {}     # UniProt → PDB IDs
RCSB_ENTRY_CACHE: Dict[str, Tuple[Optional[str], Optional[float]]] = {}  # PDB → (method, resolution)
RCSB_CHAIN_CACHE: Dict[Tuple[str, str], List[Tuple[str, List[str], List[str]]]] = {}  # (PDB, UniProt) → chain data

# ============================
# HELPERS
# ============================

def sha256_seq(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()  # hash sequence to deduplicate proteins

def mysql_connect():
    print("[DB] Connecting to MySQL")
    conn = mysql.connector.connect(
        host=MYSQL_CONFIG["host"],
        user=MYSQL_CONFIG["user"],
        password=MYSQL_CONFIG["password"],
        port=MYSQL_CONFIG["port"],
    )
    cur = conn.cursor()
    print("[DB] Creating / using database")
    cur.execute(f"CREATE DATABASE IF NOT EXISTS {MYSQL_CONFIG['database']}")
    cur.execute(f"USE {MYSQL_CONFIG['database']}")
    print("[DB] Initializing schema")
    for ddl in DDL:
        cur.execute(ddl)
    cur.close()
    return conn

# ============================
# NCBI PARSING
# ============================

@dataclass
class CDS:
    gene: str              # gene symbol
    product: str           # product/description
    start: int             # 1-based start
    end: int               # end position
    strand: int            # strand (+1/-1/0)
    geneid: Optional[int]  # NCBI GeneID if present
    refseq: Optional[str]  # RefSeq protein accession
    translation: str       # protein sequence

def fetch_genbank(acc: str):
    print(f"[NCBI] Fetching GenBank record {acc}")
    h = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
    time.sleep(NCBI_SLEEP)
    rec = SeqIO.read(h, "genbank")
    h.close()
    print(f"[NCBI] Finished fetching {acc}")
    return rec

def extract_cds(rec) -> List[CDS]:
    print("[NCBI] Extracting CDS features")
    cds = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        q = f.qualifiers
        geneid = None
        for x in q.get("db_xref", []):
            if x.startswith("GeneID:"):
                geneid = int(x.split(":")[1])
        cds.append(CDS(
            gene=q.get("gene", ["unknown"])[0],
            product=q.get("product", [""])[0],
            start=int(f.location.start) + 1,  # convert to 1-based
            end=int(f.location.end),
            strand=f.location.strand or 0,
            geneid=geneid,
            refseq=q.get("protein_id", [None])[0],
            translation=q.get("translation", [""])[0]
        ))
    print(f"[NCBI] Found {len(cds)} CDS features")
    return cds

# ============================
# UNIPROT + PDB
# ============================

UNIPROT_MAP_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_MAP_STATUS = "https://rest.uniprot.org/idmapping/status/{}"
UNIPROT_MAP_RESULTS = "https://rest.uniprot.org/idmapping/results/{}"

def batch_refseq_to_uniprot(refseq_ids: List[str]):
    todo = [r for r in refseq_ids if r and r not in REFSEQ2UNIPROT]  # only unmapped
    if not todo:
        print("[UniProt] All RefSeq IDs already cached")
        return

    print(f"[UniProt] Submitting mapping job for {len(todo)} RefSeq proteins")
    r = requests.post(
        UNIPROT_MAP_RUN,
        data={
            "from": "RefSeq_Protein",
            "to": "UniProtKB",
            "ids": ",".join(todo)
        },
        timeout=30
    )
    r.raise_for_status()
    job_id = r.json()["jobId"]
    print(f"[UniProt] Job ID: {job_id}")

    # brief status polling (informational)
    for i in range(3):
        time.sleep(2)
        try:
            s = requests.get(UNIPROT_MAP_STATUS.format(job_id), timeout=30)
            status = s.json()
            print(f"[UniProt] Status check {i+1}: {status.get('jobStatus')}")
        except Exception:
            print("[UniProt] Status check failed (ignored)")

    # fetch results with retry
    print("[UniProt] Fetching mapping results (retrying if needed)")
    max_tries = 5
    data = None
    for attempt in range(1, max_tries + 1):
        try:
            res = requests.get(UNIPROT_MAP_RESULTS.format(job_id), timeout=60)
            if res.status_code == 200:
                data = res.json()
                print(f"[UniProt] Results retrieved on attempt {attempt}")
                break
            else:
                print(f"[UniProt] Results not ready (HTTP {res.status_code}), retrying...")
        except Exception as e:
            print(f"[UniProt] Results fetch failed ({e}), retrying...")
        time.sleep(2)

    if data is None:
        # if unreachable, record empty mappings
        print("[UniProt] WARNING: mapping results unavailable, assuming no mappings")
        for r in todo:
            REFSEQ2UNIPROT.setdefault(r, [])
        return

    # parse and store mappings
    mapped = set()
    for row in data.get("results", []):
        REFSEQ2UNIPROT.setdefault(row["from"], []).append(row["to"])
        mapped.add(row["from"])
    for r in todo:
        REFSEQ2UNIPROT.setdefault(r, [])
    print(f"[UniProt] Mapping complete ({len(mapped)} mapped / {len(todo)} queried)")

def uniprot_pdbs(uniprot: str) -> List[str]:
    if uniprot in UNIPROT2PDB:
        return UNIPROT2PDB[uniprot]
    print(f"[UniProt] Fetching PDB cross-references for {uniprot}")
    j = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot}.json", timeout=30).json()
    pdbs = sorted({
        x["id"]
        for x in j.get("uniProtKBCrossReferences", [])
        if x.get("database") == "PDB"
    })
    pdbs = pdbs[:MAX_PDBS_PER_PROTEIN]  # enforce cap
    print(f"[UniProt] Using {len(pdbs)} PDB IDs for {uniprot}")
    UNIPROT2PDB[uniprot] = pdbs
    return pdbs

def rcsb_entry(pdb: str):
    if pdb in RCSB_ENTRY_CACHE:
        return RCSB_ENTRY_CACHE[pdb]
    print(f"[RCSB] Fetching entry metadata for {pdb}")
    j = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}").json()
    method = j.get("exptl", [{}])[0].get("method")
    res = j.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
    RCSB_ENTRY_CACHE[pdb] = (method, res)
    return method, res

def rcsb_chains(pdb: str, uniprot: str):
    key = (pdb, uniprot)
    if key in RCSB_CHAIN_CACHE:
        return RCSB_CHAIN_CACHE[key]
    print(f"[RCSB] Resolving chains for {pdb} ↔ {uniprot}")
    entry = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}").json()
    out = []
    for eid in entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]:
        ent = requests.get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}").json()
        refs = ent["rcsb_polymer_entity_container_identifiers"].get("reference_sequence_identifiers", [])
        if any(r["database_accession"] == uniprot for r in refs):
            ids = ent["rcsb_polymer_entity_container_identifiers"]
            out.append((eid, ids.get("asym_ids"), ids.get("auth_asym_ids")))
    print(f"[RCSB] Found {len(out)} matching entities for {pdb}")
    RCSB_CHAIN_CACHE[key] = out
    return out

# ============================
# INGESTION
# ============================

def ingest_virus(conn, name: str, acc: str):
    print(f"\n=== START {name} ===")

    rec = fetch_genbank(acc)
    cur = conn.cursor()

    # insert organism (idempotent)
    print("[DB] Inserting organism")
    cur.execute("""
        INSERT INTO organism (scientific_name, refseq_accession)
        VALUES (%s,%s)
        ON DUPLICATE KEY UPDATE scientific_name=VALUES(scientific_name)
    """, (name, acc))
    cur.execute("SELECT organism_id FROM organism WHERE refseq_accession=%s", (acc,))
    org_id = cur.fetchone()[0]

    # limit CDS count
    cds_list = extract_cds(rec)[:MAX_GENES_PER_VIRUS]
    print(f"[PIPELINE] Using {len(cds_list)} CDS features (limit applied)")

    # batch map all RefSeqs to UniProt
    refseqs = [c.refseq for c in cds_list if c.refseq]
    batch_refseq_to_uniprot(refseqs)

    for i, c in enumerate(cds_list, 1):
        print(f"\n[GENE {i}/{len(cds_list)}] {c.gene}")

        # insert gene (upsert by unique key)
        print("  [DB] Inserting gene")
        cur.execute("""
            INSERT INTO gene (organism_id, gene_symbol, gene_name,
                              start_pos, end_pos, strand, ncbi_geneid)
            VALUES (%s,%s,%s,%s,%s,%s,%s)
            ON DUPLICATE KEY UPDATE gene_name=VALUES(gene_name)
        """, (org_id, c.gene, c.product, c.start, c.end, c.strand, c.geneid))

        cur.execute("""
            SELECT gene_id FROM gene
            WHERE organism_id=%s AND gene_symbol=%s AND start_pos=%s AND end_pos=%s
        """, (org_id, c.gene, c.start, c.end))
        gene_id = cur.fetchone()[0]

        if not c.translation:
            print("  [SKIP] No translation")
            continue

        # insert protein; deduplicate by sequence hash
        print("  [DB] Inserting protein")
        seq_hash = sha256_seq(c.translation)
        cur.execute("""
            INSERT INTO protein (gene_id, refseq_protein_accession,
                                 product_name, aa_length, sequence_sha256)
            VALUES (%s,%s,%s,%s,%s)
            ON DUPLICATE KEY UPDATE gene_id=VALUES(gene_id)
        """, (gene_id, c.refseq, c.product, len(c.translation), seq_hash))

        cur.execute("SELECT protein_pk FROM protein WHERE sequence_sha256=%s", (seq_hash,))
        pk = cur.fetchone()[0]

        # choose first UniProt mapping if available
        uniprots = REFSEQ2UNIPROT.get(c.refseq, [])
        if not uniprots:
            print("  [SKIP] No UniProt mapping")
            continue

        uniprot = uniprots[0]
        print(f"  [UniProt] Using {uniprot}")
        cur.execute("UPDATE protein SET uniprot_accession=%s WHERE protein_pk=%s", (uniprot, pk))

        # attach PDB entries and chain mappings
        pdbs = uniprot_pdbs(uniprot)
        for pdb in pdbs:
            print(f"    [PDB] {pdb}")
            method, res = rcsb_entry(pdb)
            cur.execute("""
                INSERT INTO pdb_entry (pdb_id, experimental_method,
                                       resolution_angstrom, rcsb_url)
                VALUES (%s,%s,%s,%s)
                ON DUPLICATE KEY UPDATE experimental_method=VALUES(experimental_method)
            """, (pdb, method, res, f"https://www.rcsb.org/structure/{pdb}"))

            for eid, asym, auth in rcsb_chains(pdb, uniprot):
                cur.execute("""
                    INSERT IGNORE INTO protein_structure
                    (protein_pk, pdb_id, entity_id, asym_ids, auth_asym_ids, mapping_source)
                    VALUES (%s,%s,%s,%s,%s,'uniprot')
                """, (pk, pdb, eid, json.dumps(asym), json.dumps(auth)))

    conn.commit()
    cur.close()
    print(f"=== DONE {name} ===")

def main():
    conn = mysql_connect()
    for name, acc in VIRUSES:
        ingest_virus(conn, name, acc)
    conn.close()
    print("\nALL DONE")

if __name__ == "__main__":
    main()