from flask import Flask, render_template, request
import mysql.connector
from db import get_connection
from flask import jsonify

app = Flask(__name__)

VIRUSES = ["SARS-CoV-2", "SARS-CoV", "MERS-CoV"]  # supported virus names for UI

# ============================================================
# HOME
# ============================================================

@app.route("/")
def index():
    breadcrumbs = [{"label": "Home", "url": "/"}]
    # Render landing page listing viruses
    return render_template(
        "index.html",
        viruses=VIRUSES,
        breadcrumbs=breadcrumbs
    )

# ============================================================
# HTML PAGES
# ============================================================

@app.route("/genes")
def genes():
    virus = request.args.get("virus")

    # Query genes for the selected virus
    conn = get_connection()
    cur = conn.cursor(dictionary=True)
    cur.execute("""
        SELECT
            g.gene_symbol,
            g.gene_name,
            g.start_pos,
            g.end_pos,
            g.strand,
            g.ncbi_geneid
        FROM gene g
        JOIN organism o ON g.organism_id = o.organism_id
        WHERE o.scientific_name = %s
        ORDER BY g.start_pos
    """, (virus,))
    genes = cur.fetchall()
    cur.close()
    conn.close()

    breadcrumbs = [
        {"label": "Home", "url": "/"},
        {"label": virus, "url": "/"},
        {"label": "Genes", "url": f"/genes?virus={virus}"}
    ]

    # Render genes page
    return render_template(
        "genes.html",
        virus=virus,
        genes=genes,
        breadcrumbs=breadcrumbs
    )


@app.route("/proteins")
def proteins():
    virus = request.args.get("virus")

    # Query proteins for the selected virus with aggregated PDB info
    conn = get_connection()
    cur = conn.cursor(dictionary=True)
    cur.execute("""
        SELECT
            p.protein_pk,
            p.product_name,
            p.uniprot_accession,
            p.refseq_protein_accession,
            p.aa_length,
            GROUP_CONCAT(DISTINCT ps.pdb_id ORDER BY ps.pdb_id SEPARATOR ',') AS pdb_ids,
            SUM(pe.experimental_method LIKE '%cryo%') AS cryo_em,
            SUM(pe.experimental_method LIKE '%X-ray%') AS xray,
            SUM(pe.experimental_method LIKE '%NMR%') AS nmr,
            SUM(
                pe.experimental_method IS NOT NULL
                AND pe.experimental_method NOT LIKE '%cryo%'
                AND pe.experimental_method NOT LIKE '%X-ray%'
                AND pe.experimental_method NOT LIKE '%NMR%'
            ) AS other_methods
        FROM protein p
        JOIN gene g ON p.gene_id = g.gene_id
        JOIN organism o ON g.organism_id = o.organism_id
        LEFT JOIN protein_structure ps ON ps.protein_pk = p.protein_pk
        LEFT JOIN pdb_entry pe ON ps.pdb_id = pe.pdb_id
        WHERE o.scientific_name = %s
        GROUP BY p.protein_pk
    """, (virus,))
    proteins = cur.fetchall()
    cur.close()
    conn.close()

    breadcrumbs = [
        {"label": "Home", "url": "/"},
        {"label": virus, "url": "/"},
        {"label": "Proteins", "url": f"/proteins?virus={virus}"}
    ]

    # Render proteins page
    return render_template(
        "proteins.html",
        virus=virus,
        proteins=proteins,
        breadcrumbs=breadcrumbs
    )


@app.route("/pdb/<pdb_id>")
def pdb_page(pdb_id):
    # Fetch one PDB entry and its associated protein/virus
    conn = get_connection()
    cur = conn.cursor(dictionary=True)

    cur.execute("""
        SELECT
            pe.pdb_id,
            pe.experimental_method,
            pe.resolution_angstrom,
            pe.rcsb_url,
            p.product_name,
            p.uniprot_accession,
            o.scientific_name AS virus
        FROM pdb_entry pe
        JOIN protein_structure ps ON pe.pdb_id = ps.pdb_id
        JOIN protein p ON ps.protein_pk = p.protein_pk
        JOIN gene g ON p.gene_id = g.gene_id
        JOIN organism o ON g.organism_id = o.organism_id
        WHERE pe.pdb_id = %s
        LIMIT 1
    """, (pdb_id.upper(),))

    pdb = cur.fetchone()
    cur.close()
    conn.close()

    breadcrumbs = [
        {"label": "Home", "url": "/"},
        {"label": pdb["virus"], "url": "/"},
        {"label": "Proteins", "url": f"/proteins?virus={pdb['virus']}"},
        {"label": pdb["pdb_id"], "url": f"/pdb/{pdb['pdb_id']}"}
    ]

    # Render PDB details page
    return render_template(
        "pdb.html",
        pdb=pdb,
        breadcrumbs=breadcrumbs
    )

# ============================================================
# API ENDPOINTS (JSON ONLY)
# ============================================================

@app.route("/api/genes")
def api_genes():
    virus = request.args.get("virus")
    if not virus:
        return jsonify({"error": "virus parameter required"}), 400

    # Return genes for a virus as JSON
    conn = get_connection()
    cur = conn.cursor(dictionary=True)
    cur.execute("""
        SELECT
            g.gene_symbol,
            g.gene_name,
            g.start_pos,
            g.end_pos,
            g.strand,
            g.ncbi_geneid
        FROM gene g
        JOIN organism o ON g.organism_id = o.organism_id
        WHERE o.scientific_name = %s
        ORDER BY g.start_pos
    """, (virus,))
    data = cur.fetchall()
    cur.close()
    conn.close()

    return jsonify({
        "virus": virus,
        "count": len(data),
        "genes": data
    })


@app.route("/api/proteins")
def api_proteins():
    virus = request.args.get("virus")
    if not virus:
        return jsonify({"error": "virus parameter required"}), 400

    # Return proteins for a virus as JSON
    conn = get_connection()
    cur = conn.cursor(dictionary=True)
    cur.execute("""
        SELECT
            p.product_name,
            p.uniprot_accession,
            p.refseq_protein_accession,
            p.aa_length
        FROM protein p
        JOIN gene g ON p.gene_id = g.gene_id
        JOIN organism o ON g.organism_id = o.organism_id
        WHERE o.scientific_name = %s
    """, (virus,))
    data = cur.fetchall()
    cur.close()
    conn.close()

    return jsonify({
        "virus": virus,
        "count": len(data),
        "proteins": data
    })


@app.route("/api/pdbs")
def api_pdbs():
    virus = request.args.get("virus")
    if not virus:
        return jsonify({"error": "virus parameter required"}), 400

    # Return distinct PDB entries for a virus as JSON
    conn = get_connection()
    cur = conn.cursor(dictionary=True)
    cur.execute("""
        SELECT DISTINCT
            pe.pdb_id,
            pe.experimental_method,
            pe.resolution_angstrom,
            pe.rcsb_url
        FROM pdb_entry pe
        JOIN protein_structure ps ON pe.pdb_id = ps.pdb_id
        JOIN protein p ON ps.protein_pk = p.protein_pk
        JOIN gene g ON p.gene_id = g.gene_id
        JOIN organism o ON g.organism_id = o.organism_id
        WHERE o.scientific_name = %s
    """, (virus,))
    data = cur.fetchall()
    cur.close()
    conn.close()

    return jsonify({
        "virus": virus,
        "count": len(data),
        "pdbs": data
    })


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    app.run(debug=True)  # Dev server; enable debug for auto-reload