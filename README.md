# Coronavirus Genome & Protein Structure Explorer
This project is a web-based application and database pipeline for exploring severe coronavirus genomes, genes, proteins, and their 3D structures.
It integrates data from NCBI GenBank, UniProt, and the Protein Data Bank (PDB) into a local MySQL database, which can be explored via a simple Flask web interface.

# Features
## Virus-Based Exploration
## Browse genes and proteins for:
	•	SARS-CoV-2
	•	SARS-CoV
	•	MERS-CoV
## Automated Data Integration (ETL)
	•	Genome and gene extraction from NCBI GenBank
	•	Gene → protein mapping via RefSeq
	•	Protein → UniProt mapping
	•	Protein → PDB structure mapping via RCSB
## Relational MySQL Database
	•	Gene–protein–structure relationships
	•	Primary and foreign keys 
## Web Interface
## REST API Endpoints
	•	Access viruses, genes, proteins, and structures programmatically
  
# Requirements
	•	Python 3.10+
	•	MySQL Server
	•	Python packages:
	•	Flask
	•	Biopython
	•	requests
	•	mysql-connector-python

# Installation
## Clone the Repository 
git clone https://github.com/your-username/Coronavirus-Database.git
## Set UP the MySQL Databse (configurations)
## Run the ETL Pipline (extract.py)
## Run the Web Application (app.py)
## API Endpoints
• URL:/api/genes
• URL:/api/viruses
• URL:/api/proteins

# Other Data
Environment.yml
DataBase Dump
EER (Enhanced Entity–Relationship Diagram)
gitignore 
