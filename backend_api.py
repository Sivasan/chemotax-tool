# backend_api.py - Flask Backend for Multi-Annotation Chemotaxonomy Platform
# -------------------------------------------------------------------------
# Features:
# - Parses multiple annotation formats (GFF3/Prokka, RAST TSV, Bakta TSV,
#   DRAM TSV, HMMER tblout, GenBank GBK, generic TSV).
# - Detects chemotaxonomic markers (quinones, fatty acids, polar lipids,
#   cell-wall polymers, sphingolipids, hopanoids, carotenoids, etc.).
# - Optional phylum filter to limit searches to taxa-relevant markers.
# - Generates JSON results and a Markdown -> PDF/DOCX report via Pandoc.
# - “AI report” is optional; if no API key is configured, it safely falls back.
#
# Requirements (requirements.txt minimal set):
# flask
# flask-cors
# biopython
# pypandoc
# requests
# gunicorn
#
# Run (local):
#   export PORT=5000
#   python backend_api.py
#
# Docker (example):
#   docker build -t chemotax:latest .
#   docker run -p 10000:10000 -e PORT=10000 chemotax:latest
# -------------------------------------------------------------------------

from flask import Flask, request, jsonify, send_file, render_template_string
from flask_cors import CORS
from Bio import SeqIO

import io
import re
import csv
import os
import tempfile
import requests
import json
import pypandoc
from datetime import datetime
from collections import defaultdict
from typing import List, Optional

# -----------------------------------------------------------------------------
# Flask app
# -----------------------------------------------------------------------------
app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": "*"}}, expose_headers=["Content-Disposition"])

# -----------------------------------------------------------------------------
# Annotation Parsing
# -----------------------------------------------------------------------------
class AnnotationParser:
    @staticmethod
    def parse_gff3(content, source_name='gff3'):
        genes = []
        for line in content.split('\n'):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            attributes = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key.strip()] = value.strip()
            def _to_int(s):
                try:
                    return int(s)
                except Exception:
                    return None
            gene_info = {
                'seqid': fields[0],
                'source': fields[1],
                'type': fields[2],
                'start': _to_int(fields[3]),
                'end': _to_int(fields[4]),
                'strand': fields[6],
                'gene': attributes.get('gene', attributes.get('Gene', attributes.get('Name', ''))),
                'product': attributes.get('product', attributes.get('Product', '')),
                'ec_number': attributes.get('ec_number', attributes.get('eC_number', '')),
                'locus_tag': attributes.get('locus_tag', attributes.get('ID', '')),
                'db_xref': attributes.get('Dbxref', attributes.get('db_xref', '')),
                'inference': attributes.get('inference', ''),
                'note': attributes.get('note', ''),
                'format': source_name
            }
            genes.append(gene_info)
        return genes

    @staticmethod
    def parse_rast(content):
        genes = []
        reader = csv.DictReader(io.StringIO(content), delimiter='\t')
        for row in reader:
            genes.append({
                'contig': row.get('contig', row.get('Contig', '')),
                'type': row.get('feature_type', row.get('Feature Type', 'CDS')),
                'start': row.get('start', row.get('Start', '')),
                'end': row.get('end', row.get('End', '')),
                'strand': row.get('strand', row.get('Strand', '')),
                'product': row.get('function', row.get('Function', row.get('annotation', ''))),
                'gene': row.get('gene', row.get('Gene', '')),
                'ec_number': row.get('ec', row.get('EC', '')),
                'format': 'rast'
            })
        return genes

    @staticmethod
    def parse_bakta_tsv(content):
        genes = []
        reader = csv.DictReader(io.StringIO(content), delimiter='\t')
        for row in reader:
            genes.append({
                'sequence_id': row.get('#Sequence Id', ''),
                'type': row.get('Type', 'CDS'),
                'start': row.get('Start', ''),
                'end': row.get('Stop', ''),
                'strand': row.get('Strand', ''),
                'locus_tag': row.get('Locus Tag', ''),
                'gene': row.get('Gene', ''),
                'product': row.get('Product', ''),
                'db_xref': row.get('DbXrefs', ''),
                'ec_number': AnnotationParser._extract_ec_from_dbxref(row.get('DbXrefs', '')),
                'format': 'bakta'
            })
        return genes

    @staticmethod
    def _extract_ec_from_dbxref(dbxref):
        if not dbxref:
            return ''
        m = re.search(r'EC:([\d\.-]+)', dbxref)
        return m.group(1) if m else ''

    @staticmethod
    def parse_dram(content):
        genes = []
        reader = csv.DictReader(io.StringIO(content), delimiter='\t')
        for row in reader:
            pfam_hits = row.get('pfam_hits', '') or row.get('pfam', '')
            pfams = []
            if pfam_hits:
                for tok in re.split(r'[;,|\s]+', pfam_hits):
                    tok = tok.strip()
                    m = re.match(r'(PF\d{5})(?:\.\d+)?', tok)
                    if m:
                        pfams.append(m.group(1))
            genes.append({
                'gene_id': row.get('gene_id', row.get('fasta', '')),
                'scaffold': row.get('scaffold', ''),
                'start': row.get('start_position', ''),
                'end': row.get('end_position', ''),
                'product': row.get('product', row.get('kegg_hit', row.get('pfam_hits', ''))),
                'gene': row.get('gene_name', ''),
                'ec_number': row.get('ec_number', row.get('kegg_ec', '')),
                'kegg_id': row.get('kegg_id', ''),
                'pfam_list': pfams,
                'format': 'dram'
            })
        return genes

    @staticmethod
    def parse_hmmer(content):
        genes = []
        for line in content.split('\n'):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split()

            # Heuristic detection of format
            # domtblout has at least 22 columns. tblout has at least 18.
            # domtblout col 2 is tlen (int). tblout col 2 is query name (string).

            is_domtblout = False
            if len(fields) >= 22:
                try:
                    int(fields[2]) # tlen
                    is_domtblout = True
                except ValueError:
                    pass

            if is_domtblout:
                # domtblout indices
                target_name = fields[0]
                target_acc = fields[1]
                query_name = fields[3]
                try:
                    evalue = float(fields[6])
                except Exception:
                    evalue = 1.0
                try:
                    score = float(fields[7])
                except Exception:
                    score = 0.0
                product = ' '.join(fields[22:])
            else:
                # tblout indices
                if len(fields) < 18: continue
                target_name = fields[0]
                target_acc = fields[1]
                query_name = fields[2]
                try:
                    evalue = float(fields[4])
                except Exception:
                    evalue = 1.0
                try:
                    score = float(fields[5])
                except Exception:
                    score = 0.0
                product = ' '.join(fields[18:])

            pfam_acc = target_acc if target_acc.startswith('PF') else ''
            genes.append({
                'target_name': target_name,
                'target_accession': target_acc,
                'query_name': query_name,
                'evalue': evalue,
                'score': score,
                'product': product,
                'gene': query_name,
                'pfam': pfam_acc,
                'pfam_list': [pfam_acc] if pfam_acc else [],
                'format': 'hmm'
            })
        return genes

    @staticmethod
    def parse_genbank(content):
        genes = []
        try:
            handle = io.StringIO(content)
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type in ["CDS", "gene"]:
                        q = feature.qualifiers
                        start = int(feature.location.start) if feature.location else None
                        end = int(feature.location.end) if feature.location else None
                        strand = feature.location.strand if feature.location else None
                        genes.append({
                            'type': feature.type,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'gene': (q.get('gene', ['']) or [''])[0],
                            'product': (q.get('product', ['']) or [''])[0],
                            'locus_tag': (q.get('locus_tag', ['']) or [''])[0],
                            'ec_number': (q.get('EC_number', ['']) or [''])[0],
                            'format': 'genbank'
                        })
        except Exception as e:
            print(f"Error parsing GenBank: {e}")
        return genes

    @staticmethod
    def parse_generic_tsv(content):
        genes = []
        try:
            reader = csv.DictReader(io.StringIO(content), delimiter='\t')
            headers = reader.fieldnames
            if not headers:
                return genes
            lower_to_orig = {h.lower(): h for h in headers}
            header_map = {
                'gene': ['gene', 'gene_name', 'gene_id', 'gene_symbol'],
                'product': ['product', 'function', 'annotation', 'description'],
                'ec_number': ['ec', 'ec_number', 'ec_num', 'enzyme'],
                'locus_tag': ['locus_tag', 'locus', 'id', 'feature_id']
            }
            for row in reader:
                gene_info = {'format': 'tsv'}
                for std, variants in header_map.items():
                    for v in variants:
                        if v in lower_to_orig:
                            gene_info[std] = row.get(lower_to_orig[v], '')
                            break
                if not gene_info.get('gene') and not gene_info.get('product'):
                    continue
                genes.append(gene_info)
        except Exception as e:
            print(f"Error parsing TSV: {e}")
        return genes

    @staticmethod
    def parse_fasta(content):
        genes = []
        try:
            handle = io.StringIO(content)
            for record in SeqIO.parse(handle, "fasta"):
                # Header format varies:
                # >ID product
                # >ID [gene=...] [protein=...] (NCBI style)
                # >ID ... product=... (Prokka style sometimes)

                desc = record.description
                # Remove ID from description if present at start
                if desc.startswith(record.id):
                    desc = desc[len(record.id):].strip()

                gene = ''
                product = desc
                ec_number = ''

                # Try parsing attributes like [gene=X] or product=Y
                # NCBI style
                m_gene = re.search(r'\[gene=([^\]]+)\]', desc)
                if m_gene:
                    gene = m_gene.group(1)

                m_prod = re.search(r'\[protein=([^\]]+)\]', desc) # NCBI
                if m_prod:
                    product = m_prod.group(1)
                else:
                    # Prokka style: product=...
                    m_prod_pk = re.search(r'product=(.+?)(?:;|$)', desc)
                    if m_prod_pk:
                        product = m_prod_pk.group(1)

                # EC number?
                m_ec = re.search(r'\[EC_number=([^\]]+)\]', desc) # NCBI sometimes
                if not m_ec:
                    m_ec = re.search(r'EC_number=([\d\.-]+)', desc)
                if m_ec:
                    ec_number = m_ec.group(1)

                genes.append({
                    'id': record.id,
                    'gene': gene,
                    'product': product,
                    'ec_number': ec_number,
                    'format': 'fasta'
                })
        except Exception as e:
            print(f"Error parsing FASTA: {e}")
        return genes

# -----------------------------------------------------------------------------
# Marker Detection
# -----------------------------------------------------------------------------
class MarkerDetector:
    def __init__(self):
        self.marker_database = self._load_marker_database()
        self.alias_to_canon, self.canon_to_aliases = self._build_alias_maps()

    def _build_alias_maps(self):
        # Minimal alias table; extend as needed
        alias_table = {
            'menG': {'ubiE'},
            'ubiE': {'menG'},
            'idsA': {'ddsA'},
            'ddsA': {'idsA'},
            # quinone O2-independent helpers sometimes named differently
            'ubit': set(), 'ubiu': set(), 'ubiv': set(), 'ubik': set(),
        }
        alias_to_canon = {}
        canon_to_aliases = defaultdict(set)
        for canon, aliases in alias_table.items():
            canon_to_aliases[canon].update(aliases)
            for a in aliases:
                alias_to_canon[a] = canon
            alias_to_canon[canon] = canon
        return alias_to_canon, canon_to_aliases

    @staticmethod
    def _clean_token(s: str) -> str:
        return re.sub(r'[^a-z0-9]+', '', (s or '').lower())

    def _canon_gene(self, name: str, product: str = '') -> str:
        if not name and not product:
            return ''
        n = self._clean_token(name)
        if n in self.alias_to_canon:
            return self.alias_to_canon[n]
        prod = (product or '').lower()
        for alias in self.alias_to_canon.keys():
            if alias and re.search(r'\b' + re.escape(alias) + r'\b', prod):
                return self.alias_to_canon[alias]
        return n

    # -------------------- EXPANDED DATABASE --------------------
    def _load_marker_database(self):
        # Expanded marker coverage across key chemotaxonomic domains.
        # Each subcategory lists representative genes/EC/keywords/PFAMs.
        return {
            'fattyAcids': {
                'name': 'Fatty Acids & Related Lipids',
                'mycolic': {
                    'genes': ['pks13','fadD32','accD4','acpM','kasA','kasB','inhA','mabA','hadA','hadB','hadC','fbpA','fbpB','fbpC','mmpL3'],
                    'ec': ['1.3.1.9','1.1.1.100','2.3.1.41'],
                    'keywords': ['mycolic acid','FAS-II','mycoloyltransferase'],
                    'taxa': ['Actinobacteria','Actinobacteriota','Corynebacteriales'],
                    'pfam': []
                },
                'anteiso_bfa': {
                    'genes': ['bkdAA','bkdAB','bkdC','ilvE','leuA','leuB','leuC','leuD','fabH'],
                    'ec': ['1.2.4.4','2.3.1.168','1.3.1.60','2.6.1.42','2.3.3.13','4.2.1.33'],
                    'keywords': ['branched-chain fatty acid','anteiso-','BCFA','2-methylbutyryl-CoA'],
                    'taxa': ['Bacillota','Firmicutes'],
                    'pfam': []
                },
                'unsaturated_UFA': {
                    'genes': ['fabA','fabB','desA','desB'],
                    'ec': [],
                    'keywords': ['unsaturated fatty acid','UFA','desaturase'],
                    'taxa': [],
                    'pfam': []
                },
                'cyclopropane_CFA': {
                    'genes': ['cfa'],
                    'ec': ['2.1.1.-'],
                    'keywords': ['cyclopropane fatty acid','CFA synthase'],
                    'taxa': [],
                    'pfam': []
                }
            },

            'quinones': {
                'name': 'Respiratory Quinones',
                'menaquinone_classical': {
                    'genes': ['menF','menD','menH','menC','menE','menB','menI','menA','menG','ubiE','menJ'],
                    'ec': ['2.5.1.74','6.2.1.26','1.2.1.45','4.1.3.36','4.1.1.63','4.2.1.113','2.5.1.30','2.1.1.163'],
                    'keywords': ['menaquinone','vitamin K2','demethylmenaquinone'],
                    'taxa': ['Firmicutes','Bacillota','Actinobacteria','Actinobacteriota','Bacteroidetes','Bacteroidota','Proteobacteria','Pseudomonadota'],
                    'pfam': []
                },
                'menaquinone_futalosine': {
                    'genes': ['mqnA','mqnB','mqnC','mqnD','mqnE','mqnL'],
                    'ec': [],
                    'keywords': ['futalosine pathway','alternative menaquinone'],
                    'taxa': [],
                    'pfam': []
                },
                'ubiquinone': {
                    'genes': ['ubiA','ubiB','ubiC','ubiD','ubiX','ubiE','ubiF','ubiG','ubiH','ubiI','ubiJ','coq7'],
                    'ec': ['2.5.1.39','2.1.1.64','1.1.1.152','3.3.2.11','4.1.1.64','2.1.1.114','1.14.13.36'],
                    'keywords': ['ubiquinone','coenzyme Q','4-hydroxybenzoate'],
                    'taxa': ['Proteobacteria','Pseudomonadota','Acidobacteria','Chlorobi'],
                    'pfam': []
                },
                'ubiquinone_O2_independent': {
                    'genes': ['ubiT','ubiU','ubiV','ubiK'],
                    'ec': [],
                    'keywords': ['oxygen-independent ubiquinone','anaerobic ubiquinone','UbiU','UbiV','UbiT'],
                    'taxa': [],
                    'pfam': []
                }
            },

            'polarLipids': {
                'name': 'Polar & Amino Lipids',
                'phosphatidylethanolamine_PE': {
                    'genes': ['pssA','psd'],
                    'ec': ['2.7.8.8','4.1.1.65'],
                    'keywords': ['phosphatidylethanolamine','PE','phosphatidylserine decarboxylase','CDP-diacylglycerol-serine O-phosphatidyltransferase'],
                    'taxa': ['Proteobacteria','Pseudomonadota'],
                    'pfam': []
                },
                'phosphatidylglycerol_PG': {
                    'genes': ['pgsA','pgpA','pgpB'],
                    'ec': ['2.7.8.5','3.1.3.27'],
                    'keywords': ['phosphatidylglycerol','PG','phosphatidylglycerophosphate synthase'],
                    'taxa': [],
                    'pfam': []
                },
                'cardiolipin_CL': {
                    'genes': ['clsA','clsB','clsC','plsC'],
                    'ec': ['2.7.8.-','2.3.1.15','3.1.1.3'],
                    'keywords': ['cardiolipin','diphosphatidylglycerol','cardiolipin synthase'],
                    'taxa': [],
                    'pfam': []
                },
                'phosphatidylcholine_PC': {
                    'genes': ['pcs','pmtA','pmtB'],
                    'ec': [],
                    'keywords': ['phosphatidylcholine','PC','phospholipid N-methyltransferase','CDP-choline'],
                    'taxa': [],
                    'pfam': []
                },
                'ornithine_lysine_lipids': {
                    'genes': ['olsA','olsB','olsE'],
                    'ec': [],
                    'keywords': ['ornithine lipid','lysine lipid','OL','LL'],
                    'taxa': [],
                    'pfam': []
                },
                'sphingolipids': {
                    'genes': ['sptA','sptB','cerS'],
                    'ec': [],
                    'keywords': ['sphingolipid','ceramide','serine palmitoyltransferase'],
                    'taxa': ['Bacteroidetes','Bacteroidota'],
                    'pfam': []
                }
            },

            'cellWall': {
                'name': 'Cell Wall & Surface Polymers',
                'peptidoglycan_A1_gamma': {
                    'genes': ['murA','murB','murC','murD','murE','murF','murG','mraY','murJ','ddl','mepM','pbp'],
                    'ec': ['2.5.1.7','6.3.2.4','6.3.2.7','6.3.2.13','6.3.2.15','6.3.2.10','2.4.1.129','2.7.8.13','3.4.16.4'],
                    'keywords': ['peptidoglycan','murein','D-alanyl-D-alanine ligase','m-DAP'],
                    'taxa': ['Proteobacteria','Pseudomonadota'],
                    'pfam': []
                },
                'peptidoglycan_A4_alpha': {
                    'genes': ['murA','murB','murC','murD','murE','murF','murG','mraY','murJ','ddl','lysA','pbp'],
                    'ec': ['2.5.1.7','6.3.2.4','6.3.2.7','6.3.2.13','6.3.2.15','6.3.2.8','2.4.1.129','2.7.8.13','4.1.1.20'],
                    'keywords': ['peptidoglycan','murein','DAP-pathway','meso-diaminopimelate','L-Lys'],
                    'taxa': ['Bacillota','Firmicutes'],
                    'pfam': []
                },
                'teichoic_acid': {
                    'genes': ['tagA','tagB','tagD','tagE','tagF','tagO','dltA','dltB','dltC','dltD','ltaS'],
                    'ec': ['2.7.8.12','2.4.1.55','6.1.1.13','2.3.1.196'],
                    'keywords': ['teichoic acid','wall teichoic acid','lipoteichoic acid','D-alanylation'],
                    'taxa': ['Bacillota','Firmicutes'],
                    'pfam': []
                },
                'LPS_LipidA_core': {
                    'genes': ['lpxA','lpxB','lpxC','lpxD','lpxK','kdtA','lpxL','lpxM'],
                    'ec': [],
                    'keywords': ['lipid A','Kdo transferase','LPS core'],
                    'taxa': ['Proteobacteria','Pseudomonadota'],
                    'pfam': []
                },
                'LPS_O_antigen': {
                    'genes': ['waaC','waaF','waaG','waaL','wzx','wzy','wzz'],
                    'ec': [],
                    'keywords': ['O-antigen','LPS O-antigen','O-antigen polymerase'],
                    'taxa': ['Proteobacteria','Pseudomonadota'],
                    'pfam': []
                },
                'arabinogalactan_LAM': {
                    'genes': ['embA','embB','embC','glfT1','glfT2','aftA','aftB','aftC','dprE1','dprE2','ubiA'],
                    'ec': [],
                    'keywords': ['arabinogalactan','lipoarabinomannan','decaprenylphosphoryl arabinose'],
                    'taxa': ['Actinobacteria','Corynebacteriales','Actinobacteriota'],
                    'pfam': []
                }
            },

            'triterpenoids': {
                'name': 'Membrane Triterpenoids',
                'hopanoids': {
                    'genes': ['shc','hpnC','hpnD','hpnE','hpnF','hpnH','hpnI','hpnJ','hpnK'],
                    'ec': [],
                    'keywords': ['hopanoid','squalene-hopene cyclase','bacterial hopanoids'],
                    'taxa': [],
                    'pfam': []
                }
            },

            'pigments': {
                'name': 'Pigments',
                'carotenoids': {
                    'genes': ['crtE','crtB','crtI','crtY','crtZ'],
                    'ec': [],
                    'keywords': ['carotenoid','lycopene','beta-carotene'],
                    'taxa': [],
                    'pfam': []
                }
            }
        }

    def find_markers(self, all_genes, allowed_taxa: Optional[List[str]] = None):
        found_markers = defaultdict(lambda: {'subcategories': {}})
        allowed = None
        if allowed_taxa:
            allowed = {t.strip().lower() for t in allowed_taxa if t}

        gene_index = defaultdict(list)
        product_index = defaultdict(list)
        ec_index = defaultdict(list)
        pfam_index = defaultdict(list)
        product_fulltext = []

        for idx, gene in enumerate(all_genes):
            gname = (gene.get('gene') or '')
            product = gene.get('product', '') or gene.get('function', '')
            canon = self._canon_gene(gname, product)
            if canon:
                gene['gene_norm'] = canon
                gene_index[canon].append((idx, gene))
                for al in self.canon_to_aliases.get(canon, []):
                    gene_index[al].append((idx, gene))
            if product:
                tokens = re.findall(r'[\w\.-]+', product.lower())
                for token in tokens:
                    if len(token) > 3:
                        product_index[token].append((idx, gene))
                product_fulltext.append((idx, product.lower()))
            ec = gene.get('ec_number', '')
            if ec:
                for ec_num in re.findall(r'[\d\.-]+', ec):
                    if ec_num.count('.') >= 2:  # basic EC check
                        ec_index[ec_num].append((idx, gene))
            # HMMER/DRAM PFAM hits
            if gene.get('pfam_list'):
                for pf in gene.get('pfam_list'):
                    pfam_index[pf].append((idx, gene))
            if gene.get('pfam'):
                pfam_index[gene.get('pfam')].append((idx, gene))

        def _ec_matches_query(query_ec: str):
            if '-' not in query_ec and '*' not in query_ec:
                return ec_index.get(query_ec, [])
            prefix = query_ec.replace('*', '').split('-')[0].rstrip('.')
            matched = []
            for full_ec, lst in ec_index.items():
                if full_ec.startswith(prefix):
                    matched.extend(lst)
            return matched

        for category, cat_data in self.marker_database.items():
            found_markers[category]['name'] = cat_data['name']
            for subcat, marker_data in cat_data.items():
                if subcat == 'name':
                    continue
                if allowed is not None:
                    taxa_list = [t.lower() for t in marker_data.get('taxa', [])]
                    if taxa_list and allowed.isdisjoint(taxa_list):
                        continue

                found_set, sources = set(), []

                # Genes
                for marker_gene in marker_data.get('genes', []):
                    mg = self._clean_token(marker_gene)
                    names_to_try = {mg}
                    names_to_try |= self.canon_to_aliases.get(mg, set())
                    if mg in self.alias_to_canon:
                        names_to_try.add(self.alias_to_canon[mg])
                        names_to_try |= self.canon_to_aliases.get(self.alias_to_canon[mg], set())
                    for name_try in names_to_try:
                        matches = gene_index.get(name_try, [])
                        for _, gene in matches:
                            found_set.add(marker_gene)
                            sources.append({
                                'marker': marker_gene,
                                'source': gene.get('format', 'unknown'),
                                'gene': gene.get('gene') or gene.get('gene_norm', ''),
                                'product': gene.get('product', '')
                            })

                # EC
                for ec_num in marker_data.get('ec', []):
                    matches = _ec_matches_query(ec_num)
                    for _, gene in matches:
                        marker_id = f"EC:{ec_num}"
                        found_set.add(marker_id)
                        sources.append({
                            'marker': marker_id,
                            'source': gene.get('format', 'unknown'),
                            'ec': ec_num,
                            'product': gene.get('product', '')
                        })

                # Keywords
                for keyword in marker_data.get('keywords', []):
                    kw_lower = keyword.lower()
                    kw_tokens = set(re.findall(r'[\w\.-]+', kw_lower))
                    potential_matches = set()
                    for token in kw_tokens:
                        if token in product_index:
                            potential_matches.update(idx for idx, _ in product_index[token])
                    if len(kw_tokens) > 1:
                        for idx, full_product_text in product_fulltext:
                            if kw_lower in full_product_text:
                                potential_matches.add(idx)
                    for idx in potential_matches:
                        gene = all_genes[idx]
                        product_tokens = set(re.findall(r'[\w\.-]+', (gene.get('product', '') or '').lower()))
                        if kw_tokens.issubset(product_tokens) or kw_lower in (gene.get('product', '') or '').lower():
                            marker_id = f"kw:{keyword}"
                            found_set.add(marker_id)
                            sources.append({
                                'marker': marker_id,
                                'source': gene.get('format', 'unknown'),
                                'keyword': keyword,
                                'product': gene.get('product', '')
                            })

                # PFAMs
                for pf in marker_data.get('pfam', []):
                    matches = pfam_index.get(pf, [])
                    for _, gene in matches:
                        marker_id = f"pfam:{pf}"
                        found_set.add(marker_id)
                        sources.append({
                            'marker': marker_id,
                            'source': gene.get('format', 'unknown'),
                            'product': gene.get('product', ''),
                            'pfam': pf
                        })

                if found_set:
                    total_genes = set(marker_data.get('genes', []))
                    total_ec = set(marker_data.get('ec', []))
                    total_kw = set(marker_data.get('keywords', []))
                    total_pfam = set(marker_data.get('pfam', []))
                    total = len(total_genes | total_ec | total_kw | total_pfam)

                    found_markers[category]['subcategories'][subcat] = {
                        'found': list(found_set),
                        'total': total,
                        'percentage': round((len(found_set) / total * 100), 1) if total else 0.0,
                        'sources': sources
                    }

        found_markers = {k: v for k, v in found_markers.items() if v['subcategories']}
        return dict(found_markers)

# -----------------------------------------------------------------------------
# Predictions
# -----------------------------------------------------------------------------
def get_pathway_completeness(found_markers, category, subcategory):
    try:
        return float(found_markers[category]['subcategories'][subcategory]['percentage'])
    except KeyError:
        return 0.0

def _subcat_data(found_markers, cat, sub):
    return found_markers.get(cat, {}).get('subcategories', {}).get(sub, {'found': [], 'sources': []})

def generate_predictions(found_markers):
    """
    Rules support:
      - 'subcategory' (single)
      - 'subcategories' (list of (category, subcategory) tuples)
      - 'alt_groups': list of lists of (category, subcategory) tuples (use the best in each group)
    """
    rules = [
        # Corynebacteriales-type wall signature
        {
            'id': 'corynebacteriales_wall',
            'name': 'Corynebacteriales Cell Envelope (Mycolic + Arabinogalactan/LAM)',
            'ijsem': 'Cell wall with arabinogalactan–mycolic acids (Corynebacteriales-type).',
            'min_completeness': 25,
            'subcategories': [('fattyAcids', 'mycolic'), ('cellWall', 'arabinogalactan_LAM')]
        },

        # Gram-Positive (Firmicutes-type)
        {
            'id': 'gram_positive_profile',
            'name': 'Gram-Positive Profile (Firmicutes-type)',
            'ijsem': 'Typical Gram-positive cell wall polymer (teichoic acids) and MK quinones.',
            'min_completeness': 30,
            'subcategories': [('cellWall', 'peptidoglycan_A4_alpha'), ('cellWall', 'teichoic_acid')],
            'alt_groups': [[('quinones', 'menaquinone_classical'), ('quinones', 'menaquinone_futalosine')]]
        },

        # Gram-Negative (Proteobacteria-type)
        {
            'id': 'gram_negative_profile',
            'name': 'Gram-Negative Profile (Proteobacteria-type)',
            'ijsem': 'Outer membrane with LPS (lipid A) and ubiquinone quinones.',
            'min_completeness': 30,
            'subcategories': [('cellWall', 'peptidoglycan_A1_gamma'), ('cellWall', 'LPS_LipidA_core')],
            'alt_groups': [[('quinones', 'ubiquinone'), ('quinones', 'ubiquinone_O2_independent')]]
        },

        # Sphingolipid-rich Bacteroidota tendency
        {
            'id': 'bacteroidota_sphingo',
            'name': 'Sphingolipid-Positive Profile (Bacteroidota tendency)',
            'ijsem': 'Presence of bacterial sphingolipids (e.g., serine palmitoyltransferase).',
            'min_completeness': 20,
            'subcategories': [('polarLipids', 'sphingolipids')]
        },

        # Hopanoids
        {
            'id': 'hopanoid_presence',
            'name': 'Hopanoid Biosynthesis',
            'ijsem': 'Hopanoid biosynthetic capacity detected.',
            'min_completeness': 20,
            'subcategories': [('triterpenoids', 'hopanoids')]
        },

        # Carotenoids
        {
            'id': 'carotenoid_presence',
            'name': 'Carotenoid Pigments',
            'ijsem': 'Carotenoid pigment biosynthetic genes present.',
            'min_completeness': 20,
            'subcategories': [('pigments', 'carotenoids')]
        },

        # Classic mycolic-only rule for Actinobacteria signal
        {
            'id': 'actinobacteria_mycolic',
            'name': 'Mycolic Acid Profile (Actinobacteria-type)',
            'ijsem': 'Cell wall contains mycolic acids.',
            'min_completeness': 30,
            'category': 'fattyAcids',
            'subcategory': 'mycolic'
        }
    ]

    predictions = []

    for rule in rules:
        matched_markers = []
        sources = []
        total_completeness = 0.0
        group_count = 0

        # single subcategory rule
        if 'subcategory' in rule and 'category' in rule:
            completeness = get_pathway_completeness(found_markers, rule['category'], rule['subcategory'])
            if completeness > 0:
                sd = _subcat_data(found_markers, rule['category'], rule['subcategory'])
                matched_markers.extend(sd.get('found', []))
                sources.extend(sd.get('sources', []))
                total_completeness += completeness
                group_count += 1

        # multi required subcategories
        for cat, sub in rule.get('subcategories', []):
            completeness = get_pathway_completeness(found_markers, cat, sub)
            # count the group whether present or not, to penalize missing modules less aggressively we still include it,
            # but you can switch to: if completeness > 0: group_count += 1
            group_count += 1
            if completeness > 0:
                sd = _subcat_data(found_markers, cat, sub)
                matched_markers.extend(sd.get('found', []))
                sources.extend(sd.get('sources', []))
                total_completeness += completeness

        # alternative groups: choose best per group
        for group in rule.get('alt_groups', []):
            group_best = 0.0
            group_best_sd = None
            for cat, sub in group:
                c = get_pathway_completeness(found_markers, cat, sub)
                if c > group_best:
                    group_best = c
                    group_best_sd = _subcat_data(found_markers, cat, sub) if c > 0 else None
            # always count a group (present or not) for averaging
            group_count += 1
            if group_best_sd:
                matched_markers.extend(group_best_sd.get('found', []))
                sources.extend(group_best_sd.get('sources', []))
                total_completeness += group_best

        if group_count > 0:
            avg_completeness = total_completeness / group_count
            if avg_completeness >= rule['min_completeness']:
                predictions.append({
                    **{k: v for k, v in rule.items() if k not in ('subcategories', 'alt_groups', 'category', 'subcategory')},
                    'matchedMarkers': sorted(set(matched_markers)),
                    'confidence': round(avg_completeness, 1),
                    'sources': sources
                })

    return sorted(predictions, key=lambda x: x['confidence'], reverse=True)

def identify_novel_features(found_markers):
    novel, complete = [], []
    for category, cat_data in found_markers.items():
        for subcat, subcat_data in cat_data['subcategories'].items():
            if float(subcat_data['percentage']) >= 70:
                complete.append({
                    'category': cat_data['name'],
                    'subcategory': subcat,
                    'percentage': subcat_data['percentage'],
                    'markers': len(subcat_data['found'])
                })
    if complete:
        novel.append({
            'type': 'Highly Complete Pathways',
            'description': 'Pathways with ≥70% marker detection',
            'features': complete
        })
    return novel

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _normalize_phylum_param(raw_val: Optional[str]) -> Optional[List[str]]:
    if not raw_val:
        return None
    vals = [v.strip() for v in raw_val.split(',') if v.strip()]
    synonyms = {
        'actinobacteriota': 'Actinobacteria',
        'actino': 'Actinobacteria',
        'proteobacteria': 'Proteobacteria',
        'pseudomonadota': 'Proteobacteria',
        'firmicutes': 'Firmicutes',
        'bacillota': 'Firmicutes',
        'bacteroidetes': 'Bacteroidetes',
        'bacteroidota': 'Bacteroidetes',
        'cyanobacteria': 'Cyanobacteria',
        'chloroflexi': 'Chloroflexi',
        'chloroflexota': 'Chloroflexi',
        'corynebacteriales': 'Corynebacteriales'
    }
    mapped = []
    for v in vals:
        v_lower = v.lower()
        mapped.append(synonyms.get(v_lower, v.strip()))
    seen, out = set(), []
    for x in mapped:
        if x not in seen:
            out.append(x)
            seen.add(x)
    return out or None

def _analyze_uploaded(files, formats, phylum_param: Optional[str]):
    parser = AnnotationParser()
    detector = MarkerDetector()

    all_genes = []
    file_info = []

    # Store genes per file for individual analysis
    # Structure: [{'name': 'f1', 'genes': [...]}, ...]
    files_genes_map = []

    for file, format_type in zip(files, formats):
        try:
            content = file.read().decode('utf-8', errors='replace')
        except Exception as e:
            return {'error': f"Error reading file {file.filename}: {e}"}

        filename = file.filename
        genes = []
        ft = (format_type or '').lower()
        if ft in ['gff3', 'prokka']:
            genes = parser.parse_gff3(content, ft)
        elif ft == 'rast':
            genes = parser.parse_rast(content)
        elif ft == 'bakta':
            genes = parser.parse_bakta_tsv(content)
        elif ft == 'dram':
            genes = parser.parse_dram(content)
        elif ft == 'hmm':
            genes = parser.parse_hmmer(content)
        elif ft == 'gbk':
            genes = parser.parse_genbank(content)
        elif ft == 'tsv':
            genes = parser.parse_generic_tsv(content)
        elif ft in ['fasta', 'faa', 'fa']:
            genes = parser.parse_fasta(content)
        else:
            return {'error': f'Unsupported format: {format_type} for {filename}'}

        all_genes.extend(genes)
        file_info.append({'name': filename, 'format': ft, 'geneCount': len(genes)})
        files_genes_map.append({'name': filename, 'genes': genes})

    if not all_genes:
        return {'error': "No features could be parsed from the uploaded files. Please check file formats."}

    allowed_taxa = _normalize_phylum_param(phylum_param)

    # 1. Combined analysis
    found_markers = detector.find_markers(all_genes, allowed_taxa=allowed_taxa)
    predictions = generate_predictions(found_markers)
    novel_features = identify_novel_features(found_markers)

    # 2. Individual analysis stats
    individual_results = []
    # Collect all found marker IDs (e.g., 'menF', 'EC:1.2.3.4') from combined result to verify coverage
    combined_marker_ids = set()
    for cat in found_markers.values():
        for sub in cat['subcategories'].values():
            combined_marker_ids.update(sub['found'])

    for entry in files_genes_map:
        f_name = entry['name']
        f_genes = entry['genes']
        # Run detector
        f_markers = detector.find_markers(f_genes, allowed_taxa=allowed_taxa)

        # Count markers found
        f_found_ids = set()
        for cat in f_markers.values():
            for sub in cat['subcategories'].values():
                f_found_ids.update(sub['found'])

        # Comparison with combined
        # Let's just store the count and list of found markers.

        individual_results.append({
            'name': f_name,
            'marker_count': len(f_found_ids),
            'markers': list(f_found_ids)
        })

    # 3. Calculate "Cumulative Gain" or "Contribution"
    # Which markers were found by combining that wouldn't be found if we just took the best single file?
    max_single_file_count = max([r['marker_count'] for r in individual_results]) if individual_results else 0
    combined_count = len(combined_marker_ids)
    gain = combined_count - max_single_file_count

    stats = {
        'combined_count': combined_count,
        'max_single_count': max_single_file_count,
        'gain': gain,
        'individual': individual_results
    }

    result = {
        'files': file_info,
        'totalGenes': len(all_genes),
        'foundMarkers': found_markers,
        'predictions': predictions,
        'novelFeatures': novel_features,
        'timestamp': datetime.now().isoformat(),
        'phylum': allowed_taxa,
        'stats': stats
    }
    return result

# -----------------------------------------------------------------------------
# Reports
# -----------------------------------------------------------------------------
def _build_template_report(result: dict) -> str:
    ts = result.get('timestamp')
    ts_str = ''
    try:
        ts_str = datetime.fromisoformat(ts).strftime('%Y-%m-%d %H:%M')
    except Exception:
        ts_str = datetime.now().strftime('%Y-%m-%d %H:%M')

    lines = []
    lines.append(f"# Genome-resolved Chemotaxonomy Report")
    lines.append(f"_Generated: {ts_str}_")
    if result.get('phylum'):
        lines.append(f"**Scope (phylum/order filter):** {', '.join(result['phylum'])}")
    lines.append("\n---\n")

    lines.append("## 1. Introduction")
    lines.append("This report summarizes the chemotaxonomic features predicted from the provided genome annotations. The analysis searches for key genes and pathways involved in the biosynthesis of taxonomically significant compounds (respiratory quinones, fatty acids, polar lipids, cell wall polymers, sphingolipids, hopanoids, pigments).")

    lines.append("\n## 2. Methods")
    lines.append("Parsed genomic features were searched against a curated database of marker genes, EC numbers, keywords, and PFAMs. The following files were processed:")
    for f in result['files']:
        lines.append(f"- **{f['name']}** ({f['format']}): {f['geneCount']} features")
    lines.append(f"\nA total of **{result['totalGenes']} features** were parsed and analyzed.")

    # Source Comparison Section
    stats = result.get('stats')
    if stats:
        lines.append("\n### 2.1. Multi-Source Contribution Analysis")
        lines.append("By combining annotations from multiple tools, we achieved a more complete chemotaxonomic profile.")
        lines.append(f"- **Combined Unique Markers:** {stats['combined_count']}")
        lines.append(f"- **Best Single Source:** {stats['max_single_count']} markers")
        lines.append(f"- **Cumulative Gain:** +{stats['gain']} markers ({(stats['gain']/stats['max_single_count']*100 if stats['max_single_count'] else 0):.1f}% increase)")

        lines.append("\n| Source File | Marker Count | Contribution |")
        lines.append("|---|---|---|")
        for res in stats['individual']:
            # contribution = unique to this file? Or just raw count? Raw count is simpler.
            lines.append(f"| {res['name']} | {res['marker_count']} | {(res['marker_count']/stats['combined_count']*100 if stats['combined_count'] else 0):.1f}% of total |")

    lines.append("\n## 3. Results & Discussion")
    lines.append("\n### 3.1. Chemotaxonomic Predictions")
    if result['predictions']:
        for p in result['predictions']:
            lines.append(f"- **{p['name']}**")
            lines.append(f"  - **Confidence:** {p['confidence']}%")
            lines.append(f"  - **IJSEM Note:** {p['ijsem']}")
            if p.get('matchedMarkers'):
                preview = ', '.join(list(p['matchedMarkers'])[:8])
                if len(p['matchedMarkers']) > 8:
                    preview += ", …"
                lines.append(f"  - **Key Matched Markers:** {preview}")
    else:
        lines.append("- No high-confidence chemotaxonomic profiles were matched at current thresholds.")

    lines.append("\n### 3.2. Detailed Marker Findings (by category)")
    if result['foundMarkers']:
        for cat, cat_data in result['foundMarkers'].items():
            lines.append(f"\n#### {cat_data['name']}")
            for sub, sub_data in cat_data['subcategories'].items():
                lines.append(f"- **{sub}**: {sub_data['percentage']}% complete ({len(sub_data['found'])}/{sub_data['total']})")
    else:
        lines.append("- No chemotaxonomic markers from the database were detected.")

    lines.append("\n### 3.3. Notable Features")
    if result['novelFeatures']:
        for block in result['novelFeatures']:
            lines.append(f"- **{block['type']}**: {block['description']}")
            for feat in block['features']:
                lines.append(f"  - {feat['category']} / {feat['subcategory']}: {feat['percentage']}% complete")
    else:
        lines.append("- No pathways were flagged as highly complete at current thresholds.")

    lines.append("\n## 4. Conclusion")
    lines.append("This in-silico analysis provides a preliminary chemotaxonomic profile. The predictions and marker findings should be validated with laboratory-based chemical analyses (e.g., HPLC, GC-MS, LC-MS/MS).")

    return "\n".join(lines)

def _generate_ai_report(result: dict) -> str:
    """
    Optional AI narrative. If GEMINI_API_KEY is not set or request fails,
    we return the template report. Designed to never crash.
    """
    fallback_report = _build_template_report(result)
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()
    if not api_key:
        return fallback_report

    try:
        api_url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-preview-09-2025:generateContent?key={api_key}"
        system_prompt = (
            "You are a helpful assistant specialized in microbial taxonomy and chemotaxonomy. "
            "Write a professional, concise narrative report from the provided JSON. "
            "Use Markdown headings: Introduction, Methods, Results and Discussion, Conclusion. "
            "Interpret key predictions (e.g., Gram-negative with LPS & ubiquinone), "
            "discuss pathway completeness and chemotaxonomic meaning."
        )
        prompt_data = {
            'inputs': result.get('files', []),
            'phylum_filter': result.get('phylum', 'None'),
            'predictions': result.get('predictions', []),
            'marker_findings': result.get('foundMarkers', {}),
            'notable_features': result.get('novelFeatures', []),
            'stats': result.get('stats', {})
        }
        user_prompt = (
            f"Generate a full chemotaxonomy report for this JSON. "
            f"Explicitly discuss the value of the multi-annotation approach using the 'stats' block, "
            f"highlighting how combining {len(result.get('files', []))} sources improved the result (Cumulative Gain). "
            f"JSON:\n\n{json.dumps(prompt_data, indent=2)}"
        )
        payload = {
            "contents": [{"parts": [{"text": user_prompt}]}],
            "systemInstruction": {"parts": [{"text": system_prompt}]}
        }
        headers = {"Content-Type": "application/json"}
        resp = requests.post(api_url, json=payload, headers=headers, timeout=60)
        if resp.status_code != 200:
            print(f"AI report generation failed: {resp.status_code} {resp.text}")
            return fallback_report
        data = resp.json()
        text = data.get('candidates', [{}])[0].get('content', {}).get('parts', [{}])[0].get('text', '')
        return text or fallback_report
    except Exception as e:
        print(f"AI report error: {e}")
        return fallback_report

# -----------------------------------------------------------------------------
# API Endpoints
# -----------------------------------------------------------------------------
@app.route('/')
def index():
    try:
        with open('index.html', 'r', encoding='utf-8') as f:
            return render_template_string(f.read())
    except FileNotFoundError:
        return "index.html not found", 404

@app.route('/api/analyze-multi', methods=['POST'])
def analyze_multiple_annotations():
    if 'files' not in request.files:
        return jsonify({'error': 'No files uploaded'}), 400
    files = request.files.getlist('files')
    formats = request.form.getlist('formats')
    phylum_param = request.form.get('phylum', '')
    if len(files) != len(formats):
        return jsonify({'error': 'Number of files and formats must match'}), 400

    result = _analyze_uploaded(files, formats, phylum_param)
    if 'error' in result:
        return jsonify(result), 400
    return jsonify(result)

@app.route('/api/analyze-and-report', methods=['POST'])
def analyze_and_report():
    if 'files' not in request.files:
        return jsonify({'error': 'No files uploaded'}), 400

    files = request.files.getlist('files')
    formats = request.form.getlist('formats')
    phylum_param = request.form.get('phylum', '')
    # Accept format from query OR form (frontend uses query)
    req_format = (request.args.get('format') or request.form.get('format') or 'pdf').lower()
    if req_format not in ['pdf', 'docx']:
        return jsonify({'error': 'Invalid format requested. Must be pdf or docx.'}), 400

    use_ai = (request.form.get('use_ai', 'false').lower() == 'true')

    if len(files) != len(formats):
        return jsonify({'error': 'Number of files and formats must match'}), 400

    # Run the analysis
    result = _analyze_uploaded(files, formats, phylum_param)
    if 'error' in result:
        return jsonify(result), 400

    # Build Markdown (AI or template)
    print(f"Generating report. AI enabled: {use_ai}")
    markdown_string = _generate_ai_report(result) if use_ai else _build_template_report(result)
    print("Markdown generation complete.")

    md_file = None
    out_file = None
    try:
        # Temp markdown
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as md_file:
            md_file.write(markdown_string)
            md_file_path = md_file.name

        # Temp output
        with tempfile.NamedTemporaryFile(suffix=f'.{req_format}', delete=False) as out_file:
            out_file_path = out_file.name

        # Pandoc conversion
        print(f"Converting Markdown to {req_format} using Pandoc...")
        pandoc_args = []
        if req_format == 'pdf':
            # pdflatex is safest with minimal font deps
            pandoc_args = ['--pdf-engine=pdflatex', '-V', 'geometry:margin=1in']

        pypandoc.convert_file(
            md_file_path,
            req_format,
            outputfile=out_file_path,
            extra_args=pandoc_args
        )
        print("Pandoc conversion successful.")

        # Read into memory
        with open(out_file_path, 'rb') as f:
            file_data = io.BytesIO(f.read())

        mimetype_map = {
            'pdf': 'application/pdf',
            'docx': 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
        }
        fname = f"chemotax_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.{req_format}"

        return send_file(
            file_data,
            mimetype=mimetype_map[req_format],
            as_attachment=True,
            download_name=fname
        )

    except Exception as e:
        print(f"Error during file conversion or sending: {e}")
        return jsonify({'error': f'Failed to generate report file: {e}'}), 500
    finally:
        # Cleanup
        try:
            if md_file and os.path.exists(md_file.name):
                os.remove(md_file.name)
        except Exception:
            pass
        try:
            if out_file and os.path.exists(out_file.name):
                os.remove(out_file.name)
        except Exception:
            pass
        print("Temporary files cleaned up.")

@app.route('/health', methods=['GET'])
def health_check():
    return jsonify({'status': 'healthy', 'timestamp': datetime.now().isoformat()})

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    # For local development only (Docker/hosted will use CMD entrypoint)
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 5000)))
