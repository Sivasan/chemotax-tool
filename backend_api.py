# backend_api.py - Flask Backend for Multi-Annotation Chemotaxonomy Platform
from flask import Flask, request, jsonify, send_file, render_template_string
from flask_cors import CORS
from Bio import SeqIO
import io, re, csv, os, tempfile, requests, json
import pypandoc
from datetime import datetime
from collections import defaultdict
from typing import List, Optional

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": "*"}}, expose_headers=["Content-Disposition"])

# --- Annotation Parsing (No Changes) ---
class AnnotationParser:
    @staticmethod
    def parse_gff3(content, source_name='gff3'):
        genes = []
        for line in content.split('\n'):
            if line.startswith('#') or not line.strip(): continue
            fields = line.split('\t')
            if len(fields) < 9: continue
            attributes = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key.strip()] = value.strip()
            gene_info = {
                'seqid': fields[0], 'source': fields[1], 'type': fields[2],
                'start': int(fields[3]) if fields[3].isdigit() else None,
                'end': int(fields[4]) if fields[4].isdigit() else None,
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
                'ec_number': row.get('ec', row.get('EC', '')), 'format': 'rast'
            })
        return genes

    @staticmethod
    def parse_bakta_tsv(content):
        genes = []
        reader = csv.DictReader(io.StringIO(content), delimiter='\t')
        for row in reader:
            genes.append({
                'sequence_id': row.get('#Sequence Id', ''), 'type': row.get('Type', 'CDS'),
                'start': row.get('Start', ''), 'end': row.get('Stop', ''), 'strand': row.get('Strand', ''),
                'locus_tag': row.get('Locus Tag', ''), 'gene': row.get('Gene', ''), 'product': row.get('Product', ''),
                'db_xref': row.get('DbXrefs', ''),
                'ec_number': AnnotationParser._extract_ec_from_dbxref(row.get('DbXrefs', '')),
                'format': 'bakta'
            })
        return genes

    @staticmethod
    def _extract_ec_from_dbxref(dbxref):
        if not dbxref: return ''
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
                    if m: pfams.append(m.group(1))
            genes.append({
                'gene_id': row.get('gene_id', row.get('fasta', '')), 'scaffold': row.get('scaffold', ''),
                'start': row.get('start_position', ''), 'end': row.get('end_position', ''),
                'product': row.get('product', row.get('kegg_hit', row.get('pfam_hits', ''))),
                'gene': row.get('gene_name', ''),
                'ec_number': row.get('ec_number', row.get('kegg_ec', '')), 'kegg_id': row.get('kegg_id', ''),
                'pfam_list': pfams, 'format': 'dram'
            })
        return genes

    @staticmethod
    def parse_hmmer(content):
        genes = []
        for line in content.split('\n'):
            if line.startswith('#') or not line.strip(): continue
            fields = line.strip().split(None, 22) # Split max 22 times
            if len(fields) < 23: continue
            target_acc = fields[1]
            pfam_acc = target_acc if target_acc.startswith('PF') else ''
            genes.append({
                'target_name': fields[0], 'target_accession': target_acc, 'query_name': fields[3],
                'evalue': float(fields[6]), 'score': float(fields[7]),
                'product': fields[22] if len(fields) > 22 else '', # Product is the rest of the line
                'gene': fields[3], # Use query_name as gene
                'pfam': pfam_acc, 'format': 'hmm'
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
                        genes.append({
                            'type': feature.type, 'start': int(feature.location.start), 'end': int(feature.location.end),
                            'strand': feature.location.strand,
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
            if not headers: return genes
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
                    continue # Skip rows with no identifiable gene or product
                genes.append(gene_info)
        except Exception as e:
            print(f"Error parsing TSV: {e}")
        return genes

# --- Marker Detection (Expanded Database) ---
class MarkerDetector:
    def __init__(self):
        self.marker_database = self._load_marker_database()
        self.alias_to_canon, self.canon_to_aliases = self._build_alias_maps()

    def _build_alias_maps(self):
        alias_table = { 'menG': {'ubiE'}, 'ubiE': {'menG'}, 'idsA': {'ddsA'}, 'ddsA': {'idsA'} }
        alias_to_canon = {}
        canon_to_aliases = defaultdict(set)
        for canon, aliases in alias_table.items():
            canon_to_aliases[canon].update(aliases)
            for a in aliases: alias_to_canon[a] = canon
            alias_to_canon[canon] = canon
        return alias_to_canon, canon_to_aliases

    @staticmethod
    def _clean_token(s: str) -> str:
        return re.sub(r'[^a-z0-9]+', '', (s or '').lower())

    def _canon_gene(self, name: str, product: str = '') -> str:
        if not name and not product: return ''
        n = self._clean_token(name)
        if n in self.alias_to_canon: return self.alias_to_canon[n]
        prod = (product or '').lower()
        for alias in self.alias_to_canon.keys():
            if alias and re.search(r'\b' + re.escape(alias) + r'\b', prod):
                return self.alias_to_canon[alias]
        return n

    def _load_marker_database(self):
        return {
            'fattyAcids': {
                'name': 'Fatty Acids & Mycolic Acids',
                'mycolic': {
                    'genes': ['pks13','fadD32','accD4','acpM','kasA','kasB','inhA','mabA','hadA','hadB','hadC', 'fbpA','fbpB','fbpC','mmpL3'],
                    'ec': ['1.3.1.9','1.1.1.100','2.3.1.41'],
                    'keywords': ['mycolic acid','FAS-II','mycoloyltransferase'],
                    'taxa': ['Actinobacteria','Actinobacteriota','Corynebacteriales']
                },
                'anteiso_bfa': {
                    'genes': ['bkdAA', 'bkdAB', 'bkdC', 'ilvE', 'leuA', 'leuB', 'leuC', 'leuD', 'fabH'],
                    'ec': ['1.2.4.4', '2.3.1.168', '1.3.1.60', '2.6.1.42', '2.3.3.13', '4.2.1.33'],
                    'keywords': ['branched-chain fatty acid', 'anteiso-', 'BCFA', '2-methylbutyryl-CoA'],
                    'taxa': ['Bacillota', 'Firmicutes']
                }
            },
            'quinones': {
                'name': 'Respiratory Quinones',
                'menaquinone_classical': {
                    'genes': ['menF','menD','menH','menC','menE','menB','menI','menA','menG','ubiE','menJ'],
                    'ec': ['2.5.1.74','6.2.1.26', '1.2.1.45', '4.1.3.36', '4.1.1.63', '4.2.1.113', '2.5.1.30', '2.1.1.163'],
                    'keywords': ['menaquinone','vitamin K2','demethylmenaquinone'],
                    'taxa': ['Firmicutes','Bacillota','Actinobacteria','Actinobacteriota','Bacteroidetes','Bacteroidota','Proteobacteria','Pseudomonadota']
                },
                'ubiquinone': {
                    'genes': ['ubiA','ubiB','ubiC','ubiD','ubiX','ubiE','ubiF','ubiG','ubiH','ubiI','ubiJ','coq7'],
                    'ec': ['2.5.1.39', '2.1.1.64', '1.1.1.152', '3.3.2.11', '4.1.1.64', '2.1.1.114', '1.14.13.36'],
                    'keywords': ['ubiquinone','coenzyme Q','4-hydroxybenzoate', 'O2-independent ubiquinone'],
                    'taxa': ['Proteobacteria','Pseudomonadota','Acidobacteria','Chlorobi']
                }
            },
            'polarLipids': {
                'name': 'Polar Lipids',
                'phosphatidylethanolamine_PE': {
                    'genes': ['pssA', 'psd'],
                    'ec': ['2.7.8.8', '4.1.1.65'],
                    'keywords': ['phosphatidylethanolamine', 'PE', 'phosphatidylserine decarboxylase', 'CDP-diacylglycerol-serine O-phosphatidyltransferase'],
                    'taxa': ['Proteobacteria', 'Pseudomonadota']
                },
                'phosphatidylglycerol_PG': {
                    'genes': ['pgsA', 'pgpA', 'pgpB'],
                    'ec': ['2.7.8.5', '3.1.3.27'],
                    'keywords': ['phosphatidylglycerol', 'PG', 'phosphatidylglycerophosphate synthase'],
                    'taxa': []
                },
                'cardiolipin_CL': {
                    'genes': ['clsA', 'clsB', 'clsC', 'plsC'],
                    'ec': ['2.7.8.-', '2.3.1.15', '3.1.1.3'],
                    'keywords': ['cardiolipin', 'diphosphatidylglycerol', 'cardiolipin synthase'],
                    'taxa': []
                }
            },
            'cellWall': {
                'name': 'Cell Wall Components',
                'peptidoglycan_A1_gamma': {
                    'genes': ['murA', 'murB', 'murC', 'murD', 'murE', 'murF', 'murG', 'mraY', 'murJ', 'ddl', 'mepM', 'pbp'],
                    'ec': ['2.5.1.7', '6.3.2.4', '6.3.2.7', '6.3.2.13', '6.3.2.15', '6.3.2.10', '2.4.1.129', '2.7.8.13', '3.4.16.4'],
                    'keywords': ['peptidoglycan', 'murein', 'D-alanyl-D-alanine ligase', 'L-Lysine'],
                    'taxa': ['Proteobacteria', 'Pseudomonadota']
                },
                'peptidoglycan_A4_alpha': {
                    'genes': ['murA', 'murB', 'murC', 'murD', 'murE', 'murF', 'murG', 'mraY', 'murJ', 'ddl', 'lysA', 'pbp'],
                    'ec': ['2.5.1.7', '6.3.2.4', '6.3.2.7', '6.3.2.13', '6.3.2.15', '6.3.2.8', '2.4.1.129', '2.7.8.13', '4.1.1.20'],
                    'keywords': ['peptidoglycan', 'murein', 'DAP-pathway', 'meso-diaminopimelate', 'diaminopimelic acid'],
                    'taxa': ['Bacillota', 'Firmicutes']
                },
                'teichoic_acid': {
                    'genes': ['tagA', 'tagB', 'tagD', 'tagE', 'tagF', 'tagO', 'dltA', 'dltB', 'dltC', 'dltD'],
                    'ec': ['2.7.8.12', '2.4.1.55', '6.1.1.13', '2.3.1.196'],
                    'keywords': ['teichoic acid', 'wall teichoic acid', 'lipoteichoic acid', 'D-alanylation'],
                    'taxa': ['Bacillota', 'Firmicutes']
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
                    if ec_num.count('.') >= 2: # Basic EC number validation
                        ec_index[ec_num].append((idx, gene))
            if gene.get('pfam_list'):
                for pf in gene.get('pfam_list'):
                    pfam_index[pf].append((idx, gene))

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
                if subcat == 'name': continue
                if allowed is not None:
                    taxa_list = [t.lower() for t in marker_data.get('taxa', [])]
                    if taxa_list and allowed.isdisjoint(taxa_list):
                        continue
                found_set, sources = set(), []

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
                                'marker': marker_gene, 'source': gene.get('format', 'unknown'),
                                'gene': gene.get('gene') or gene.get('gene_norm', ''), 'product': gene.get('product', '')
                            })

                for ec_num in marker_data.get('ec', []):
                    matches = _ec_matches_query(ec_num)
                    for _, gene in matches:
                        marker_id = f"EC:{ec_num}"
                        found_set.add(marker_id)
                        sources.append({
                            'marker': marker_id, 'source': gene.get('format', 'unknown'),
                            'ec': ec_num, 'product': gene.get('product', '')
                        })
                
                for keyword in marker_data.get('keywords', []):
                    kw_lower = keyword.lower()
                    kw_tokens = set(re.findall(r'[\w\.-]+', kw_lower))
                    
                    # Search token index
                    potential_matches = set()
                    for token in kw_tokens:
                        if token in product_index:
                            potential_matches.update(idx for idx, _ in product_index[token])

                    # Search full-text for multi-word keywords
                    if len(kw_tokens) > 1:
                        for idx, full_product_text in product_fulltext:
                            if kw_lower in full_product_text:
                                potential_matches.add(idx)
                    
                    for idx in potential_matches:
                        gene = all_genes[idx]
                        # Post-filter to ensure all tokens are present for multi-word keywords
                        product_tokens = set(re.findall(r'[\w\.-]+', (gene.get('product', '') or '').lower()))
                        if kw_tokens.issubset(product_tokens) or kw_lower in (gene.get('product', '') or '').lower():
                            marker_id = f"kw:{keyword}"
                            found_set.add(marker_id)
                            sources.append({
                                'marker': marker_id, 'source': gene.get('format', 'unknown'),
                                'keyword': keyword, 'product': gene.get('product', '')
                            })
                            
                for pf in marker_data.get('pfam', []):
                    matches = pfam_index.get(pf, [])
                    for _, gene in matches:
                        marker_id = f"pfam:{pf}"
                        found_set.add(marker_id)
                        sources.append({
                            'marker': marker_id, 'source': gene.get('format', 'unknown'),
                            'product': gene.get('product', ''), 'pfam': pf
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

        found_markers = {k:v for k,v in found_markers.items() if v['subcategories']}
        return dict(found_markers)

# --- Prediction Logic (Expanded) ---
def get_pathway_completeness(found_markers, category, subcategory):
    try:
        return found_markers[category]['subcategories'][subcategory]['percentage']
    except KeyError:
        return 0.0

def generate_predictions(found_markers):
    rules = [
        {'id':'actinobacteria_mycolic','name':'Mycolic Acid Profile (Actinobacteria-type)','category':'fattyAcids','subcategory':'mycolic','min_completeness':30,'ijsem':'Cell wall contains mycolic acids.'},
        {'id':'gram_positive_profile','name':'Gram-Positive Profile (Firmicutes-type)','min_completeness':30,'ijsem':'Typical Gram-Positive cell wall and lipid profile.',
         'subcategories': [('cellWall', 'teichoic_acid'), ('cellWall', 'peptidoglycan_A4_alpha'), ('quinones', 'menaquinone_classical')]},
        {'id':'gram_negative_profile','name':'Gram-Negative Profile (Proteobacteria-type)','min_completeness':30,'ijsem':'Typical Gram-Negative cell wall and lipid profile.',
         'subcategories': [('cellWall', 'peptidoglycan_A1_gamma'), ('quinones', 'ubiquinone'), ('polarLipids', 'phosphatidylethanolamine_PE')]}
    ]
    predictions = []
    
    for rule in rules:
        matched_markers = []
        sources = []
        total_completeness = 0
        num_subcats = 0
        
        if 'subcategory' in rule: # Single-category rule
            completeness = get_pathway_completeness(found_markers, rule['category'], rule['subcategory'])
            if completeness >= rule['min_completeness']:
                subcat_data = found_markers[rule['category']]['subcategories'][rule['subcategory']]
                matched_markers = subcat_data['found']
                sources = subcat_data['sources']
                total_completeness = completeness
                num_subcats = 1
        
        elif 'subcategories' in rule: # Multi-category rule
            for cat, subcat in rule['subcategories']:
                completeness = get_pathway_completeness(found_markers, cat, subcat)
                if completeness > 0:
                    total_completeness += completeness
                    subcat_data = found_markers[cat]['subcategories'][subcat]
                    matched_markers.extend(subcat_data['found'])
                    sources.extend(subcat_data['sources'])
                num_subcats += 1
        
        if num_subcats > 0:
            avg_completeness = total_completeness / num_subcats
            if avg_completeness >= rule['min_completeness']:
                predictions.append({
                    **rule,
                    'matchedMarkers': list(set(matched_markers)),
                    'confidence': round(avg_completeness, 1),
                    'sources': sources
                })
                
    return sorted(predictions, key=lambda x: x['confidence'], reverse=True)

def identify_novel_features(found_markers):
    novel, complete = [], []
    for category, cat_data in found_markers.items():
        for subcat, subcat_data in cat_data['subcategories'].items():
            if float(subcat_data['percentage']) >= 70:
                complete.append({'category': cat_data['name'],'subcategory': subcat,'percentage': subcat_data['percentage'],'markers': len(subcat_data['found'])})
    if complete:
        novel.append({'type': 'Highly Complete Pathways','description': 'Pathways with ≥70% marker detection','features': complete})
    return novel

def _normalize_phylum_param(raw_val: Optional[str]) -> Optional[List[str]]:
    if not raw_val: return None
    vals = [v.strip() for v in raw_val.split(',') if v.strip()]
    synonyms = {'actinobacteriota':'Actinobacteria','pseudomonadota':'Proteobacteria','bacillota':'Firmicutes','bacteroidota':'Bacteroidetes'}
    mapped = []
    for v in vals:
        v_lower = v.lower()
        mapped.append(synonyms.get(v_lower, v_lower.capitalize()))
    return list(set(mapped)) or None

# --- Main Analysis Function ---
def _analyze_uploaded(files, formats, phylum_param: Optional[str]):
    parser = AnnotationParser()
    detector = MarkerDetector()
    all_genes, file_info = [], []
    for file, format_type in zip(files, formats):
        try:
            content = file.read().decode('utf-8', errors='replace')
        except Exception as e:
            return {'error': f"Error reading file {file.filename}: {e}"}
            
        filename = file.filename
        genes = []
        if format_type in ['gff3','prokka']:
            genes = parser.parse_gff3(content, format_type)
        elif format_type == 'rast':
            genes = parser.parse_rast(content)
        elif format_type == 'bakta':
            genes = parser.parse_bakta_tsv(content)
        elif format_type == 'dram':
            genes = parser.parse_dram(content)
        elif format_type == 'hmm':
            genes = parser.parse_hmmer(content)
        elif format_type == 'gbk':
            genes = parser.parse_genbank(content)
        elif format_type == 'tsv':
            genes = parser.parse_generic_tsv(content)
        else:
            raise ValueError(f'Unsupported format: {format_type}')
        all_genes.extend(genes)
        file_info.append({'name': filename, 'format': format_type, 'geneCount': len(genes)})
    
    if not all_genes:
        return {'error': "No features could be parsed from the uploaded files. Please check file formats."}

    allowed_taxa = _normalize_phylum_param(phylum_param)
    found_markers = detector.find_markers(all_genes, allowed_taxa=allowed_taxa)
    predictions = generate_predictions(found_markers)
    novel_features = identify_novel_features(found_markers)
    
    result = {
        'files': file_info,'totalGenes': len(all_genes),
        'foundMarkers': found_markers,'predictions': predictions,
        'novelFeatures': novel_features,'timestamp': datetime.now().isoformat(),
        'phylum': allowed_taxa
    }
    return result

# --- Report Generation (Template & AI) ---

def _build_template_report(result: dict) -> str:
    lines = []
    lines.append(f"# Genome-resolved Chemotaxonomy Report")
    lines.append(f"_Generated: {datetime.fromisoformat(result.get('timestamp','')).strftime('%Y-%m-%d %H:%M')}_")
    if result.get('phylum'): lines.append(f"**Scope (phylum):** {', '.join(result['phylum'])}")
    lines.append("\n---\n")
    
    lines.append("## 1. Introduction")
    lines.append("This report summarizes the chemotaxonomic features predicted from the provided genome annotations. The analysis searches for key genes and pathways involved in the biosynthesis of taxonomically significant compounds, such as respiratory quinones, fatty acids, polar lipids, and cell wall components.")
    
    lines.append("\n## 2. Methods")
    lines.append("The analysis was performed by searching the parsed genomic features against a curated database of marker genes, EC numbers, and keywords. The following files were processed:")
    for f in result['files']:
        lines.append(f"- **{f['name']}** ({f['format']}): {f['geneCount']} features")
    lines.append(f"\nA total of **{result['totalGenes']} features** were parsed and analyzed.")
    
    lines.append("\n## 3. Results & Discussion")
    lines.append("\n### 3.1. Chemotaxonomic Predictions")
    if result['predictions']:
        for p in result['predictions']:
            lines.append(f"- **{p['name']}**")
            lines.append(f"  - **Confidence:** {p['confidence']}%")
            lines.append(f"  - **IJSEM Note:** {p['ijsem']}")
            lines.append(f"  - **Key Matched Markers:** {', '.join(list(set(p['matchedMarkers']))[:5])}...")
    else:
        lines.append("- No high-confidence chemotaxonomic profiles were matched.")
    
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
    lines.append("This in-silico analysis provides a preliminary chemotaxonomic profile. The predictions and marker findings should be validated with laboratory-based chemical analysis (e.g., HPLC, mass spectrometry).")
    
    return "\n".join(lines)

def _generate_ai_report(result: dict) -> str:
    # This is a fallback in case the AI call fails
    fallback_report = _build_template_report(result)
    
    try:
        # NOTE: This uses the Gemini 2.5 Flash model. No API key is needed when run in this environment.
        api_url = "https://generativelanguage.googleapis.com/v1beta/models/gemini-2.5-flash-preview-09-2025:generateContent?key="
        
        system_prompt = (
            "You are a helpful assistant specialized in microbial taxonomy and chemotaxonomy. "
            "Your task is to write a professional, narrative report based on JSON data provided by a user. "
            "The report should be structured with Markdown headings (##) for 'Introduction', 'Methods', 'Results and Discussion', and 'Conclusion'. "
            "In the 'Results and Discussion' section, interpret the findings. Explain what the predictions (e.g., 'Gram-Negative Profile') mean. "
            "Discuss the completeness of the pathways (e.g., 'The high completeness of the ubiquinone pathway...'). "
            "Do not just list the data; synthesize it into a coherent scientific summary. Be formal and concise."
        )
        
        # Prepare a simpler version of the results for the prompt
        prompt_data = {
            'inputs': result.get('files', []),
            'phylum_filter': result.get('phylum', 'None'),
            'predictions': result.get('predictions', []),
            'marker_findings': result.get('foundMarkers', {}),
            'notable_features': result.get('novelFeatures', [])
        }
        
        user_prompt = f"Please generate a full chemotaxonomy report based on the following JSON data:\n\n{json.dumps(prompt_data, indent=2)}"

        payload = {
            "contents": [{"parts": [{"text": user_prompt}]}],
            "systemInstruction": {"parts": [{"text": system_prompt}]}
        }
        
        headers = {"Content-Type": "application/json"}
        
        # Retry logic for API calls
        response = None
        for attempt in range(3):
            try:
                response = requests.post(api_url, json=payload, headers=headers, timeout=60)
                if response.status_code == 200:
                    break
            except requests.exceptions.RequestException as e:
                print(f"AI report generation attempt {attempt+1} failed: {e}")
                if attempt == 2:
                    raise # Re-raise on last attempt
        
        if not response or response.status_code != 200:
            print(f"AI report generation failed with status {response.status_code}: {response.text}")
            return fallback_report
            
        data = response.json()
        text = data.get('candidates', [{}])[0].get('content', {}).get('parts', [{}])[0].get('text', '')
        
        if not text:
            print("AI report generation returned empty text.")
            return fallback_report
            
        return text

    except Exception as e:
        print(f"Error during AI report generation: {e}")
        return fallback_report # Return template report on any failure


# --- API Endpoints ---

@app.route('/')
def index():
    # This serves the index.html file from the same directory
    try:
        with open('index.html', 'r') as f:
            return render_template_string(f.read())
    except FileNotFoundError:
        return "index.html not found", 404

@app.route('/api/analyze-multi', methods=['POST'])
def analyze_multiple_annotations():
    if 'files' not in request.files:
        return jsonify({'error': 'No files uploaded'}), 400
    files = request.files.getlist('files')
    formats = request.form.getlist('formats')
    phylum_param = request.form.get('phylum')
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
    phylum_param = request.form.get('phylum')
    
    # Get the requested format (pdf, docx)
    req_format = request.args.get('format', 'pdf')
    if req_format not in ['pdf', 'docx']:
        return jsonify({'error': 'Invalid format requested. Must be pdf or docx.'}), 400

    # Get the AI toggle state
    use_ai = request.form.get('use_ai', 'false').lower() == 'true'

    if len(files) != len(formats):
        return jsonify({'error': 'Number of files and formats must match'}), 400
    
    # Run the analysis
    result = _analyze_uploaded(files, formats, phylum_param)
    if 'error' in result:
        return jsonify(result), 400

    # Generate the Markdown report string
    print(f"Generating report. AI enabled: {use_ai}")
    if use_ai:
        markdown_string = _generate_ai_report(result)
    else:
        markdown_string = _build_template_report(result)
    print("Markdown generation complete.")

    # --- NEW: Pandoc Conversion ---
    md_file = None
    out_file = None
    try:
        # 1. Create a temporary markdown file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as md_file:
            md_file.write(markdown_string)
            md_file_path = md_file.name
        
        # 2. Create a temporary output file
        with tempfile.NamedTemporaryFile(suffix=f'.{req_format}', delete=False) as out_file:
            out_file_path = out_file.name

        # 3. Run Pandoc
        print(f"Converting Markdown to {req_format} using Pandoc...")
        pandoc_args = []
        if req_format == 'pdf':
            pandoc_args = ['--pdf-engine=pdflatex'] # Use pdflatex for PDF generation
            
        pypandoc.convert_file(
            md_file_path,
            req_format,
            outputfile=out_file_path,
            extra_args=pandoc_args
        )
        print("Pandoc conversion successful.")

        # 4. Prepare file for sending
        mimetype_map = {
            'pdf': 'application/pdf',
            'docx': 'application/vnd.openxmlformats-officedocument.wordprocessingml.document'
        }
        
        fname = f"chemotax_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.{req_format}"
        
        return send_file(
            out_file_path,
            mimetype=mimetype_map[req_format],
            as_attachment=True,
            download_name=fname
        )
        
    except Exception as e:
        print(f"Error during file conversion or sending: {e}")
        return jsonify({'error': f'Failed to generate report file: {e}'}), 500
    finally:
        # 5. Clean up temporary files
        if md_file and os.path.exists(md_file.name):
            os.remove(md_file.name)
        if out_file and os.path.exists(out_file.name):
            os.remove(out_file.name)
        print("Temporary files cleaned up.")


@app.route('/health', methods=['GET'])
def health_check():
    return jsonify({'status': 'healthy', 'timestamp': datetime.now().isoformat()})

if __name__ == '__main__':
    # This is for local development only
    # Render uses the CMD in the Dockerfile
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 5000)))
