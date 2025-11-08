# backend_api.py - Flask Backend for Multi-Annotation Chemotaxonomy Platform
from flask import Flask, request, jsonify, send_file, send_from_directory
from flask_cors import CORS
from Bio import SeqIO
import io, re, csv
from datetime import datetime
from collections import defaultdict
from typing import List, Optional

app = Flask(__name__)
# CORS is still good to have, just in case.
CORS(app, resources={r"/api/*": {"origins": "*"}}, expose_headers=["Content-Disposition"])

class AnnotationParser:
    """Parser for multiple genome annotation formats"""
    @staticmethod
    def parse_gff3(content, source_name='gff3'):
        genes = []
        for line in content.split('\n'):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9: continue
            attributes = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key.strip()] = value.strip()
            gene_info = {
                'seqid': fields[0],
                'source': fields[1],
                'type': fields[2],
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
        if not dbxref: return ''
        m = re.search(r'EC:(\d+\.\d+\.\d+\.\d+)', dbxref)
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
                'gene_id': row.get('gene_id', row.get('fasta', '')),
                'scaffold': row.get('scaffold', ''),
                'gene_position': row.get('gene_position', ''),
                'start': row.get('start_position', ''),
                'end': row.get('end_position', ''),
                'strandedness': row.get('strandedness', ''),
                'product': row.get('product', row.get('kegg_hit', row.get('pfam_hits', ''))),
                'gene': row.get('gene_name', ''),
                'ec_number': row.get('ec_number', row.get('kegg_ec', '')),
                'kegg_id': row.get('kegg_id', ''),
                'ko_id': row.get('ko_id', ''),
                'pfam_list': pfams,
                'cazy': row.get('cazy_hits', ''),
                'format': 'dram'
            })
        return genes

    @staticmethod
    def parse_hmmer(content):
        genes = []
        for line in content.split('\n'):
            if line.startswith('#') or not line.strip(): continue
            
            # Use re.split with maxsplit=22 to keep the description intact
            fields = re.split(r'\s+', line.strip(), maxsplit=22)
            
            if len(fields) < 23: continue 

            target_acc = fields[1]
            pfam_acc = target_acc if target_acc.startswith('PF') else ''
            
            try:
                evalue = float(fields[6])
                score = float(fields[7])
                bias = float(fields[8])
                domain_evalue = float(fields[12])
                domain_score = float(fields[13])
            except ValueError:
                continue

            genes.append({
                'target_name': fields[0],
                'target_accession': target_acc,
                'query_name': fields[3],
                'query_accession': fields[4],
                'evalue': evalue,
                'score': score,
                'bias': bias,
                'domain_evalue': domain_evalue,
                'domain_score': domain_score,
                'product': fields[22], # field 22 is the full description
                'gene': fields[0], # Use target_name as gene identifier
                'pfam': pfam_acc,
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
                        genes.append({
                            'type': feature.type,
                            'location': str(feature.location),
                            'start': int(feature.location.start),
                            'end': int(feature.location.end),
                            'strand': feature.location.strand,
                            'gene': (q.get('gene', ['']) or [''])[0],
                            'product': (q.get('product', ['']) or [''])[0],
                            'locus_tag': (q.get('locus_tag', ['']) or [''])[0],
                            'protein_id': (q.get('protein_id', ['']) or [''])[0],
                            'ec_number': (q.get('EC_number', ['']) or [''])[0],
                            'note': (q.get('note', ['']) or [''])[0],
                            'translation': (q.get('translation', ['']) or [''])[0],
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
                    if std not in gene_info:
                        gene_info[std] = ''
                
                gene_info.update(row)
                genes.append(gene_info)
        except Exception as e:
            print(f"Error parsing TSV: {e}")
        return genes

class MarkerDetector:
    def __init__(self):
        self.marker_database = self._load_marker_database()
        self.alias_to_canon, self.canon_to_aliases = self._build_alias_maps()

    def _build_alias_maps(self):
        alias_table = {
            'meng': {'ubie'},
            'ubie': {'meng'},
            'idsa': {'ddsa'},
            'ddsa': {'idsa'},
            'fbpa': {'antigen85a', 'ag85a'},
            'fbpb': {'antigen85b', 'ag85b'},
            'fbpc': {'antigen85c', 'ag85c'},
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
                    'genes': [
                        'pks13','fadD32','accD4','acpM','kasA','kasB','inhA','mabA','hadA','hadB','hadC',
                        'mmaA1','mmaA2','mmaA3','mmaA4','cmaA1','cmaA2','pcaA','fbpA','fbpB','fbpC','mmpL3'
                    ],
                    'ec': ['1.3.1.9','1.1.1.100','2.3.1.41','2.3.1.-'],
                    'keywords': [
                        'mycolic acid','FAS-II','trehalose monomycolate','trehalose dimycolate',
                        'antigen 85','mycoloyltransferase','MmpL3 transporter','condensing enzyme'
                    ],
                    'taxa': ['Actinobacteria','Actinobacteriota','Corynebacteriales']
                }
            },
            'quinones': {
                'name': 'Respiratory Quinones',
                'menaquinone_classical': {
                    'genes': ['menF','menD','menH','menC','menE','menB','menI','menA','menG','ubiE','menJ'],
                    'ec': ['2.5.1.74','6.2.1.26'],
                    'keywords': ['menaquinone','vitamin K2','demethylmenaquinone methyltransferase','1,4-dihydroxy-2-naphthoate'],
                    'taxa': ['Firmicutes','Bacillota','Actinobacteria','Actinobacteriota','Bacteroidetes','Bacteroidota','Proteobacteria','Pseudomonadota']
                },
                'menaquinone_futalosine': {
                    'genes': ['mqnA','mqnB','mqnC','mqnD','mqnE','mqnP','mqnL','mqnX'],
                    'ec': [],
                    'keywords': ['futalosine pathway','aminofutalosine','MK biosynthesis futalosine'],
                    'taxa': ['Actinobacteria','Bacteroidetes','Thermus','Campylobacter','Helicobacter']
                },
                'ubiquinone': {
                    'genes': ['ubiA','ubiB','ubiC','ubiD','ubiX','ubiE','ubiF','ubiG','ubiH','ubiI','ubiJ','ubiK','ubiU','ubiV','ubiT','coq7'],
                    'ec': ['2.5.1.39'],
                    'keywords': ['ubiquinone','coenzyme Q','4-hydroxybenzoate','O2-independent ubiquinone','anaerobic UQ biosynthesis'],
                    'taxa': ['Proteobacteria','Pseudomonadota','Acidobacteria','Chlorobi']
                },
                'polyprenyl_side_chain': {
                    'genes': ['ispB','hepS','hepT','idsA','ddsA'],
                    'ec': [],
                    'keywords': ['octaprenyl diphosphate synthase','heptaprenyl diphosphate synthase','polyprenyl diphosphate synthase'],
                    'taxa': ['Proteobacteria','Firmicutes','Bacillota','Actinobacteria','Bacteroidetes']
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
                tokens = re.findall(r'\w+', product.lower())
                for token in tokens:
                    if len(token) > 3:
                        product_index[token].append((idx, gene))
                product_fulltext.append((idx, product.lower()))
            ec = gene.get('ec_number', '')
            if ec:
                for ec_num in re.findall(r'\d+\.\d+\.\d+\.\d+', ec):
                    ec_index[ec_num].append((idx, gene))
            if gene.get('pfam'):
                m = re.match(r'(PF\d{5})', str(gene['pfam']))
                if m: pfam_index[m.group(1)].append((idx, gene))
            if gene.get('pfam_list'):
                for pf in gene['pfam_list']:
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
                found_gene_indices = set()

                for marker_gene in marker_data.get('genes', []):
                    mg = self._clean_token(marker_gene)
                    names_to_try = {mg}
                    names_to_try |= self.canon_to_aliases.get(mg, set())
                    if mg in self.alias_to_canon:
                        names_to_try.add(self.alias_to_canon[mg])
                        names_to_try |= self.canon_to_aliases.get(self.alias_to_canon[mg], set())
                    
                    for name_try in names_to_try:
                        matches = gene_index.get(name_try, [])
                        for idx, gene in matches:
                            if idx not in found_gene_indices:
                                found_set.add(marker_gene)
                                sources.append({
                                    'marker': marker_gene,
                                    'source': gene.get('format', 'unknown'),
                                    'gene': gene.get('gene') or gene.get('gene_norm', ''),
                                    'product': gene.get('product', '')
                                })
                                found_gene_indices.add(idx)

                for ec_num in marker_data.get('ec', []):
                    matches = _ec_matches_query(ec_num)
                    for idx, gene in matches:
                        if idx not in found_gene_indices:
                            marker_id = f"EC:{ec_num}"
                            found_set.add(marker_id)
                            sources.append({
                                'marker': marker_id,
                                'source': gene.get('format', 'unknown'),
                                'ec': gene.get('ec_number',''),
                                'product': gene.get('product', '')
                            })
                            found_gene_indices.add(idx)

                for keyword in marker_data.get('keywords', []):
                    kw = keyword.lower()
                    tokens_to_check = set()
                    for t in kw.split():
                        if len(t) > 3:
                            tokens_to_check.add(t)
                    
                    for token in tokens_to_check:
                         matches = product_index.get(token, [])
                         for idx, gene in matches:
                            if kw in gene.get('product','').lower() and idx not in found_gene_indices:
                                marker_id = f"kw:{keyword}"
                                found_set.add(marker_id)
                                sources.append({
                                    'marker': marker_id,
                                    'source': gene.get('format', 'unknown'),
                                    'keyword': keyword,
                                    'product': gene.get('product', '')
                                })
                                found_gene_indices.add(idx)

                    for idx, full in product_fulltext:
                        if kw in full and idx not in found_gene_indices:
                            gene = all_genes[idx]
                            marker_id = f"kw:{keyword}"
                            found_set.add(marker_id)
                            sources.append({
                                'marker': marker_id,
                                'source': gene.get('format', 'unknown'),
                                'keyword': keyword,
                                'product': gene.get('product', '')
                            })
                            found_gene_indices.add(idx)

                for pf in marker_data.get('pfam', []):
                    matches = pfam_index.get(pf, [])
                    for idx, gene in matches:
                        if idx not in found_gene_indices:
                            marker_id = f"pfam:{pf}"
                            found_set.add(marker_id)
                            sources.append({
                                'marker': marker_id,
                                'source': gene.get('format', 'unknown'),
                                'product': gene.get('product', ''),
                                'pfam': pf
                            })
                            found_gene_indices.add(idx)

                if found_set:
                    total = len(marker_data.get('genes', [])) + len(marker_data.get('ec', [])) + \
                            len(marker_data.get('keywords', [])) + len(marker_data.get('pfam', []))
                    found_markers[category]['subcategories'][subcat] = {
                        'found': list(found_set),
                        'total': total,
                        'percentage': round((len(found_set) / total * 100), 1) if total else 0.0,
                        'sources': sources
                    }

        found_markers = {k:v for k,v in found_markers.items() if v['subcategories']}
        return dict(found_markers)

def generate_predictions(found_markers):
    rules = [
        {'id':'actinobacteria_mycolic','name':'Actinobacteria (Mycolic Acid Type)','category':'fattyAcids','subcategory':'mycolic','min_markers':2,'ijsem':'Cell wall contains mycolic acids'},
        {'id':'firmicutes_menaquinone','name':'Firmicutes Type','category':'quinones','subcategory':'menaquinone_classical','min_markers':3,'ijsem':'Major quinone: Menaquinone MK-7, MK-8, or MK-9'}
    ]
    predictions = []
    for rule in rules:
        if rule['category'] in found_markers:
            subcats = found_markers[rule['category']]['subcategories']
            if rule['subcategory'] in subcats:
                subcat_data = subcats[rule['subcategory']]
                if len(subcat_data['found']) >= rule['min_markers']:
                    confidence = min((len(subcat_data['found']) / subcat_data['total']) * 100, 100)
                    predictions.append({**rule,'matchedMarkers': subcat_data['found'],'confidence': round(confidence,1),'sources': subcat_data['sources']})
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
        mapped.append(synonyms.get(v_lower, v.title()))
    return list(set(mapped)) or None # Return unique values

def _analyze_uploaded(files, formats, phylum_param: Optional[str]):
    parser = AnnotationParser()
    detector = MarkerDetector()
    all_genes, file_info = [], []
    for file, format_type in zip(files, formats):
        content = file.read().decode('utf-8', errors='replace')
        filename = file.filename
        if format_type in ['gff3','prokka']:
            genes = parser.parse_gff3(content, format_type)
        elif format_type == 'bakta':
            if filename.endswith('.gff') or filename.endswith('.gff3'):
                genes = parser.parse_gff3(content, 'bakta')
            else:
                genes = parser.parse_bakta_tsv(content)
        elif format_type == 'rast':
            genes = parser.parse_rast(content)
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
    allowed_taxa = _normalize_phylum_param(phylum_param)
    found_markers = detector.find_markers(all_genes, allowed_taxa=allowed_taxa)
    predictions = generate_predictions(found_markers)
    novel_features = identify_novel_features(found_markers)
    result = {'files': file_info,'totalGenes': len(all_genes),'foundMarkers': found_markers,'predictions': predictions,'novelFeatures': novel_features,'timestamp': datetime.now().isoformat(),'phylum': allowed_taxa}
    return result

def _build_markdown_report(result: dict) -> bytes:
    lines = []
    lines.append(f"# Genome-resolved Chemotaxonomy Report")
    lines.append(f"_Generated: {result.get('timestamp','')}_")
    if result.get('phylum'): lines.append(f"**Scope (phylum):** {', '.join(result['phylum'])}")
    lines.append("")
    lines.append("## Inputs")
    for f in result['files']:
        lines.append(f"- **{f['name']}** ({f['format']}): {f['geneCount']} features")
    lines.append(f"- **Total parsed features:** {result['totalGenes']}")
    lines.append("")
    lines.append("## Predictions")
    if result['predictions']:
        for p in result['predictions']:
            lines.append(f"- **{p['name']}** — confidence {p['confidence']}%")
            lines.append(f"  - IJSEM note: {p['ijsem']}")
            lines.append(f"  - Matched markers: {', '.join(sorted(list(set(p['matchedMarkers']))))}")
    else:
        lines.append("- No rules matched with sufficient confidence.")
    lines.append("")
    lines.append("## Marker Findings (by category)")
    for cat, cat_data in result['foundMarkers'].items():
        lines.append(f"### {cat_data['name']}")
        for sub, sub_data in cat_data['subcategories'].items():
            lines.append(f"- **{sub}** — {sub_data['percentage']}% ({len(sub_data['found'])}/{sub_data['total']})")
    lines.append("")
    lines.append("## Notable Features")
    if result['novelFeatures']:
        for block in result['novelFeatures']:
            lines.append(f"- **{block['type']}** — {block['description']}")
            for feat in block['features']:
                lines.append(f"  - {feat['category']} / {feat['subcategory']}: {feat['percentage']}%")
    else:
        lines.append("- None flagged at current thresholds.")
    lines.append("")
    return ("\n".join(lines)).encode("utf-8")

# --- API ROUTES ---
# These handle the '/api/...' calls from the frontend

@app.route('/api/analyze-multi', methods=['POST'])
def analyze_multiple_annotations():
    if 'files' not in request.files:
        return jsonify({'error': 'No files uploaded'}), 400
    files = request.files.getlist('files')
    formats = request.form.getlist('formats')
    phylum_param = request.form.get('phylum')
    if len(files) != len(formats):
        return jsonify({'error': 'Number of files and formats must match'}), 400
    try:
        result = _analyze_uploaded(files, formats, phylum_param)
        return jsonify(result)
    except Exception as e:
        app.logger.error(f"Error in /api/analyze-multi: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500

@app.route('/api/analyze-and-report', methods=['POST'])
def analyze_and_report():
    if 'files' not in request.files:
        return jsonify({'error': 'No files uploaded'}), 400
    files = request.files.getlist('files')
    formats = request.form.getlist('formats')
    phylum_param = request.form.get('phylum')
    if len(files) != len(formats):
        return jsonify({'error': 'Number of files and formats must match'}), 400
    try:
        result = _analyze_uploaded(files, formats, phylum_param)
        md_bytes = _build_markdown_report(result)
        fname = f"chemotax_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
        return send_file(io.BytesIO(md_bytes), mimetype='text/markdown; charset=utf-8', as_attachment=True, download_name=fname)
    except Exception as e:
        app.logger.error(f"Error in /api/analyze-and-report: {e}", exc_info=True)
        return jsonify({'error': str(e)}), 500

@app.route('/health', methods=['GET'])
def health_check():
    return jsonify({'status': 'healthy', 'timestamp': datetime.now().isoformat()})

# --- NEW FRONTEND ROUTE ---
# This new route will serve your index.html file when
# someone visits the main URL (e.g., "[my-tool.onrender.com/](https://my-tool.onrender.com/)")

@app.route('/')
def serve_index():
    # 'send_from_directory' serves files from a directory.
    # '.' means the current directory (where backend_api.py is)
    return send_from_directory('.', 'index.html')

if __name__ == '__main__':
    # This now only runs if you execute 'python backend_api.py' locally
    app.run(debug=True, host='0.0.0.0', port=5000)
