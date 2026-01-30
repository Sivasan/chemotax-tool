import pytest
from backend_api import AnnotationParser

def test_parse_fasta_ncbi_style():
    content = """>WP_12345.1 [gene=menF] [protein=menaquinone biosynthesis protein MenF] [EC_number=1.2.3.4]
MKLTV...
"""
    genes = AnnotationParser.parse_fasta(content)
    assert len(genes) == 1
    g = genes[0]
    assert g['id'] == 'WP_12345.1'
    assert g['gene'] == 'menF'
    assert g['product'] == 'menaquinone biosynthesis protein MenF'
    assert g['ec_number'] == '1.2.3.4'

def test_parse_fasta_prokka_style():
    content = """>seq1 product=Isochorismate synthase MenF;EC_number=5.4.4.2
MKL...
"""
    genes = AnnotationParser.parse_fasta(content)
    assert len(genes) == 1
    g = genes[0]
    assert g['id'] == 'seq1'
    # Prokka often doesn't have gene=, just product.
    # Our detector relies on product name too.
    assert g['product'] == 'Isochorismate synthase MenF'
    assert g['ec_number'] == '5.4.4.2'

def test_parse_fasta_simple():
    content = """>seq2 menaquinone biosynthesis protein
MKL...
"""
    genes = AnnotationParser.parse_fasta(content)
    assert len(genes) == 1
    g = genes[0]
    # Should use description as product if no tags
    assert 'menaquinone biosynthesis protein' in g['product']
