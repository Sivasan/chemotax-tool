import pytest
from backend_api import AnnotationParser

def test_parse_hmmer_domtblout():
    content = """#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
# ------------------ ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
target1        -            100 queryA      PF00001      200   1.2e-30  100.5   0.1   1   1   1.2e-30   1.2e-30  100.5   0.1     1   200     1   100     1   100 0.99 description of target1
"""
    genes = AnnotationParser.parse_hmmer(content)
    assert len(genes) == 1
    g = genes[0]
    assert g['target_name'] == 'target1'
    assert g['query_name'] == 'queryA'
    assert g['evalue'] == 1.2e-30
    assert g['product'] == 'description of target1'

def test_parse_hmmer_tblout():
    # Valid tblout line with all columns
    content = """#                                                               --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   query name           accession   E-value  score  bias   E-value  score  bias   exp  N  Model    Description
# ------------------ ----------  -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --  -------- -----------
target2        -           queryB      PF00002     2.5e-10  50.0   0.0   2.5e-10  50.0   0.0   1.0  1  0  0  0  0  0  0  description of target2
"""
    genes = AnnotationParser.parse_hmmer(content)
    assert len(genes) == 1
    g = genes[0]
    assert g['target_name'] == 'target2'
    assert g['query_name'] == 'queryB'
    assert g['evalue'] == 2.5e-10
    assert g['product'] == 'description of target2'
