import pytest
from backend_api import app
import io

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_analyze_multi_gff3(client):
    data = {
        'files': (io.BytesIO(b"##gff-version 3\nseq1\tprokka\tCDS\t1\t100\t.\t+\t0\tID=gene1;product=test product;gene=menF"), 'test.gff3'),
        'formats': 'gff3'
    }
    rv = client.post('/api/analyze-multi', data=data, content_type='multipart/form-data')
    assert rv.status_code == 200
    json_data = rv.get_json()
    assert len(json_data['files']) == 1
    assert json_data['totalGenes'] == 1
    # menF should be found in quinones
    assert 'menF' in json_data['foundMarkers']['quinones']['subcategories']['menaquinone_classical']['found']

def test_analyze_comparative(client):
    # File 1 has menF
    file1 = (io.BytesIO(b">seq1 [gene=menF] [protein=MenF]\nATG..."), 'file1.faa')
    # File 2 has ubiE
    file2 = (io.BytesIO(b">seq2 [gene=ubiE] [protein=UbiE]\nATG..."), 'file2.faa')

    data = {
        'files': [file1, file2],
        'formats': ['fasta', 'fasta']
    }
    rv = client.post('/api/analyze-multi', data=data, content_type='multipart/form-data')
    assert rv.status_code == 200
    res = rv.get_json()

    # Combined should have both
    markers = []
    for cat in res['foundMarkers'].values():
        for sub in cat['subcategories'].values():
            markers.extend(sub['found'])

    assert 'menF' in markers
    # ubiE is an alias for menG, or sometimes separate. Let's check logic.
    # In MarkerDetector, menG aliases include ubiE. So 'menG' might be the marker name found.
    # Let's check 'individual' stats

    stats = res['stats']
    assert stats['combined_count'] >= 2

    # Check individual counts
    f1_stats = next(x for x in stats['individual'] if x['name'] == 'file1.faa')
    f2_stats = next(x for x in stats['individual'] if x['name'] == 'file2.faa')

    assert f1_stats['marker_count'] >= 1
    assert f2_stats['marker_count'] >= 1

    # Cumulative gain should be positive (combining gives more than max of single)
    # If file1 has 1, file2 has 1, and they are different -> combined is 2. Max single is 1. Gain is 1.
    assert stats['gain'] >= 1

def test_analyze_and_report_pdf(client):
    data = {
        'files': (io.BytesIO(b"##gff-version 3\nseq1\tprokka\tCDS\t1\t100\t.\t+\t0\tID=gene1;product=test product;gene=menF"), 'test.gff3'),
        'formats': 'gff3',
        'use_ai': 'false'
    }

    # Passing query param format=pdf
    rv = client.post('/api/analyze-and-report?format=pdf', data=data, content_type='multipart/form-data')
    assert rv.status_code == 200
    assert rv.headers['Content-Type'] == 'application/pdf'
    assert len(rv.data) > 0
