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
