"""Validation regression expected to fail at product 1ccbd030.

Run with pinned product src/ in PYTHONPATH and fastapi/httpx installed.
This is an in-process test, with no external request or account action.
"""
from fastapi.testclient import TestClient
from omega.api import app
from sugarcode.modules.ecosystem import generate_sdk_client

def test_generated_sdk_route_exists_in_shipped_app():
    client = TestClient(app)
    assert client.get('/modules/crispr_opt').status_code == 200
    sdk = generate_sdk_client('crispr_opt', 'design_guides')
    assert client.post(sdk['path'], json={'sequence': 'ATGC'}).status_code != 404
