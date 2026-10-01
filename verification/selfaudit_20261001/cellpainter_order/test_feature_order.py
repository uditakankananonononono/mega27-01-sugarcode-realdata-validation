"""Validation regression expected to fail at product 1ccbd030."""
from sugarcode.modules.cellpainter import nearest_mechanism, profile_perturbation

def test_feature_mapping_insertion_order_does_not_change_similarity():
    query = profile_perturbation('dna_damage')
    reference = profile_perturbation('dna_damage')
    reordered = {**reference, 'signature_vector': dict(reversed(list(reference['signature_vector'].items())))}
    assert reference['signature_vector'] == reordered['signature_vector']
    original = nearest_mechanism(query, [reference])[0]['similarity']
    changed = nearest_mechanism(query, [reordered])[0]['similarity']
    assert abs(original - changed) < 1e-12
