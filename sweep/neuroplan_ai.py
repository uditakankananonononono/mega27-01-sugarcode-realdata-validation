"""Sweep: neuroplan_ai vs admissible-A* oracle (same cost), exact corridor/NeuroTwin/segmentation/tractography recomputation."""
import json, math, os, sys, heapq
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from sugarcode.modules.neuroplan_ai import core as N
res = {'module': 'neuroplan_ai', 'sources': ['module-internal oracle (admissible A*, same cost function)', 'numpy/scipy recomputation']}

# 1. A* optimality: module uses 50x inflated heuristic but docstring claims admissible.
# Oracle: identical step cost with weight-1 euclidean heuristic (admissible: cost >= step).
def oracle(tumor, entry, image_size=(64,64,64), risk_lambda=8.0):
    center = [float(c) for c in tumor['center']]
    regions = N._place_regions(image_size, [int(c) for c in center], tumor.get('radius_mm', 20))
    risk = N._risk_field_fn(regions)
    start = tuple(int(round(e)) for e in entry); goal = tuple(int(round(c)) for c in center)
    neigh = [(dx,dy,dz) for dx in (-1,0,1) for dy in (-1,0,1) for dz in (-1,0,1) if (dx,dy,dz)!=(0,0,0)]
    dist = {start: 0.0}; pq = [(0.0, start)]
    while pq:
        g, p = heapq.heappop(pq)
        if g > dist.get(p, float('inf')): continue
        if p == goal: return g
        for d in neigh:
            q = (p[0]+d[0], p[1]+d[1], p[2]+d[2])
            if not all(0 <= q[k] < image_size[k] for k in range(3)): continue
            step = math.sqrt(d[0]**2+d[1]**2+d[2]**2)
            mid = tuple((p[k]+q[k])/2 for k in range(3))
            ng = g + step * (1.0 + risk_lambda * min(risk(mid), 1.25))
            if ng < dist.get(q, float('inf')):
                dist[q] = ng; heapq.heappush(pq, (ng, q))
    return dist.get(goal, float('inf'))
def path_cost(path, risk, risk_lambda=8.0):
    tot = 0.0
    for i in range(len(path)-1):
        p, q = path[i], path[i+1]
        step = math.sqrt(sum((p[k]-q[k])**2 for k in range(3)))
        mid = tuple((p[k]+q[k])/2 for k in range(3))
        tot += step * (1.0 + risk_lambda * min(risk(mid), 1.25))
    return tot
tumor = {'center': [32, 32, 18], 'radius_mm': 8, 'type': 'glioma'}
entry = [10, 50, 63]  # oblique entry so the risk field matters
r = N.plan_path_astar(tumor, entry=entry, image_size=(64,64,64), max_expand=400000)
regions = N._place_regions((64,64,64), tumor['center'], 8)
risk = N._risk_field_fn(regions)
opt = oracle(tumor, entry)
got = path_cost(r['path'], risk)
cont = all(max(abs(r['path'][i+1][k]-r['path'][i][k]) for k in range(3)) == 1 for i in range(len(r['path'])-1))
res['astar_optimality'] = {'oracle_cost': round(opt, 4), 'module_path_cost': round(got, 4),
    'optimality_gap_pct': round(100*(got-opt)/opt, 2), 'path_continuous_26conn': cont,
    'docstring_claim': 'euclidean heuristic (admissible since risk term is >= 1)',
    'code_reality': 'frontier priority is g + 50*h: inflated 50x, greedy weighted A*, not admissible'}
# 2. corridor math recomputation
plan = N.plan_surgery(tumor, (64,64,64))
ok_c = []
for c in plan['corridor_options']:
    d = [x for x in c['direction']]; norm = math.sqrt(sum(x*x for x in d))
    risk_sum = 0.0
    for rn, rg in regions.items():
        rel = [rg['position'][i]-tumor['center'][i] for i in range(3)]
        proj = sum(rel[i]*d[i]/norm for i in range(3)) if abs(norm-1) > 1e-9 else sum(rel[i]*d[i] for i in range(3))
        dn = [x/norm for x in d]
        proj = sum(rel[i]*dn[i] for i in range(3))
        closest = [proj*dn[i] for i in range(3)]
        perp = math.sqrt(max(0.0, sum((rel[i]-closest[i])**2 for i in range(3))))
        if proj > 0: risk_sum += rg['risk_weight']/(1+0.1*perp)
    ok_c.append(abs(risk_sum - c['risk']) < 1e-3)
res['corridor_recompute'] = {'all_match': all(ok_c), 'n': len(ok_c)}
# 3. simulate_neurotwin recomputation
twin = N.simulate_neurotwin(tumor, plan['recommended_corridor'], regions, resection_fraction=.8)
ok_t = []
for d in twin['functional_risks']:
    rg = regions[d['region']]
    clearance = max(0, rg['distance_to_tumor_mm'] - 8 - 2.0)
    cf = 1.5 if plan['recommended_corridor'].get('nearest_eloquent') == d['region'] else 1
    p = min(.95, rg['risk_weight']*cf*.8/(1+clearance/5))
    ok_t.append(abs(p - d['probability']) < 1e-6)
agg = 1 - math.prod(1-d['probability'] for d in twin['functional_risks'])
res['neurotwin_recompute'] = {'probabilities_match': all(ok_t), 'aggregate_ok': abs(agg - twin['aggregate_deficit_risk']) < 1e-6,
                              'residual_ok': abs(twin['residual_tumor_volume_mm3'] - .2*4/3*math.pi*8**3) < 1e-3}
# 4. segment_mri on a synthetic bright sphere
rng = np.random.default_rng(0)
vol = rng.normal(10, 1, (32,32,32))
yy, xx, zz = np.mgrid[0:32,0:32,0:32]
sphere = (xx-16)**2 + (yy-16)**2 + (zz-16)**2 <= 6.0**2
vol[sphere] += 30
seg = N.segment_mri(vol, voxel_size_mm=(2.0,2.0,2.0))
res['segment_mri'] = {'threshold': seg['threshold'], 'threshold_recomputed': float(vol.mean()+1.5*vol.std()),
    'voxels': seg['tumor_voxels'], 'sphere_voxels': int(sphere.sum()),
    'voxel_ratio': round(seg['tumor_voxels']/sphere.sum(), 3), 'volume_ok': abs(seg['tumor_volume_mm3'] - seg['tumor_voxels']*8) < 1e-6,
    'centroid_mm': [round(c,1) for c in seg['centroid_mm']]}
# 5. analyze_tractography recomputation
sl = [[[0,0,0],[10,0,0],[20,0,0]], [[5,9,0],[15,9,0],[25,9,0]]]
tr = N.analyze_tractography(sl, [10,0,0], 3.0, critical_distance_mm=5.0)
exp = [ (0, -3.0, True), (9.0, 6.0, False) ]
res['tractography'] = {'rows_post_fix': [(x['length_mm'], x['clearance_mm'], x['at_risk'], x['intersects_tumor']) for x in tr['streamlines']],
                       'expected_clearances': [-3.0, 6.0],
                       'pre_fix': 'vertex-only distance: streamline 2 clearance 7.29563 (BUG 47); segment-based fix gives 6.0'}
# 6. atlas anatomy ordering sanity
pos = {k: v['position'] for k, v in N._place_regions((128,128,128), [64,64,64], 10).items()}
res['atlas_ordering'] = {'brainstem_lowest_z': pos['brainstem'][2] == min(p[2] for p in pos.values()),
                         'motor_most_superior_z': pos['motor_cortex'][2] == max(p[2] for p in pos.values()),
                         'visual_most_posterior_y': pos['visual_cortex'][1] == min(p[1] for p in pos.values())}
fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_neuroplan_ai.json')
json.dump(res, open(fn, 'w'), indent=1, default=str)
print(json.dumps(res, default=str))
