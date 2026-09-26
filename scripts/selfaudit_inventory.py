#!/usr/bin/env python3
"""Sugarcode self-audit inventory v2: import-based module <-> test mapping.
Per-module pass/fail attribution uses the test file's own imports of
sugarcode.modules.<name>; test files importing no module (meta gates, CLI,
bio toolkit) are reported separately. Honest counts, no inference."""
import ast, json, pathlib, re, collections
root = pathlib.Path('/tmp/sc_audit')
mods = sorted(p.name for p in (root/'src/sugarcode/modules').iterdir() if p.is_dir() and not p.name.startswith('__'))
specs = {p.stem for p in (root/'spec').glob('*.md')}
log = (root/'pytest_full.log').read_text()
file_counts = collections.defaultdict(lambda: {'passed':0,'failed':0,'skipped':0})
for line in log.splitlines():
    m = re.match(r'(PASSED|FAILED|SKIPPED)(?:\s+\[\d+\])?\s+(tests/[^:]+)::', line)
    if m: file_counts[m.group(2)][m.group(1).lower()] += 1
inv = {m: {'spec': m in specs, 'test_files': [], 'passed':0, 'failed':0, 'skipped':0} for m in mods}
meta_files, suite_totals = [], {'passed':0,'failed':0,'skipped':0}
for tf in sorted((root/'tests').rglob('test_*.py')):
    tree = ast.parse(tf.read_text())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith('sugarcode.modules.'):
            imported.add(node.module.split('.')[2])
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.startswith('sugarcode.modules.'): imported.add(a.name.split('.')[2])
    rel = str(tf.relative_to(root))
    c = file_counts.get(rel)
    if c is None: continue
    for k in ('passed','failed','skipped'): suite_totals[k] += c[k]
    if not imported:
        meta_files.append(tf.name); continue
    for m in imported:
        if m in inv:
            inv[m]['test_files'].append(tf.name)
            for k in ('passed','failed','skipped'): inv[m][k] += c[k]
zero = [m for m,v in inv.items() if not v['test_files']]
summary = {'modules_total': len(mods), 'modules_with_spec': sum(v['spec'] for v in inv.values()),
           'specs_total': len(specs),
           'modules_with_import_linked_tests': len(mods)-len(zero),
           'modules_without_import_linked_tests': zero,
           'meta_or_nonmodule_test_files': len(meta_files),
           'suite_totals': suite_totals, 'per_module': inv}
json.dump(summary, open('/home/sandbox/wave1/logs/selfaudit_inventory.json','w'), indent=1)
print(json.dumps({k:v for k,v in summary.items() if k!='per_module'}, indent=1))
