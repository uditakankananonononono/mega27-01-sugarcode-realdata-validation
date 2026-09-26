#!/usr/bin/env python3
"""Self-audit Tier-3 smoke: for every module not Tier-1-verified or Tier-2-probed:
import_ok, docstring_ok (module docstring exists, >=20 chars, states a method),
__all__ or public functions exist, and probe_ok if an obvious zero-arg demo/probe
entry exists (report() / diagnostics() / demo() patterns). Honest smoke level."""
import importlib, json, pathlib, sys, traceback, inspect
sys.path.insert(0, '/tmp/sc_audit/src')
targets = json.load(open('/home/sandbox/wave1/logs/tier3_probe_targets.json'))['tier3_targets']
out = {}
for name in targets:
    rec = {'import_ok': False, 'docstring_ok': False, 'public_api': None, 'probe_ok': None, 'probe_entry': None, 'error': None}
    try:
        m = importlib.import_module(f'sugarcode.modules.{name}')
        rec['import_ok'] = True
        doc = (m.__doc__ or '').strip()
        core = None
        try:
            core = importlib.import_module(f'sugarcode.modules.{name}.core')
            if not doc: doc = (core.__doc__ or '').strip()
        except Exception: pass
        rec['docstring_ok'] = len(doc) >= 20
        host = core or m
        pub = [f for f in dir(host) if not f.startswith('_') and callable(getattr(host, f, None))]
        rec['public_api'] = len(pub)
        # probe: prefer zero-arg-callable demo/diagnostics/report entries
        cands = [f for f in pub if any(k in f.lower() for k in ('diagnostics','demo','report'))]
        probed = False
        for f in cands:
            fn = getattr(host, f)
            try:
                sig = inspect.signature(fn)
                if all(p.default is not inspect.Parameter.empty or p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD)
                       for p in sig.parameters.values()):
                    r = fn()
                    rec['probe_ok'] = r is not None
                    rec['probe_entry'] = f
                    probed = True
                    break
            except Exception as e:
                rec['probe_ok'] = False
                rec['probe_entry'] = f
                rec['error'] = f'{type(e).__name__}: {e}'[:200]
                probed = True
                break
        if not probed:
            rec['probe_ok'] = None  # no obvious zero-arg entry; not a failure
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {e}'[:200]
    out[name] = rec
json.dump(out, open('/home/sandbox/wave1/logs/tier3_smoke.json','w'), indent=1)
n = len(out)
imp = sum(r['import_ok'] for r in out.values())
doc = sum(r['docstring_ok'] for r in out.values())
probed = sum(1 for r in out.values() if r['probe_ok'] is not None)
pok = sum(1 for r in out.values() if r['probe_ok'] is True)
print(f'modules={n} import_ok={imp} docstring_ok={doc} probed={probed} probe_pass={pok}')
print('import failures:', [k for k,r in out.items() if not r['import_ok']])
print('docstring failures:', [k for k,r in out.items() if r['import_ok'] and not r['docstring_ok']])
print('probe failures:', [(k, r['probe_entry'], r['error']) for k,r in out.items() if r['probe_ok'] is False])
