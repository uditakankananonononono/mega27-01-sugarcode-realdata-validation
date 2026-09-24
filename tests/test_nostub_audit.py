import subprocess, sys, pathlib
def test_no_stubs_repo_wide():
    root = pathlib.Path(__file__).resolve().parent.parent
    r = subprocess.run([sys.executable, str(root / "nostub_audit.py")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout
