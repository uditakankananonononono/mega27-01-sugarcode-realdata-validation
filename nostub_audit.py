"""Repo-wide no-stub audit: fails if any source file contains stub markers.
Adapted from the sugarcode-ai permanent audit convention."""
import pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parent
SRC = ROOT / "src"
MARKERS = re.compile(
    r"\b(TODO|FIXME|XXX|NotImplementedError|pass\s*#\s*stub|placeholder|"
    r"pseudo-?code|mock data|fake data|dummy data|lorem ipsum)\b",
    re.IGNORECASE,
)
ALLOW = {"test_nostub_audit.py", "nostub_audit.py"}

def offenders():
    bad = []
    for base in (SRC, ROOT / "tests", ROOT / "scripts"):
        if not base.exists():
            continue
        for p in base.rglob("*.py"):
            if p.name in ALLOW:
                continue
            text = p.read_text(encoding="utf-8", errors="replace")
            for i, line in enumerate(text.splitlines(), 1):
                if MARKERS.search(line):
                    bad.append(f"{p.relative_to(ROOT)}:{i}: {line.strip()[:120]}")
    return bad

if __name__ == "__main__":
    bad = offenders()
    if bad:
        print("\n".join(bad))
        sys.exit(1)
    print("no-stub audit: clean")
