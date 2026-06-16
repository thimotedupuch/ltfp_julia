import shutil
from pathlib import Path

src = Path(".")
dst = src / "docs"

if dst.exists():
    shutil.rmtree(dst)

for html in src.rglob("*.html"):
    if html.parts[0] == "docs":
        continue

    out = dst / html.relative_to(src)
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(html, out)
