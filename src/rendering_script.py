import subprocess
from pathlib import Path

qmd_files = sorted(Path(".").rglob("*.qmd"))

failed = []

for qmd in qmd_files:
    print(f"Rendering {qmd}...")
    try:
        subprocess.run(
            ["quarto", "render", str(qmd)],
            check=True,
        )
    except subprocess.CalledProcessError:
        print(f"FAILED: {qmd}")
        failed.append(qmd)

print(f"\nFinished. Rendered {len(qmd_files) - len(failed)}/{len(qmd_files)} files.")

if failed:
    print("\nFailed files:")
    for f in failed:
        print(" ", f)
