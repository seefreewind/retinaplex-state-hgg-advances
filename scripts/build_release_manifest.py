from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "metadata" / "release_file_manifest.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    rows = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == OUTPUT or path.name.startswith("._") or path.name == ".DS_Store":
            continue
        rows.append((path.relative_to(ROOT).as_posix(), path.stat().st_size, sha256(path)))
    OUTPUT.write_text(
        "path\tsize_bytes\tsha256\n"
        + "".join(f"{path}\t{size}\t{digest}\n" for path, size, digest in rows),
        encoding="utf-8",
    )
    print(f"Wrote {OUTPUT} with {len(rows)} entries")


if __name__ == "__main__":
    main()
