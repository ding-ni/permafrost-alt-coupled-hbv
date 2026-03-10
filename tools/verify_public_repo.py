# -*- coding: utf-8 -*-
import json
from pathlib import Path

try:
    from jsonschema import validate as jsonschema_validate
except ImportError:  # pragma: no cover - fallback for bare environments
    jsonschema_validate = None


REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = REPO_ROOT / "config" / "basin.schema.json"
TEMPLATE_PATH = REPO_ROOT / "config" / "basin.template.json"
SMOKE_CONFIG_PATH = REPO_ROOT / "examples" / "synthetic_basin" / "smoke_config.json"
FORBIDDEN_TERMS = [
    "F:\\Hapi",
    "Tuotuohe",
    "tuotuohe",
    "沱沱河",
    "S32",
    "S40",
    "post2018",
]
TEXT_SUFFIXES = {".py", ".md", ".json", ".yml", ".yaml", ".cff", ".txt", ".gitignore"}


def iter_text_files():
    for path in REPO_ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.suffix.lower() in TEXT_SUFFIXES or path.name == ".gitignore":
            yield path


def validate_configs():
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    for path in [TEMPLATE_PATH, SMOKE_CONFIG_PATH]:
        if not path.exists():
            raise FileNotFoundError(f"Missing config: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if jsonschema_validate is not None:
            jsonschema_validate(instance=payload, schema=schema)
            continue

        required = schema.get("required", [])
        missing = [key for key in required if key not in payload]
        if missing:
            raise ValueError(f"{path.name} is missing required keys: {missing}")


def scan_forbidden_terms():
    hits = []
    for path in iter_text_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        for term in FORBIDDEN_TERMS:
            if term in text:
                hits.append((path, term))
    return hits


def main():
    validate_configs()
    hits = scan_forbidden_terms()
    if hits:
        for path, term in hits:
            print(f"FORBIDDEN\t{term}\t{path.relative_to(REPO_ROOT)}")
        raise SystemExit(1)
    print("Repository verification passed.")


if __name__ == "__main__":
    main()
