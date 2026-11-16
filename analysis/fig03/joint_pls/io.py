"""Lightweight IO shared by computation and standalone document export."""
from pathlib import Path
import hashlib
import json

MODELS = ("M0_covariates", "M1_PFC", "M2_distributed432", "M3_nonPFC")
SCHEMES = ("joint", "shared")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1048576), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(table, path):
    table.to_csv(path, index=False, encoding="utf-8-sig", float_format="%.17g")


def own_output(output, kind):
    output = Path(output).resolve()
    marker = output / "joint_pls_owner.json"
    if output.exists() and any(output.iterdir()):
        if not marker.is_file() or json.loads(marker.read_text())["kind"] != kind:
            raise ValueError(f"Nonempty output is not owned by {kind}: {output}")
    output.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({"kind": kind}), encoding="utf-8")
