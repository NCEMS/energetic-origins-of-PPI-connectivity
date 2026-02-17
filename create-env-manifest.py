#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path
import csv

SNAKEFILE_NAME = "Snakefile"

RULE_RE = re.compile(r'^\s*rule\s+([A-Za-z0-9_]+)\s*:\s*$')

# Matches conda: <something> on the same line (quoted or bare token)
CONDA_INLINE_RE = re.compile(r'^\s*conda\s*:\s*(?:"([^"]+)"|\'([^\']+)\'|([^\s#]+))\s*$')

# Matches a line that starts a conda block: "conda:" and nothing else
CONDA_BLOCK_START_RE = re.compile(r'^\s*conda\s*:\s*$')

# Within any python expression block, try to find a quoted .yml/.yaml substring
YAML_QUOTED_RE = re.compile(r'["\']([^"\']+\.ya?ml)["\']')

# Snakemake section keys that typically end a conda block.
# (This doesn't need to be perfect; it just prevents swallowing the whole rule.)
SECTION_KEY_RE = re.compile(
    r'^\s*(input|output|params|threads|resources|log|benchmark|message|priority|'
    r'group|shell|run|script|notebook|wrapper|container|conda|envmodules)\s*:\s*'
)

def indent_level(line: str) -> int:
    return len(line) - len(line.lstrip(" "))

def extract_env_from_inline_match(m: re.Match) -> str:
    return next(g for g in m.groups() if g is not None)

def extract_yaml_from_block(block_lines: list[str]) -> str | None:
    """
    Given lines inside a conda: block (already indented), return a YAML path if we can find one.
    We look for any quoted substring that ends in .yml/.yaml.
    """
    blob = "\n".join(block_lines)
    m = YAML_QUOTED_RE.search(blob)
    if m:
        return m.group(1)
    # Fallback: sometimes people put a bare token like env/foo.yml without quotes
    m2 = re.search(r'(^|\s)([^"\']+\.(?:yml|yaml))(\s|$)', blob)
    if m2:
        return m2.group(2).strip()
    return None

def main(repo_root: Path) -> None:
    rows: list[dict[str, str]] = []

    for snakefile in repo_root.rglob(SNAKEFILE_NAME):
        if any(part.startswith(".") for part in snakefile.parts):
            continue

        rel_snakefile = str(snakefile.relative_to(repo_root))
        lines = snakefile.read_text(errors="ignore").splitlines()

        current_rule: str | None = None
        i = 0
        while i < len(lines):
            line = lines[i]

            m_rule = RULE_RE.match(line)
            if m_rule:
                current_rule = m_rule.group(1)
                i += 1
                continue

            # inline conda
            m_inline = CONDA_INLINE_RE.match(line)
            if m_inline:
                env = extract_env_from_inline_match(m_inline)
                rows.append({"snakefile": rel_snakefile, "rule": current_rule or "", "env": env})
                i += 1
                continue

            # block conda
            if CONDA_BLOCK_START_RE.match(line):
                start_indent = indent_level(line)
                block_lines: list[str] = []
                i += 1
                while i < len(lines):
                    nxt = lines[i]
                    # stop if indentation returns to start_indent or less AND it's a new section key or rule
                    if indent_level(nxt) <= start_indent and (SECTION_KEY_RE.match(nxt) or RULE_RE.match(nxt)):
                        break
                    # also stop at an unindented blank line that often separates rules
                    if indent_level(nxt) == 0 and nxt.strip() == "":
                        break
                    block_lines.append(nxt)
                    i += 1

                env = extract_yaml_from_block(block_lines)
                if env:
                    rows.append({"snakefile": rel_snakefile, "rule": current_rule or "", "env": env})
                continue

            i += 1

    # de-duplicate
    uniq = {(r["snakefile"], r["rule"], r["env"]) for r in rows}
    rows = [{"snakefile": a, "rule": b, "env": c} for (a, b, c) in sorted(uniq)]

    out = repo_root / "env-manifest.tsv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["snakefile", "rule", "env"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {out} with {len(rows)} entries.")
    if len(rows) == 0:
        print("No conda envs detected. If you use 'conda:' via variables with no quoted .yml/.yaml anywhere, "
              "we'll need a different strategy (see notes).")

if __name__ == "__main__":
    main(Path(".").resolve())
