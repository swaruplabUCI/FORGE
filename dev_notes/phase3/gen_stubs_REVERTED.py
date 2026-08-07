#!/usr/bin/env python3
"""
Generate `stub:` blocks for FORGE Nextflow modules.

For each process lacking a stub, we:
  1. parse its `output:` section and collect every `path` pattern,
  2. copy the Groovy prelude from `script:` (the lines before the opening `\"\"\"`)
     so that interpolated output names like ${cell_type_safe} resolve identically,
  3. emit a stub that mkdir/touches each declared output,
  4. insert the stub immediately before `script:`.

Anything we cannot translate confidently is reported and left untouched rather
than guessed at.

Usage: gen_stubs.py <repo_root> [--apply]
"""
import re
import sys
from pathlib import Path

SECTION = re.compile(r'^\s*(input|output|script|stub|shell|exec|when):\s*$')
# `path "x"` / `path("x")` / `path 'x'` / `path('x')`
PATH_PAT = re.compile(r'\bpath\s*\(?\s*(["\'])(.+?)\1')
TRIPLE = re.compile(r'^\s*"""\s*$')


def find_processes(lines):
    """Yield (name, start_idx, end_idx_exclusive) for each top-level process."""
    out = []
    for i, ln in enumerate(lines):
        m = re.match(r'^process\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{', ln)
        if not m:
            continue
        depth = 0
        for j in range(i, len(lines)):
            depth += lines[j].count('{') - lines[j].count('}')
            if depth == 0 and j > i:
                out.append((m.group(1), i, j + 1))
                break
    return out


def section_bounds(lines, lo, hi, name):
    """Return (start, end) line range of a `name:` section within [lo, hi)."""
    start = None
    for i in range(lo, hi):
        m = SECTION.match(lines[i])
        if m and m.group(1) == name:
            start = i + 1
            break
    if start is None:
        return None
    for j in range(start, hi):
        if SECTION.match(lines[j]):
            return (start, j)
    return (start, hi)


def shell_for(pattern):
    """Translate one declared output pattern into a shell command, or None."""
    p = pattern.strip()
    if not p or p in ('*', '**', '.', './'):
        return None                      # too ambiguous to fabricate
    if p.startswith('/'):
        return None                      # absolute path: not a task output
    is_dir = p.endswith('/')
    p = p.rstrip('/')
    # Concretize globs: enhancer_*.h5ad -> enhancer_stub.h5ad
    p = p.replace('**', 'stub').replace('*', 'stub')
    p = re.sub(r'\?', 'x', p)
    if is_dir:
        return f'mkdir -p "{p}"'
    parent = str(Path(p).parent)
    cmds = []
    if parent not in ('.', ''):
        cmds.append(f'mkdir -p "{parent}"')
    cmds.append(f'touch "{p}"')
    return ' && '.join(cmds)


def build_stub(lines, lo, hi, indent='    '):
    """Return (stub_lines, problems) for the process spanning [lo, hi)."""
    ob = section_bounds(lines, lo, hi, 'output')
    if ob is None:
        return None, ['no output: section']
    body = ''.join(lines[ob[0]:ob[1]])
    patterns = [m.group(2) for m in PATH_PAT.finditer(body)]
    if not patterns:
        return None, ['output: section declares no path outputs']

    cmds, problems = [], []
    for pat in patterns:
        sh = shell_for(pat)
        if sh is None:
            problems.append(f'untranslatable output pattern: {pat!r}')
        elif sh not in cmds:
            cmds.append(sh)
    if not cmds:
        return None, problems or ['no translatable outputs']

    # Groovy prelude from script:, so ${...} in output names resolves the same way.
    sb = section_bounds(lines, lo, hi, 'script')
    prelude, depth = [], 0
    if sb is not None:
        for i in range(sb[0], sb[1]):
            if TRIPLE.match(lines[i]):
                break
            if lines[i].strip():
                prelude.append(lines[i].rstrip('\n'))
                depth += lines[i].count('{') - lines[i].count('}')
        else:
            problems.append('script: has no triple-quoted body (template/exec?)')
            return None, problems

    # A branching script (`script: if (...) { """..."""} else {...}`) leaves the
    # prelude with an open brace. Copying it verbatim would emit invalid Groovy,
    # so drop the prelude entirely — safe only if no output name needs a variable
    # the prelude defined.
    if depth != 0:
        needs_var = [p for p in patterns if '${' in p]
        if needs_var:
            problems.append(
                'branching script: AND interpolated output(s) '
                f'{needs_var} — needs a hand-written stub')
            return None, problems
        problems.append('branching script: — prelude dropped (no interpolated outputs)')
        prelude = []

    stub = [f'{indent}stub:\n']
    stub += [pl + '\n' for pl in prelude]
    stub.append(f'{indent}"""\n')
    for c in cmds:
        stub.append(f'{indent}{c}\n')
    stub.append(f'{indent}"""\n')
    stub.append('\n')
    return stub, problems


def main():
    root = Path(sys.argv[1])
    apply = '--apply' in sys.argv
    added = skipped = 0
    report = []

    for f in sorted(root.glob('modules/**/*.nf')):
        lines = f.read_text().splitlines(keepends=True)
        procs = find_processes(lines)
        edits = []
        for name, lo, hi in procs:
            if section_bounds(lines, lo, hi, 'stub') is not None:
                continue
            sb = section_bounds(lines, lo, hi, 'script')
            if sb is None:
                report.append(f'SKIP {f.relative_to(root)}::{name} — no script: section')
                skipped += 1
                continue
            stub, problems = build_stub(lines, lo, hi)
            if stub is None:
                report.append(f'SKIP {f.relative_to(root)}::{name} — {"; ".join(problems)}')
                skipped += 1
                continue
            for p in problems:
                report.append(f'WARN {f.relative_to(root)}::{name} — {p}')
            insert_at = sb[0] - 1              # the `script:` line itself
            edits.append((insert_at, stub))
            added += 1

        if edits and apply:
            for insert_at, stub in sorted(edits, reverse=True):
                lines[insert_at:insert_at] = stub
            f.write_text(''.join(lines))

    print('\n'.join(report) if report else '(no problems)')
    print(f'\n--- stubs {"written" if apply else "planned"}: {added} | skipped: {skipped} ---')


if __name__ == '__main__':
    main()
