"""Generate the <module>_reset_state routines that `qdyn --serve` calls between stages.

A persistent qdyn process runs every stage through the same setup path a fresh process
takes. Before each stage, every module variable must therefore hold the value it has
when a process starts: its declared initial value, or zero (module variables without
an initial value live in zero-filled static storage). Allocated memory is freed.

The routine of each module is generated from that module's declarations and written
between marker lines into its source:

    python3 tools/reset_state.py          # rewrite the routines in place
    python3 tools/reset_state.py --check  # exit 1 if any routine is out of date

test/q6/test_persistent_qdyn.py runs --check, so a variable added to a module without
regenerating fails the tests instead of leaking state from one stage into the next.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]  # src/q6
BEGIN = "! BEGIN GENERATED reset_state (tools/reset_state.py; do not edit by hand)"
END = "! END GENERATED reset_state"


@dataclass
class Module:
    file: str
    name: str
    #: Statements run first, e.g. closing what the module has open.
    before: list[str] = field(default_factory=list)
    #: Variables never reset (they must survive between stages).
    keep: set[str] = field(default_factory=set)
    #: Pointers that only point into memory owned elsewhere: nullified, never freed.
    borrowed: set[str] = field(default_factory=set)
    #: Statements run last.
    after: list[str] = field(default_factory=list)


MODULES = [
    Module("md.f90", "md", keep={"serve_mode", "serve_dir", "serve_has_state", "serve_x", "serve_v",
                                 "serve_has_theta", "serve_theta_corr", "serve_boxlength",
                                 "serve_boxcenter", "serve_topology", "serve_input",
                                 "serve_from_memory_count"}),
    Module("topo.f90", "topo"),
    Module("qatom.f90", "qatom"),
    Module("trj.f90", "trj", before=["call trj_close"]),
    Module("mask.f90", "atom_mask"),
    Module("prmfile.f90", "prmfile", before=["call prm_close", "reset_ok = prm_get_field(reset_field, skip=.true.)"],
           borrowed={"first_sec", "current_sec", "current_lp"}),
    Module("index.f90", "indexer", before=["call index_clear"], borrowed={"ndx"}),
    Module("mpiglob.f90", "mpiglob", keep={"nodeid", "numnodes"}),
    Module("lincs.f90", "lincs"),
    Module("lincs_mpi.f90", "lincs_mpi"),
    Module("settle.f90", "settle"),
]


@dataclass
class Variable:
    name: str
    decl: str  # attributes before ::, lower case
    init: str
    pp: tuple  # enclosing preprocessor conditions


def _logical_lines(text: str):
    """(physical line number, statement) with comments removed and continuations joined."""
    buffer, start = "", 0
    for number, raw in enumerate(text.splitlines(), 1):
        if raw.lstrip().startswith("#"):
            yield number, raw.strip()
            continue
        line, quote = raw, None
        for i, c in enumerate(raw):
            if quote:
                if c == quote:
                    quote = None
            elif c in "'\"":
                quote = c
            elif c == "!":
                line = raw[:i]
                break
        line = line.rstrip()
        if not buffer:
            start = number
        if line.endswith("&"):
            buffer += line[:-1] + " "
            continue
        buffer += line.lstrip().lstrip("&")
        yield start, buffer
        buffer = ""


def _split(names: str) -> list[str]:
    parts, depth, current = [], 0, ""
    for c in names:
        depth += c in "(["
        depth -= c in ")]"
        if c == "," and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += c
    return [p.strip() for p in [*parts, current] if p.strip()]


def module_variables(text: str, module: str) -> list[Variable]:
    variables, inside, depth, interface, pp = [], False, 0, 0, []
    for _, statement in _logical_lines(text):
        s = statement.strip()
        if not s:
            continue
        if s.startswith("#"):
            if re.match(r"#\s*if", s):
                pp.append(s)
            elif re.match(r"#\s*endif", s):
                pp.pop()
            elif re.match(r"#\s*else", s):
                pp[-1] = pp[-1] + " !else"
            continue
        if re.match(rf"module\s+{module}\s*$", s, re.I):
            inside = True
            continue
        if not inside:
            continue
        if re.match(r"contains\s*$", s, re.I):
            break
        if re.match(r"(abstract\s+)?interface\b", s, re.I):
            interface += 1
            continue
        if re.match(r"end\s*interface", s, re.I):
            interface -= 1
            continue
        if interface:
            continue
        if re.match(r"type\s*(,\s*[\w()]+\s*)*(::)?\s*\w+\s*$", s, re.I) and not re.match(r"type\s*\(", s, re.I):
            depth += 1
            continue
        if re.match(r"end\s*type", s, re.I):
            depth -= 1
            continue
        if depth or "::" not in s:
            continue
        decl, names = s.split("::", 1)
        lowered = decl.lower()
        if "parameter" in lowered or re.match(r"\s*(use|procedure|public|private)\b", lowered):
            continue
        for item in _split(names):
            match = re.match(r"(\w+)\s*(\([^=]*?\))?\s*(=>?\s*(.*))?$", item)
            if match is None:
                raise SystemExit(f"{module}: cannot read declaration {item!r}")
            init = (match.group(3) or "").strip()  # "= value" or "=> target"
            variables.append(Variable(match.group(1), lowered.strip(), init, tuple(pp)))
    return variables


def _reset(module: Module, variable: Variable, fresh: dict) -> list[str]:
    name, decl = variable.name, variable.decl
    if name in module.borrowed:
        return [f"nullify({name})"]
    if "pointer" in decl:
        return [f"deallocate({name}, stat=reset_status)", f"nullify({name})"]
    if "allocatable" in decl:
        return [f"if (allocated({name})) deallocate({name})"]
    if variable.init:
        if variable.init.startswith("=>"):
            raise SystemExit(f"{module.name}: {name} is initialised as a pointer; not supported")
        return [f"{name} = {variable.init.lstrip('=').strip()}"]
    kind = decl.split(",")[0].strip()
    if kind.startswith(("type", "class")):
        type_name = re.match(r"(type|class)\s*\(\s*(\w+)\s*\)", kind).group(2)
        fresh.setdefault(type_name, (f"fresh_{type_name.lower()}", variable.pp))
        return [f"{name} = {fresh[type_name][0]}"]
    if kind.startswith("logical"):
        return [f"{name} = .false."]
    if kind.startswith("character"):
        return [f"{name} = repeat(achar(0), len({name}))"]
    return [f"{name} = 0"]


def routine(module: Module, text: str) -> list[str]:
    variables = [v for v in module_variables(text, module.name) if v.name not in module.keep]
    fresh: dict = {}
    body, pp = [], ()
    for statement in module.before:
        body.append(f"  {statement}")
    for variable in variables:
        lines = _reset(module, variable, fresh)
        if variable.pp != pp:
            # close the conditions that end here, open the new ones
            common = 0
            while common < min(len(pp), len(variable.pp)) and pp[common] == variable.pp[common]:
                common += 1
            body += ["#endif"] * (len(pp) - common)
            body += [c.replace(" !else", "") for c in variable.pp[common:]]
            pp = variable.pp
        body += [f"  {line}" for line in lines]
    body += ["#endif"] * len(pp)
    for statement in module.after:
        body.append(f"  {statement}")
    preprocessed = any(line.startswith("#") for line in body)
    if preprocessed and module.file != "md.f90":
        raise SystemExit(f"{module.file} is not preprocessed: it cannot have #if around its variables")
    declarations = ["  integer :: reset_status"]
    if module.name == "prmfile":  # prm_get_field keeps a pointer into the input between calls
        declarations += ["  logical :: reset_ok", "  character(len=1) :: reset_field"]
    for type_name, (local, conditions) in fresh.items():
        declarations += [c.replace(" !else", "") for c in conditions]
        declarations.append(f"  type({type_name}), save :: {local}")
        declarations += ["#endif"] * len(conditions)
    return [
        BEGIN,
        f"subroutine {module.name}_reset_state",
        "!! Every variable of this module back to its value at program start (qdyn --serve).",
        *declarations,
        "  reset_status = 0",
        *body,
        f"end subroutine {module.name}_reset_state",
        END,
    ]


PUBLIC = "  public :: {name}_reset_state ! generated by tools/reset_state.py"


def place(module: Module) -> tuple[str, str]:
    """(current text, text with an up-to-date routine) of a module's source."""
    path = HERE / module.file
    text = path.read_text()
    lines = text.splitlines()
    if BEGIN in lines:
        start, end = lines.index(BEGIN), lines.index(END)
        rest = lines[:start] + lines[end + 1:]
    else:
        rest = lines
        contains = [i for i, line in enumerate(rest) if re.match(r"\s*contains\s*$", line, re.I)]
        if contains:
            start = contains[0] + 1
        else:  # a module without procedures gets a contains section
            end_module = next(i for i, line in enumerate(rest) if re.match(r"\s*end\s+module", line, re.I))
            rest = rest[:end_module] + ["contains"] + rest[end_module:]
            start = end_module + 1
    new = rest[:start] + routine(module, "\n".join(rest)) + rest[start:]
    public = PUBLIC.format(name=module.name)
    if public not in new:
        # modules whose default is private export it from right after `private`
        for i, line in enumerate(new):
            if re.match(r"\s*private\s*$", line, re.I):
                new.insert(i + 1, public)
                break
    return text, "\n".join(new) + ("\n" if text.endswith("\n") else "")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="Exit 1 if a routine is missing or out of date.")
    args = parser.parse_args(argv)
    stale = []
    for module in MODULES:
        current, updated = place(module)
        if current != updated:
            stale.append(module.file)
            if not args.check:
                (HERE / module.file).write_text(updated)
    if args.check and stale:
        print("Out of date: " + ", ".join(stale) + ". Run python3 src/q6/tools/reset_state.py", file=sys.stderr)
        return 1
    if not args.check:
        print("Updated: " + (", ".join(stale) or "nothing"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
