"""Which distance restraints between the two ligand copies to keep.

QligFEP ties every mapped heavy atom of one ligand copy to its partner in the other. For atoms
past a rotatable bond (the methyl of a methoxy group, say) this couples the two copies'
rotations: the ghost copy, which has only bonded terms, can stay frozen in its starting
orientation and hold the real ligand there through the restraint. The end states then depend on
the starting pose and the error does not cancel between the legs.

release_rotors() keeps a restraint only when both of its atoms are rigid with respect to a ring:
in a ring or bonded to a ring atom.
"""
from pathlib import Path


def topology_bonds(path):
    """Bonds (0-based atom pairs) of a Q topology file."""
    lines = Path(path).read_text().splitlines()
    for at, line in enumerate(lines):
        if "No. of bonds" in line:
            count = int(line.split()[0])
            numbers = []
            for row in lines[at + 1:]:
                numbers.extend(int(x) for x in row.split())
                if len(numbers) >= 3 * count:
                    break
            return [(numbers[3 * k] - 1, numbers[3 * k + 1] - 1) for k in range(count)]
    raise ValueError(f"{path}: no bond list")


def ring_rigid_atoms(bonds, atoms):
    """Of `atoms` (0-based), those in a ring or bonded to a ring atom, with the bonds among them."""
    atoms = set(atoms)
    neighbours = {a: set() for a in atoms}
    for i, j in bonds:
        if i in atoms and j in atoms:
            neighbours[i].add(j)
            neighbours[j].add(i)

    def in_ring(i, j):     # is j reachable from i without the bond i-j?
        seen, todo = {i}, [i]
        while todo:
            a = todo.pop()
            for b in neighbours[a]:
                if (a, b) in ((i, j), (j, i)) or b in seen:
                    continue
                if b == j:
                    return True
                seen.add(b)
                todo.append(b)
        return False

    ring = {a for a in atoms for b in neighbours[a] if in_ring(a, b)}
    return ring | {a for a in atoms if neighbours[a] & ring}


def release_rotors(restraints, bonds, q_atoms):
    """(kept, dropped): the distance restraints whose two atoms are both ring-rigid, and the rest."""
    rigid = ring_rigid_atoms(bonds, q_atoms)
    kept = [r for r in restraints if r.ai - 1 in rigid and r.aj - 1 in rigid]
    dropped = [r for r in restraints if r not in kept]
    return kept, dropped


def topology_masses(path):
    """Mass of every atom (by its type code) from a Q topology file."""
    lines = Path(path).read_text().splitlines()

    def block(at, count, kind):
        values = []
        for row in lines[at + 1:]:
            values.extend(kind(x) for x in row.split())
            if len(values) >= count:
                return values[:count]

    iac, masses, n_types = None, None, None
    for at, line in enumerate(lines):
        if "No. of integer atom codes" in line:
            iac = block(at, int(line.split()[0]), int)
        elif "No. of atom types" in line and line.rstrip().endswith("types"):
            n_types = int(line.split()[0])
        elif line.strip() == "Masses:" and n_types:
            masses = block(at, n_types, float)
    return [masses[code - 1] for code in iac]


def ring_systems(bonds, atoms):
    """Ring atoms of `atoms`, grouped into ring systems (rings that share or bond ring atoms)."""
    rigid = ring_rigid_atoms(bonds, atoms)
    neighbours = {a: set() for a in atoms}
    for i, j in bonds:
        if i in neighbours and j in neighbours:
            neighbours[i].add(j)
            neighbours[j].add(i)
    ring = {a for a in rigid if len(neighbours[a] & rigid) >= 2 and
            any(_in_cycle(neighbours, a, b) for b in neighbours[a])}
    systems, seen = [], set()
    for start in ring:
        if start in seen:
            continue
        group, todo = {start}, [start]
        while todo:
            a = todo.pop()
            for b in neighbours[a] & ring:
                if b not in group and _in_cycle(neighbours, a, b):
                    group.add(b)
                    todo.append(b)
        seen |= group
        systems.append(group)
    return systems, neighbours


def _in_cycle(neighbours, i, j):
    seen, todo = {i}, [i]
    while todo:
        a = todo.pop()
        for b in neighbours[a]:
            if (a, b) in ((i, j), (j, i)) or b in seen:
                continue
            if b == j:
                return True
            seen.add(b)
            todo.append(b)
    return False


def restrain_core(restraints, bonds, q_atoms, masses):
    """(kept, dropped): keep only restraints inside the rigid core around the perturbation.

    The core of each ligand copy is the ring system that carries an atom whose element changes
    between the copies (a restrained pair of different mass, here O -> S), plus the atoms bonded to
    that ring system. Everything past a rotatable bond from it (other rings included) is released,
    so the two copies' flexible parts no longer move as one.
    """
    changed = {a for r in restraints if abs(masses[r.ai - 1] - masses[r.aj - 1]) > 0.5
               for a in (r.ai - 1, r.aj - 1)}
    systems, neighbours = ring_systems(bonds, q_atoms)
    core = set()
    for group in systems:
        touching = group | {b for a in group for b in neighbours[a]}
        if touching & changed:
            core |= touching
    if not core:     # no element change: fall back to the ring-rigid rule
        return release_rotors(restraints, bonds, q_atoms)
    kept = [r for r in restraints if r.ai - 1 in core and r.aj - 1 in core]
    return kept, [r for r in restraints if r not in kept]


def release_past_joints(restraints, bonds, q_atoms, masses):
    """(kept, dropped): release the restraints on atoms past an sp3 joint of the ligand.

    A joint is a divalent O or S or a carbon with four neighbours (an sp3 centre: the bonds around
    it rotate freely). Starting from the ring system that carries the changed element (or from
    every ring when no element changes), the kept region grows through bonds but does not pass a
    joint: a joint bonded to the region is kept, the atoms behind it are released. Planar links
    such as an aryl N-H stay kept, so rigidly joined rings (the hinge-binding purine) stay tied.
    """
    atoms = set(q_atoms)
    neighbours = {a: set() for a in atoms}
    for i, j in bonds:
        if i in atoms and j in atoms:
            neighbours[i].add(j)
            neighbours[j].add(i)

    def joint(a):
        n, m = len(neighbours[a]), masses[a]
        return (abs(m - 12.0) < 0.5 and n == 4) or (abs(m - 16.0) < 0.5 and n == 2) or (abs(m - 32.06) < 0.6 and n == 2)

    seeds = perturbed_rings(restraints, neighbours, bonds, q_atoms, masses)
    kept_region, todo = set(seeds), list(seeds)
    while todo:
        a = todo.pop()
        for b in neighbours[a]:
            if b in kept_region:
                continue
            kept_region.add(b)
            if not joint(b) or b in seeds:      # a joint is kept but not passed
                todo.append(b)
    kept = [r for r in restraints if r.ai - 1 in kept_region and r.aj - 1 in kept_region]
    return kept, [r for r in restraints if r not in kept]


def perturbed_rings(restraints, neighbours, bonds, q_atoms, masses):
    """Atoms of the ring system(s) nearest to the perturbation.

    The perturbation is where the two ligands differ: atoms whose restrained partner has another
    element, and mapped heavy atoms bonded to a heavy atom that has no partner (a group present in
    one ligand only). The ring systems at the smallest bond distance from those atoms are returned;
    with no perturbation found, all ring systems.
    """
    mapped = {a for r in restraints for a in (r.ai - 1, r.aj - 1)}
    heavy = {a for a in neighbours if masses[a] > 2.0}
    perturbed = {a for r in restraints if abs(masses[r.ai - 1] - masses[r.aj - 1]) > 0.5
                 for a in (r.ai - 1, r.aj - 1)}
    perturbed |= {a for a in mapped for b in neighbours[a] if b in heavy and b not in mapped}
    systems, _ = ring_systems(bonds, q_atoms)
    if not systems:
        return set()
    if not perturbed:
        return set().union(*systems)
    distance, todo = {a: 0 for a in perturbed}, list(perturbed)
    while todo:                                    # breadth-first distances over the bonds
        nxt = []
        for a in todo:
            for b in neighbours[a]:
                if b not in distance:
                    distance[b] = distance[a] + 1
                    nxt.append(b)
        todo = nxt
    near = [min(distance.get(a, 10 ** 6) for a in g) for g in systems]
    best = min(near)
    seeds = set().union(*[g for g, d in zip(systems, near) if d == best])
    # the same rings in the other ligand copy: a restraint needs both of its atoms kept
    partner = {}
    for r in restraints:
        partner[r.ai - 1], partner[r.aj - 1] = r.aj - 1, r.ai - 1
    across = {partner[a] for a in seeds if a in partner}
    return seeds | set().union(*[g for g in systems if g & across]) if across else seeds
