#!/usr/bin/env python3
"""
cremer_pople.py

Cremer-Pople puckering parameters and endocyclic torsion angles of a
six-membered ring, from Cartesian coordinates.

Input   one or more coordinate files in XYZ format (Angstrom).
Output  Q, theta, phi (and q2, q3), the six endocyclic torsion angles and
        the nearest canonical conformations. For a transition structure in
        which an oxygen nucleophile approaches the ring, the forming O...C
        and breaking C...N distances and the dihedral angles between the
        forming bond and the vicinal bonds are reported as well.

Usage
    python cremer_pople.py structure.xyz --ring 1,2,3,4,5,6
    python cremer_pople.py *.xyz --csv puckering.csv

Ring numbering
    --ring  six 1-based atom indices, in ring order, atom 1 first. This
            works for any six-membered ring, carbocyclic or heterocyclic.

    Without --ring the script numbers the ring itself, for one class of
    molecule only: a six-membered carbocycle carrying a side chain that ends
    in an amide, in which two adjacent ring carbons are bonded to (or in
    close contact with) nitrogen or oxygen. C1 is the ring carbon bearing
    the side chain; C3 and C4 are the two heteroatom-bearing carbons, C3
    being two bonds from C1. Files that do not fit are skipped with a reason.

Conventions
    Puckering coordinates follow Cremer and Pople, J. Am. Chem. Soc. 1975,
    97, 1354-1358 (equation numbers below refer to that paper). Canonical
    conformations follow Boeyens, J. Cryst. Mol. Struct. 1978, 8, 317-320.
    Torsion angles carry the Klyne-Prelog sign.

Cautions
    theta and phi depend on the numbering: moving the origin by one atom, or
    reversing the direction, converts theta to 180 - theta and shifts phi.
    The mirror image of a structure, numbered in the same way, also gives
    180 - theta. Compare values only between structures numbered alike and
    belonging to the same enantiomeric series. phi is poorly defined when
    theta is close to 0 or 180 degrees.

Requires numpy.
"""

import argparse
import csv
import math
import os
import sys
from collections import deque

import numpy as np

# Covalent radii (Angstrom) used to perceive bonds.
RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05,
         "Cl": 1.02, "Br": 1.20, "I": 1.39}
BOND_FACTOR = 1.25      # bonded if d < factor * (r_i + r_j)
CONTACT_FACTOR = 1.40   # looser cutoff: keeps the breaking C...N bond of a TS
NUCLEOPHILE_MAX = 2.80  # longest forming O...C distance accepted (Angstrom)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def read_xyz(path):
    """Return (symbols, coordinates) from an XYZ file."""
    with open(path) as handle:
        lines = handle.read().splitlines()
    natoms = int(lines[0].split()[0])
    symbols, coords = [], []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        symbols.append(parts[0].capitalize())
        coords.append([float(v) for v in parts[1:4]])
    if len(symbols) != natoms:
        raise ValueError(f"{path}: expected {natoms} atoms, found {len(symbols)}")
    unknown = sorted(set(symbols) - set(RADII))
    if unknown:
        raise ValueError(f"{path}: no covalent radius for {', '.join(unknown)}")
    return symbols, np.array(coords)


def distance(xyz, i, j):
    return float(np.linalg.norm(xyz[i] - xyz[j]))


def dihedral(xyz, a, b, c, d):
    """Signed dihedral angle a-b-c-d in degrees (Klyne-Prelog convention)."""
    b1 = xyz[b] - xyz[a]
    b2 = xyz[c] - xyz[b]
    b3 = xyz[d] - xyz[c]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    sine = float(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)))
    return math.degrees(math.atan2(sine, float(np.dot(n1, n2))))


def neighbours(symbols, xyz, factor=BOND_FACTOR):
    """Adjacency list from interatomic distances."""
    n = len(symbols)
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if distance(xyz, i, j) < factor * (RADII[symbols[i]] + RADII[symbols[j]]):
                adj[i].add(j)
                adj[j].add(i)
    return adj


def path_length(adj, start, goal, blocked=()):
    """Number of bonds on the shortest path from start to goal."""
    queue, seen = deque([(start, 0)]), {start}
    while queue:
        atom, depth = queue.popleft()
        if atom == goal:
            return depth
        for nxt in adj[atom]:
            if nxt not in seen and nxt not in blocked:
                seen.add(nxt)
                queue.append((nxt, depth + 1))
    return None


def fragment(adj, start):
    """All atoms connected to start."""
    queue, seen = deque([start]), {start}
    while queue:
        for nxt in adj[queue.popleft()]:
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return seen


# --------------------------------------------------------------------------
# Cremer-Pople puckering coordinates for a six-membered ring
# --------------------------------------------------------------------------

def cremer_pople(xyz, ring):
    """
    Puckering coordinates of a six-membered ring.

    ring: six atom indices in order (atom 1 first).
    Returns Q, theta, phi, q2, q3 and the displacements z_j from the mean
    plane; angles in degrees, lengths in Angstrom.
    """
    n = 6
    pos = xyz[ring] - xyz[ring].mean(axis=0)                  # eq 4
    angle = [2.0 * math.pi * j / n for j in range(n)]
    r1 = sum(pos[j] * math.sin(angle[j]) for j in range(n))   # eq 8
    r2 = sum(pos[j] * math.cos(angle[j]) for j in range(n))   # eq 9
    normal = np.cross(r1, r2)
    normal /= np.linalg.norm(normal)                          # eq 10
    z = pos @ normal                                          # eq 11

    q2cos = math.sqrt(2.0 / n) * sum(z[j] * math.cos(2 * angle[j]) for j in range(n))   # eq 12
    q2sin = -math.sqrt(2.0 / n) * sum(z[j] * math.sin(2 * angle[j]) for j in range(n))  # eq 13
    q3 = math.sqrt(1.0 / n) * sum(z[j] * (-1) ** j for j in range(n))                   # eq 14

    q2 = math.hypot(q2cos, q2sin)
    total = math.hypot(q2, q3)                                # eq 17
    theta = math.degrees(math.atan2(q2, q3))                  # eq 22, 0-180
    phi = math.degrees(math.atan2(q2sin, q2cos)) % 360.0
    return total, theta, phi, q2, q3, z


def ring_torsions(xyz, ring):
    """Endocyclic torsion angles about C1-C2, C2-C3, ... C6-C1."""
    return [dihedral(xyz, ring[k - 1], ring[k], ring[(k + 1) % 6], ring[(k + 2) % 6])
            for k in range(6)]


def canonical_forms():
    """Canonical conformations on the Cremer-Pople sphere (Boeyens 1978)."""
    half_chair = math.degrees(math.atan(math.sqrt(1.5)))     # 50.77
    envelope = math.degrees(math.atan(math.sqrt(2.0)))       # 54.74
    screw_boat = math.degrees(math.atan(1.0 + math.sqrt(2.0)))  # 67.50
    forms = [("chair", 0.0, 0.0), ("chair", 180.0, 0.0)]
    for k in range(6):
        forms.append(("boat", 90.0, 60.0 * k))
        forms.append(("twist-boat", 90.0, 30.0 + 60.0 * k))
        for theta in (envelope, 180.0 - envelope):
            forms.append(("envelope", theta, 60.0 * k))
        for theta in (half_chair, 180.0 - half_chair):
            forms.append(("half-chair", theta, 30.0 + 60.0 * k))
        for theta in (screw_boat, 180.0 - screw_boat):
            forms.append(("screw-boat", theta, 30.0 + 60.0 * k))
    return forms


def nearest_forms(theta, phi, count=3):
    """Nearest canonical forms by angular separation on the sphere."""
    t, p = math.radians(theta), math.radians(phi)
    found = []
    for name, ct, cp in canonical_forms():
        a, b = math.radians(ct), math.radians(cp)
        cosine = math.cos(t) * math.cos(a) + math.sin(t) * math.sin(a) * math.cos(p - b)
        found.append((math.degrees(math.acos(max(-1.0, min(1.0, cosine)))), name, ct, cp))
    found.sort()
    return found[:count]


# --------------------------------------------------------------------------
# Ring perception and numbering
# --------------------------------------------------------------------------

def carbocycles(symbols, adj):
    """All six-membered rings composed of carbon only (each listed once)."""
    rings, seen = [], set()

    def walk(path):
        if len(path) == 6:
            if path[0] in adj[path[-1]] and frozenset(path) not in seen:
                seen.add(frozenset(path))
                rings.append(list(path))
            return
        for nxt in adj[path[-1]]:
            if symbols[nxt] == "C" and nxt not in path:
                walk(path + [nxt])

    for start, symbol in enumerate(symbols):
        if symbol == "C":
            walk([start])
    return rings


def number_ring(symbols, xyz, adj):
    """
    Return the carbocycle as [C1 ... C6] together with the amide nitrogen.

    C1 bears the side chain leading to the amide carbonyl carbon; C3 and C4
    are the ring carbons in contact with a heteroatom, C3 being two bonds
    from C1.
    """
    nitrogens = [i for i, s in enumerate(symbols) if s == "N"]
    if len(nitrogens) != 1:
        raise ValueError("expected exactly one nitrogen; use --ring")
    nitrogen = nitrogens[0]

    carbonyls = [c for c in adj[nitrogen] if symbols[c] == "C"
                 and any(symbols[o] == "O" and distance(xyz, c, o) < 1.30 for o in adj[c])]
    if len(carbonyls) != 1:
        raise ValueError("amide carbonyl carbon not found; use --ring")
    carbonyl = carbonyls[0]

    rings = carbocycles(symbols, adj)
    if len(rings) != 1:
        raise ValueError(f"found {len(rings)} six-membered carbocycles; use --ring")
    ring = rings[0]

    # C1: shortest exocyclic path to the carbonyl carbon, avoiding N and the ring.
    def side_chain(atom):
        blocked = (set(ring) - {atom}) | {nitrogen}
        length = path_length(adj, atom, carbonyl, blocked)
        return length if length is not None else 99

    c1 = min(ring, key=side_chain)
    if side_chain(c1) == 99:
        raise ValueError("no ring carbon is linked to the carbonyl carbon; use --ring")
    start = ring.index(c1)
    ring = ring[start:] + ring[:start]

    # C3 and C4: ring carbons with a nitrogen or oxygen within the contact cutoff.
    def hetero(atom):
        return any(symbols[j] in ("N", "O") and j != atom and
                   distance(xyz, atom, j) < CONTACT_FACTOR * (RADII["C"] + RADII[symbols[j]])
                   for j in range(len(symbols)))

    if not (hetero(ring[2]) and hetero(ring[3])):
        ring = [ring[0]] + ring[1:][::-1]
    if not (hetero(ring[2]) and hetero(ring[3])):
        raise ValueError("C3/C4 could not be assigned; use --ring")
    return ring, nitrogen


# --------------------------------------------------------------------------
# Transition-structure details
# --------------------------------------------------------------------------

def nucleophile_details(symbols, xyz, adj, ring, nitrogen):
    """
    Forming and breaking bonds of a transition structure with an O nucleophile.

    The nucleophile is an oxygen of a separate fragment lying within
    NUCLEOPHILE_MAX of C3 or C4. Returns None when there is none.
    """
    cation = fragment(adj, nitrogen)
    candidates = [(distance(xyz, carbon, o), o, carbon)
                  for o, s in enumerate(symbols) if s == "O" and o not in cation
                  for carbon in (ring[2], ring[3])]
    if not candidates:
        return None
    forming, oxygen, attacked = min(candidates)
    if forming > NUCLEOPHILE_MAX:
        return None

    position = ring.index(attacked)                 # 2 (C3) or 3 (C4)
    partner = ring[3] if position == 2 else ring[2]     # the other heteroatom-bearing carbon
    vicinal_c = ring[1] if position == 2 else ring[4]   # C2 or C5
    substituents = sorted(adj[vicinal_c] - {attacked})
    labels = {atom: f"C{k + 1}" for k, atom in enumerate(ring)}
    vicinal = [(labels.get(s, symbols[s]), dihedral(xyz, oxygen, attacked, vicinal_c, s))
               for s in substituents]
    return {
        "attacked": labels[attacked],
        "site": "C3" if position == 2 else "C4",
        "forming": forming,
        "breaking": distance(xyz, attacked, nitrogen),
        "retained": distance(xyz, partner, nitrogen),
        "o_c_c_n": dihedral(xyz, oxygen, attacked, partner, nitrogen),
        "about": f"{labels[attacked]}-{labels[vicinal_c]}",
        "vicinal": vicinal,
        "closest": min(abs(angle) for _, angle in vicinal),
    }


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def analyse(path, manual_ring=None):
    symbols, xyz = read_xyz(path)
    adj = neighbours(symbols, xyz)
    if manual_ring:
        ring, nitrogen = manual_ring, None
        if len(ring) != 6 or max(ring) >= len(symbols) or min(ring) < 0:
            raise ValueError("--ring needs six valid atom indices")
    else:
        ring, nitrogen = number_ring(symbols, xyz, adj)

    total, theta, phi, q2, q3, z = cremer_pople(xyz, ring)
    result = {
        "file": os.path.basename(path),
        "ring": [i + 1 for i in ring],
        "Q": total, "theta": theta, "phi": phi, "q2": q2, "q3": q3,
        "z": list(z),
        "torsions": ring_torsions(xyz, ring),
        "nearest": nearest_forms(theta, phi),
        "ts": None,
    }
    if nitrogen is not None:
        result["ts"] = nucleophile_details(symbols, xyz, adj, ring, nitrogen)
    return result


def report(result):
    print(f"\n{result['file']}")
    print("  ring atoms C1-C6 (file indices): " + ", ".join(str(i) for i in result["ring"]))
    print(f"  Q = {result['Q']:.3f} A   theta = {result['theta']:.1f}   phi = {result['phi']:.1f}"
          f"   (q2 = {result['q2']:.3f}, q3 = {result['q3']:+.3f})")
    names = ("C1-C2", "C2-C3", "C3-C4", "C4-C5", "C5-C6", "C6-C1")
    print("  endocyclic torsions: " +
          "  ".join(f"{n} {t:+.1f}" for n, t in zip(names, result["torsions"])))
    print("  nearest canonical forms: " +
          "; ".join(f"{name} ({ct:.1f}, {cp:.0f}) at {sep:.0f} deg"
                    for sep, name, ct, cp in result["nearest"]))
    ts = result["ts"]
    if ts:
        print(f"  nucleophile attacks {ts['attacked']}: "
              f"O...C {ts['forming']:.3f} A, C...N {ts['breaking']:.3f} A, "
              f"other C-N {ts['retained']:.3f} A, O-C-C-N {ts['o_c_c_n']:+.1f}")
        print(f"  forming bond viewed along {ts['about']}: " +
              ", ".join(f"{label} {angle:+.1f}" for label, angle in ts["vicinal"]) +
              f"   (closest {ts['closest']:.1f})")


def write_csv(results, path):
    header = ["file", "ring_atoms", "Q_A", "theta_deg", "phi_deg", "q2_A", "q3_A",
              "t_C1C2", "t_C2C3", "t_C3C4", "t_C4C5", "t_C5C6", "t_C6C1",
              "nearest_form", "separation_deg", "attack_site", "O_C_A", "C_N_A",
              "closest_vicinal_deg"]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for r in results:
            ts = r["ts"] or {}
            writer.writerow(
                [r["file"], " ".join(str(i) for i in r["ring"]),
                 f"{r['Q']:.3f}", f"{r['theta']:.1f}", f"{r['phi']:.1f}",
                 f"{r['q2']:.3f}", f"{r['q3']:.3f}"]
                + [f"{t:.1f}" for t in r["torsions"]]
                + [r["nearest"][0][1], f"{r['nearest'][0][0]:.1f}",
                   ts.get("site", ""),
                   f"{ts['forming']:.3f}" if ts else "",
                   f"{ts['breaking']:.3f}" if ts else "",
                   f"{ts['closest']:.1f}" if ts else ""])


def main():
    parser = argparse.ArgumentParser(
        description="Cremer-Pople puckering parameters of a six-membered carbocycle.")
    parser.add_argument("files", nargs="+", help="XYZ coordinate files")
    parser.add_argument("--ring", help="six 1-based atom indices, C1...C6, "
                                       "to override automatic numbering")
    parser.add_argument("--csv", help="also write the results to this CSV file")
    args = parser.parse_args()

    manual = [int(v) - 1 for v in args.ring.split(",")] if args.ring else None
    results, failed = [], 0
    for path in args.files:
        try:
            result = analyse(path, manual)
        except (ValueError, OSError) as error:
            print(f"\n{os.path.basename(path)}\n  skipped: {error}", file=sys.stderr)
            failed += 1
            continue
        results.append(result)
        report(result)
    if args.csv and results:
        write_csv(results, args.csv)
        print(f"\nwrote {args.csv}")
    return 1 if failed and not results else 0


if __name__ == "__main__":
    sys.exit(main())
