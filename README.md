# cremer-pople

A single-file Python script that computes Cremer–Pople puckering parameters and endocyclic torsion angles for a six-membered ring from Cartesian coordinates, and names the nearest canonical conformations.

It was written to compare the ring conformations of computed stationary points (minima and transition structures), where the question is how far a ring has moved from one conformation toward another, and "chair" or "boat" is too coarse an answer.

## What it reports

For each structure:

- **Q, θ, φ** — the total puckering amplitude (Å) and the two angles that locate the conformation on the Cremer–Pople sphere, together with q₂ and q₃.
- **Endocyclic torsion angles** — about each ring bond, with the Klyne–Prelog sign.
- **Nearest canonical forms** — the three closest of the chair, half-chair, envelope, screw-boat, boat and twist-boat forms, each with its angular separation on the sphere.
- **Transition-structure metrics** (optional, see below) — forming and breaking bond lengths and the dihedral angles between the forming bond and the vicinal bonds.

## Requirements

Python 3.8 or later and NumPy.

```
pip install numpy
```

## Usage

Give the six ring atoms explicitly, as 1-based indices in ring order:

```
python cremer_pople.py examples/cyclohexane_chair.xyz --ring 1,2,3,4,5,6
```

```
cyclohexane_chair.xyz
  ring atoms C1-C6 (file indices): 1, 2, 3, 4, 5, 6
  Q = 0.562 A   theta = 180.0   phi = 93.6   (q2 = 0.000, q3 = -0.562)
  endocyclic torsions: C1-C2 +54.9  C2-C3 -54.9  C3-C4 +54.9  C4-C5 -54.9  C5-C6 +54.9  C6-C1 -54.9
  nearest canonical forms: chair (180.0, 0) at 0 deg; half-chair (129.2, 90) at 51 deg; half-chair (129.2, 30) at 51 deg
```

Several files can be processed at once, and the results written to a CSV file:

```
python cremer_pople.py *.xyz --csv puckering.csv
```

`--ring` applies to every file on the command line, so use it when the files share one atom ordering.

### Automatic numbering

Without `--ring`, the script numbers the ring itself for one class of molecule: a six-membered carbocycle carrying a side chain that ends in an amide, in which two adjacent ring carbons are bonded to, or in close contact with, nitrogen or oxygen. C1 is the ring carbon bearing the side chain; C3 and C4 are the two heteroatom-bearing carbons, C3 being two bonds from C1. The file indices chosen are printed, so the numbering can be checked. A file that does not fit this pattern is skipped with a reason, and `--ring` should be used instead.

For such a structure, if an oxygen atom of a separate fragment lies within 2.8 Å of C3 or C4, it is treated as a nucleophile and the transition-structure metrics are reported.

## Reading the output

| θ | Conformation |
|---|---|
| 0° or 180° | Chair (the two chairs are related by ring inversion) |
| 50.8° or 129.2° | Half-chair |
| 54.7° or 125.3° | Envelope |
| 67.5° or 112.5° | Screw-boat |
| 90° | Boat (φ = 0, 60, 120 …) and twist-boat (φ = 30, 90, 150 …) |

The intermediate values are those of Cremer and Pople and of Boeyens, derived for small displacements of a regular hexagon. They are reference points, not exact positions for a real ring. Most real structures fall between canonical forms, which is why the separations are reported.

## Cautions

1. **θ and φ depend on the numbering.** Moving the origin by one atom, or reversing the direction of numbering, converts θ to 180° − θ and shifts φ. The same chair gives θ = 180° with `--ring 1,2,3,4,5,6` and θ = 0° with `--ring 2,3,4,5,6,1`. Compare values only between structures numbered in the same way, and state the numbering whenever a value is reported.
2. **Enantiomers.** The mirror image of a structure, numbered in the same way, also gives 180° − θ. Check that the structures being compared belong to the same enantiomeric series; the signs of the torsion angles show this at a glance.
3. **φ near the poles.** φ is poorly defined when θ is close to 0° or 180°, as in the example above.
4. **Bond perception.** Bonds are assigned from interatomic distances (1.25 × the sum of covalent radii). This is used only for automatic numbering and for the transition-structure metrics; the puckering parameters themselves depend only on the six ring atoms.

## Method

The ring atoms are referred to their geometric centre. The mean plane is defined by eqs 8–10 of Cremer and Pople, the displacements z_j from that plane are obtained by eq 11, and q₂, φ₂ and q₃ follow from eqs 12–14. These are reported as the spherical polar set Q, θ, φ, with Q = (q₂² + q₃²)^½ and tan θ = q₂/q₃ (0° ≤ θ ≤ 180°). Each step in the source is commented with its equation number.

## References

1. Cremer, D.; Pople, J. A. General Definition of Ring Puckering Coordinates. *J. Am. Chem. Soc.* **1975**, *97*, 1354–1358. https://doi.org/10.1021/ja00839a011
2. Boeyens, J. C. A. The Conformation of Six-Membered Rings. *J. Cryst. Mol. Struct.* **1978**, *8*, 317–320. https://doi.org/10.1007/BF01200485

If you use the script, please cite these two papers for the method and the nomenclature.

## Licence

MIT. See `LICENSE`.
