# Feature tables

Thirteen parameter tables ship with DNAflexpy. Each maps every k-mer to a
number; a profile is the sliding-window mean of those numbers along a
sequence.

Everything on this page is read from the packaged table itself, so it cannot
drift out of date:

```python
from DNAflexpy import default_table

table = default_table()
print(len(table.features), "features")
print(table.kmer_len("DNaseI"), "-mer for DNaseI")
```

```
13 features
3 -mer for DNaseI
```

## Trinucleotide tables (k = 3, 64 entries)

| Feature | Range | Lowest | Highest | Values |
|---|---|---|---|---|
| `DNaseI` | −0.28 … 0.194 | `AAT` | `TCA` | float |
| `NPP` | 2 … 45 | `CAG` | `GCC` | integer |
| `bendabilityDNase` | 0 … 10 | `AAT` | `TCA` | mixed |
| `bendabilityConcensus` | 0.05 … 9.1 | `AAA` | `GCC` | mixed |

## Dinucleotide tables (k = 2, 16 entries)

| Feature | Range | Lowest | Highest | Values |
|---|---|---|---|---|
| `wedge` | 0.9 … 8.4 | `TA` | `AG` | float |
| `prop` | −17.3 … −6.7 | `AA` | `AC` | float |
| `freeen` | −2.24 … −0.58 | `GC` | `TA` | float |
| `gc` | 0.0 … 1.0 | `AA` | `CC` | float |
| `twistDisp` | 3.9 … 6.5 | `AA` | `CA` | float |
| `stiffness` | 1230 … 5500 | `GA` | `AT` | integer |
| `bendingStiffness` | 20 … 130 | `AT` | `CC` | integer |
| `mechen` | −3.42 … −0.81 | `GC` | `CA` | float |
| `trx` | 0 … 43 | `AT` | `CG` | integer |

!!! note "Add your primary references"
    The repository does not record where each table's numbers came from. The
    names point at the standard literature — DNase I cutting frequency,
    nucleosome positioning preference, wedge angle, propeller twist, stacking
    free energy, and so on — but this documentation will not put a citation
    next to a number it cannot verify. If you know the source paper for each
    table, adding it here would make the package considerably more useful to
    someone else. See [Citation](citation.md).

## What the ranges mean in practice

**Scales differ by four orders of magnitude.** `stiffness` reaches 5500 while
`gc` tops out at 1.0. That is why:

- a heatmap draws **one feature per figure** — they cannot share a colour scale
- `encode(normalize=True)` scales **each feature block separately** — one
  shared scale would flatten `gc` to nothing

**Four tables take negative values**: `DNaseI`, `prop`, `freeen` and `mechen`.
That is why `heatmap(order_rows="cv")` refuses them — the coefficient of
variation is `std/mean`, which inverts when the mean is negative.

**Six tables contain integers.** Four hold nothing else — `NPP`, `stiffness`,
`bendingStiffness`, `trx` — and two are mixed, holding both integers and
floats: `bendabilityDNase` and `bendabilityConcensus` (the "mixed" rows above).

The type is preserved exactly as written, never coerced. Writing `18.0` where
the original wrote `18` would break byte-equality with the frozen reference
implementation, and `assert a == b` would never catch it, because `18 == 18.0`
is `True` in Python. Only comparing the serialised bytes exposes it.

```python
gc = default_table().table("gc")
npp = default_table().table("NPP")
print(type(gc["AA"]).__name__, type(npp["AAA"]).__name__)
```

```
float int
```

## k-mer length is inferred, never declared

There is no feature → k mapping to maintain. `FeatureTable` reads the length
from the table's own keys and rejects any table whose keys are not all the same
length:

```python
print(default_table().kmer_len("gc"), default_table().kmer_len("NPP"))
```

```
2 3
```

The original package kept a hand-maintained map, and it had drifted: `freeen`,
`gc` and `mechen` were in the table but missing from the map, so they silently
produced rows of `None`. Inferring k removed the second place to update and
unlocked those three.

**To add a feature, add it to `DNAflexpy/data/lookupNEW.yaml` and nothing
else.**

## Using your own table

Pass a dict, a path to a YAML file, or a `FeatureTable`:

```python
from DNAflexpy import FlexProfiler

my_table = {"stickiness": {a + b: (i % 7) / 2 for i, (a, b) in enumerate(
    (x, y) for x in "ACGT" for y in "ACGT")}}

p = FlexProfiler("stickiness", window_size=4, lookup=my_table)
print(p.profile("ACGTACGTAC").round(3))
```

```
[1.833 2.5   1.667 2.    1.833 2.5   1.667]
```

Validation happens as the table loads:

- **mixed k-mer lengths inside one feature are rejected** — this is what makes
  inferring k safe
- keys outside A, C, G, T are rejected
- non-numeric values are rejected
- a table with fewer than 4<sup>k</sup> entries is allowed but **warns**;
  missing k-mers resolve to `NaN` and are masked

## Inspecting a table

```python
table = default_table()

print(sorted(table.features)[:4])
print(table.table("gc")["GC"], table.table("gc")["AT"])
```

```
['DNaseI', 'NPP', 'bendabilityConcensus', 'bendabilityDNase']
1.0 0.0
```

`.table(feature)` returns a read-only view. That is deliberate: the packaged
table is memoised and shared across every profiler in the process, so handing
out a mutable reference would let one caller corrupt every other.
