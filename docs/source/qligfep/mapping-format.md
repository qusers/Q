# Mapping JSON format

Network generators write a JSON object with `nodes` and `edges`. `setupFEP`
requires the `edges` array; visualizers and analyzers use additional metadata.

## Minimal input

```json
{
  "nodes": [
    {"name": "ligandA"},
    {"name": "ligandB"}
  ],
  "edges": [
    {
      "from": "ligandA",
      "to": "ligandB",
      "same_charge": true
    }
  ]
}
```

`from` and `to` names must match the ligand filenames and molecule names used by
the setup directory.

## Direction convention

An edge `ligandA` → `ligandB` represents the state change from A to B.
Experimental relative values stored as `ddg_value` follow:

$$
\Delta\Delta G_{A\rightarrow B} = G_B - G_A.
$$

Reversing an edge changes the sign. Keep direction unchanged between setup,
simulation, and analysis.

## Common node fields

| Field | Meaning |
|---|---|
| `name` | Canonical ligand identifier |
| `smiles` | Structure string |
| `formal_charge` | Molecular formal charge |
| `dg_value` | Optional experimental absolute binding free energy |

Network tools may preserve arbitrary SDF properties as extra node fields.

## Common edge fields

| Field | Meaning |
|---|---|
| `from`, `to` | Directed endpoint names |
| `weight` or score fields | Mapping/network score |
| `same_charge` | Whether endpoint formal charges match |
| `strict` | LoMap strict transformation flag, when available |
| `ddg_value` | Optional experimental to-minus-from value |
| `Q_ddG_avg` | Calculated relative binding free energy after analysis |
| `Q_ddG_sem` | Standard error written by analysis |
| `Q_ddG_std` | Standard deviation written by analysis |

## Discovery by `setupFEP`

If `-j` is omitted, setup checks in this order:

1. `mapping.json`;
2. `lomap.json`;
3. exactly one other `*.json` in the current directory.

Pass `-j` explicitly in reproducible scripts, especially after analysis has
created additional JSON files.

## Validation

Before setup, verify:

- every edge endpoint has `.sdf`, `.pdb`, `.lib`, and `.prm` files;
- names and case match exactly;
- there are no duplicate directed edges or accidental self-edges;
- intended ligands form the expected connected component(s);
- formal charge and `same_charge` agree;
- experimental values use one unit and direction convention;
- network score does not conceal chemically implausible mappings.

The analyzers write a copy rather than modifying the original mapping. Archive
the original input and the annotated output separately.
