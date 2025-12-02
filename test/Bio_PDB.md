# Documentation for `Bio.PDB`
**File Path:** `/home/nicai_zht/.local/share/mamba/envs/idr_llm/lib/python3.13/site-packages/Bio/PDB/__init__.py`

## Module Docstring
```text
Classes that deal with macromolecular crystal structures.

Includes: PDB and mmCIF parsers, a Structure class, a module to keep a local
copy of the PDB up-to-date, selective IO of PDB files, etc.

Original Author: Thomas Hamelryck.
Contributions by:
- Peter Cock
- Joe Greener
- Rob Miller
- Lenna X. Peterson
- Joao Rodrigues
- Kristian Rother
- Eric Talevich
- and many others.
```

## 🚦 Navigator: How to Drive
This section helps you understand how to run this library from the command line or entry points.
- ℹ️ **No Direct Entry Point**: This module seems to be a library intended for import, not direct execution.

### 🐍 Python API Usage (Inferred)
Since no CLI entry point was found, here are the likely **Python API entry points** for your script:

#### 🚀 Top 20 Recommended Entry Points
| Type | API | Description |
| :--- | :--- | :--- |
| `ƒ` | **Bio.PDB.m2rotaxis**(**m**) | Return angles, axis pair that corresponds to rotation matrix m. |
| `ƒ` | **Bio.PDB.make_dssp_dict**(**filename**) | DSSP dictionary mapping identifiers to properties. |
| `ƒ` | **Bio.PDB.parse_pdb_header**(**infile**) | Return the header lines of a pdb file as a dictionary. |
| `C` | **Bio.PDB.FastMMCIFParser**(structure_builder, auth_chains, auth_residues, QUIET) | Parse an MMCIF file and return a Structure object. |
| `C` | **Bio.PDB.MMCIFParser**(structure_builder, auth_chains, auth_residues, QUIET) | Parse a mmCIF file and return a Structure object. |
| `C` | **Bio.PDB.PDBMLParser**() | A parser for PDBML (PDB XML) files. See https://pdbml.wwpdb.org/. |
| `C` | **Bio.PDB.PDBParser**(PERMISSIVE, get_header, structure_builder, QUIET, is_pqr) | Parse a PDB file and return a Structure object. |
| `ƒ` | **Bio.PDB.extract**(**structure**, **chain_id**, **start**, **end**, **filename**) | Write out selected portion to filename. |
| `ƒ` | **Bio.PDB.is_nucleic**(**residue**, standard) | Return True if residue object/string is a nucleic acid. |
| `ƒ` | **Bio.PDB.rotaxis2m**(**theta**, **vector**) | Calculate left multiplying rotation matrix. |
| `ƒ` | **Bio.PDB.rotmat**(**p**, **q**) | Return a (left multiplying) matrix that rotates p onto q. |
| `ƒ` | **Bio.PDB.vector_to_axis**(**line**, **point**) | Vector to axis method. |
| `C` | **Bio.PDB.CEAligner**(window_size, max_gap) | Protein Structure Alignment by Combinatorial Extension. |
| `C` | **Bio.PDB.DSSP**(**model**, **in_file**, dssp, acc_array, file_type) | Run DSSP and parse secondary structure and accessibility. |
| `C` | **Bio.PDB.ExposureCN**(**model**, radius, offset) | Residue exposure as number of CA atoms around its CA atom. |
| `C` | **Bio.PDB.FragmentMapper**(**model**, lsize, flength, fdir) | Map polypeptides in a model to lists of representative fragments. |
| `C` | **Bio.PDB.HSExposureCA**(**model**, radius, offset) | Class to calculate HSE based on the approximate CA-CB vectors. |
| `C` | **Bio.PDB.MMCIFIO**() | Write a Structure object or a mmCIF dictionary as a mmCIF file. |
| `C` | **Bio.PDB.NeighborSearch**(**atom_list**, bucket_size) | Class for neighbor searching. |
| `C` | **Bio.PDB.PDBIO**(use_model_flag, is_pqr) | Write a Structure object (or a subset of a Structure object) as a PDB or PQR file. |

> **Note:** Bold parameters are required. Others are optional.

#### 🧩 Code Snippets (Auto-Generated)
```python
import Bio.PDB

# --- Top 20 Ranked Functions ---
# 1. m2rotaxis
result_1 = Bio.PDB.m2rotaxis(m=...)

# 2. make_dssp_dict
result_2 = Bio.PDB.make_dssp_dict(filename=...)

# 3. parse_pdb_header
result_3 = Bio.PDB.parse_pdb_header(infile=...)

# 4. extract
result_4 = Bio.PDB.extract(
    structure=...,
    chain_id=...,
    start=...,
    end=...,
    filename=...
)

# 5. is_nucleic
result_5 = Bio.PDB.is_nucleic(residue=...)

# 6. rotaxis2m
result_6 = Bio.PDB.rotaxis2m(theta=..., vector=...)

# 7. rotmat
result_7 = Bio.PDB.rotmat(p=..., q=...)

# 8. vector_to_axis
result_8 = Bio.PDB.vector_to_axis(line=..., point=...)

# 9. calc_angle
result_9 = Bio.PDB.calc_angle(v1=..., v2=..., v3=...)

# 10. get_surface
result_10 = Bio.PDB.get_surface(model=...)

# 11. refmat
result_11 = Bio.PDB.refmat(p=..., q=...)

# 12. rotaxis
result_12 = Bio.PDB.rotaxis(theta=..., vector=...)

# 13. calc_dihedral
result_13 = Bio.PDB.calc_dihedral(
    v1=...,
    v2=...,
    v3=...,
    v4=...
)

# 14. is_aa
result_14 = Bio.PDB.is_aa(residue=...)

# --- Top 20 Core Classes Initialization ---
# 1. FastMMCIFParser
fastmmcifparser = Bio.PDB.FastMMCIFParser()

# 2. MMCIFParser
mmcifparser = Bio.PDB.MMCIFParser()

# 3. PDBMLParser
pdbmlparser = Bio.PDB.PDBMLParser()

# 4. PDBParser
pdbparser = Bio.PDB.PDBParser()

# 5. CEAligner
cealigner = Bio.PDB.CEAligner()

# 6. DSSP
dssp = Bio.PDB.DSSP(model=..., in_file=...)

# 7. ExposureCN
exposurecn = Bio.PDB.ExposureCN(model=...)

# 8. FragmentMapper
fragmentmapper = Bio.PDB.FragmentMapper(model=...)

# 9. HSExposureCA
hsexposureca = Bio.PDB.HSExposureCA(model=...)

# 10. MMCIFIO
mmcifio = Bio.PDB.MMCIFIO()

# 11. NeighborSearch
neighborsearch = Bio.PDB.NeighborSearch(atom_list=...)

# 12. PDBIO
pdbio = Bio.PDB.PDBIO()

# 13. PDBList
pdblist = Bio.PDB.PDBList()

# 14. ShrakeRupley
shrakerupley = Bio.PDB.ShrakeRupley()

# 15. StructureAlignment
structurealignment = Bio.PDB.StructureAlignment(fasta_align=..., m1=..., m2=...)

# 16. Vector
vector = Bio.PDB.Vector(x=...)

# 17. CaPPBuilder
cappbuilder = Bio.PDB.CaPPBuilder()

# 18. HSExposureCB
hsexposurecb = Bio.PDB.HSExposureCB(model=...)

# 19. PPBuilder
ppbuilder = Bio.PDB.PPBuilder()

# 20. ResidueDepth
residuedepth = Bio.PDB.ResidueDepth(model=...)

```

_No explicit `argparse` configuration detected in the main module._


## 📊 Network & Architecture Analysis
### 🌍 Top 20 External Dependencies
| Library | Usage Count |
| :--- | :--- |
| **_frozen_importlib_external** | 43 |
| **_frozen_importlib** | 43 |
| **Bio** | 14 |
| **collections** | 8 |
| **urllib** | 5 |
| **_io** | 3 |
| **re** | 2 |
| **copy** | 1 |
| **concurrent** | 1 |
| **xml** | 1 |
| **datetime** | 1 |
| **itertools** | 1 |
| **numpy** | 1 |
| **numbers** | 1 |


### 🕸️ Network Metrics (Advanced)
#### 👑 Top 20 Modules by PageRank (Authority)
| Rank | Module | Score | Type | Role |
| :--- | :--- | :--- | :--- | :--- |
| 1 | `_frozen_importlib_external` | 0.1311 | External | External Lib |
| 2 | `_frozen_importlib` | 0.1311 | External | External Lib |
| 3 | `Bio` | 0.0481 | External | External Lib |
| 4 | `PDBExceptions` | 0.0460 | Internal | Utility / Core |
| 5 | `Entity` | 0.0262 | Internal | Data Processing |
| 6 | `collections` | 0.0243 | External | External Lib |
| 7 | `vectors` | 0.0204 | Internal | Data Processing |
| 8 | `Polypeptide` | 0.0198 | Internal | Utility / Core |
| 9 | `StructureBuilder` | 0.0193 | Internal | Data Processing |
| 10 | `internal_coords` | 0.0171 | Internal | Data Processing |
| 11 | `Atom` | 0.0160 | Internal | Data Processing |
| 12 | `Structure` | 0.0157 | Internal | Utility / Core |
| 13 | `PDBIO` | 0.0156 | Internal | Utility / Core |
| 14 | `AbstractPropertyMap` | 0.0154 | Internal | Utility / Core |
| 15 | `Residue` | 0.0150 | Internal | Utility / Core |
| 16 | `copy` | 0.0130 | External | External Lib |
| 17 | `urllib` | 0.0127 | External | External Lib |
| 18 | `parse_pdb_header` | 0.0122 | Internal | Utility / Core |
| 19 | `PDBParser` | 0.0121 | Internal | Data Processing |
| 20 | `re` | 0.0120 | External | External Lib |


### 🗺️ Dependency & Architecture Map
```mermaid
graph TD
    classDef core fill:#f96,stroke:#333,stroke-width:2px;
    classDef external fill:#9cf,stroke:#333,stroke-width:1px;
    id_33["PDB"] --> id_32["cealign"]
    class id_33 core;
    class id_32 core;
    id_33["PDB"] --> id_5["Polypeptide"]
    class id_33 core;
    class id_5 core;
    id_33["PDB"] --> id_12["DSSP"]
    class id_33 core;
    class id_12 core;
    id_33["PDB"] --> id_52["HSExposure"]
    class id_33 core;
    class id_52 core;
    id_33["PDB"] --> id_41["MMCIFParser"]
    class id_33 core;
    class id_41 core;
    id_33["PDB"] --> id_18["FragmentMapper"]
    class id_33 core;
    class id_18 core;
    id_33["PDB"] --> id_11["mmcifio"]
    class id_33 core;
    class id_11 core;
    id_33["PDB"] --> id_45["NeighborSearch"]
    class id_33 core;
    class id_45 core;
    id_33["PDB"] --> id_48["PDBIO"]
    class id_33 core;
    class id_48 core;
    id_33["PDB"] --> id_34["PDBList"]
    class id_33 core;
    class id_34 core;
    id_33["PDB"] --> id_3["PDBMLParser"]
    class id_33 core;
    class id_3 core;
    id_33["PDB"] --> id_36["PDBParser"]
    class id_33 core;
    class id_36 core;
    id_33["PDB"] --> id_27["ResidueDepth"]
    class id_33 core;
    class id_27 core;
    id_33["PDB"] --> id_30["SASA"]
    class id_33 core;
    class id_30 core;
    id_33["PDB"] --> id_53["StructureAlignment"]
    class id_33 core;
    class id_53 core;
    id_33["PDB"] --> id_42["Superimposer"]
    class id_33 core;
    class id_42 core;
    id_33["PDB"] --> id_55["vectors"]
    class id_33 core;
    class id_55 core;
    id_33["PDB"] -.-> id_46["_frozen_importlib_external"]
    class id_33 core;
    class id_46 external;
    id_33["PDB"] -.-> id_13["_frozen_importlib"]
    class id_33 core;
    class id_13 external;
    id_33["PDB"] --> id_15["Dice"]
    class id_33 core;
    class id_15 core;
    id_33["PDB"] --> id_40["parse_pdb_header"]
    class id_33 core;
    class id_40 core;
    id_32["cealign"] --> id_19["PDBExceptions"]
    class id_32 core;
    class id_19 core;
    id_32["cealign"] --> id_22["qcprot"]
    class id_32 core;
    class id_22 core;
    id_32["cealign"] -.-> id_46["_frozen_importlib_external"]
    class id_32 core;
    class id_46 external;
    id_32["cealign"] -.-> id_13["_frozen_importlib"]
    class id_32 core;
    class id_13 external;
    id_32["cealign"] --> id_0["ccealign"]
    class id_32 core;
    class id_0 core;
    id_5["Polypeptide"] --> id_19["PDBExceptions"]
    class id_5 core;
    class id_19 core;
    id_5["Polypeptide"] -.-> id_6["Bio"]
    class id_5 core;
    class id_6 external;
    id_5["Polypeptide"] -.-> id_46["_frozen_importlib_external"]
    class id_5 core;
    class id_46 external;
    id_5["Polypeptide"] -.-> id_13["_frozen_importlib"]
    class id_5 core;
    class id_13 external;
    id_5["Polypeptide"] --> id_55["vectors"]
    class id_5 core;
    class id_55 core;
    id_12["DSSP"] --> id_31["AbstractPropertyMap"]
    class id_12 core;
    class id_31 core;
    id_12["DSSP"] --> id_43["MMCIF2Dict"]
    class id_12 core;
    class id_43 core;
    id_12["DSSP"] --> id_19["PDBExceptions"]
    class id_12 core;
    class id_19 core;
    id_12["DSSP"] --> id_36["PDBParser"]
    class id_12 core;
    class id_36 core;
    id_12["DSSP"] -.-> id_54["_io"]
    class id_12 core;
    class id_54 external;
    id_12["DSSP"] -.-> id_46["_frozen_importlib_external"]
    class id_12 core;
    class id_46 external;
    id_12["DSSP"] -.-> id_13["_frozen_importlib"]
    class id_12 core;
    class id_13 external;
    id_12["DSSP"] -.-> id_2["re"]
    class id_12 core;
    class id_2 external;
    id_52["HSExposure"] --> id_31["AbstractPropertyMap"]
    class id_52 core;
    class id_31 core;
    id_52["HSExposure"] --> id_5["Polypeptide"]
    class id_52 core;
    class id_5 core;
    id_52["HSExposure"] -.-> id_46["_frozen_importlib_external"]
    class id_52 core;
    class id_46 external;
    id_52["HSExposure"] -.-> id_13["_frozen_importlib"]
    class id_52 core;
    class id_13 external;
    id_52["HSExposure"] --> id_55["vectors"]
    class id_52 core;
    class id_55 core;
    id_41["MMCIFParser"] --> id_43["MMCIF2Dict"]
    class id_41 core;
    class id_43 core;
    id_41["MMCIFParser"] --> id_19["PDBExceptions"]
    class id_41 core;
    class id_19 core;
    id_41["MMCIFParser"] --> id_35["StructureBuilder"]
    class id_41 core;
    class id_35 core;
    id_41["MMCIFParser"] -.-> id_46["_frozen_importlib_external"]
    class id_41 core;
    class id_46 external;
    id_41["MMCIFParser"] -.-> id_13["_frozen_importlib"]
    class id_41 core;
    class id_13 external;
    id_41["MMCIFParser"] -.-> id_6["Bio"]
    class id_41 core;
    class id_6 external;
    id_18["FragmentMapper"] --> id_19["PDBExceptions"]
    class id_18 core;
    class id_19 core;
    id_18["FragmentMapper"] --> id_5["Polypeptide"]
    class id_18 core;
    class id_5 core;
    id_18["FragmentMapper"] -.-> id_6["Bio"]
    class id_18 core;
    class id_6 external;
    id_18["FragmentMapper"] -.-> id_46["_frozen_importlib_external"]
    class id_18 core;
    class id_46 external;
    id_18["FragmentMapper"] -.-> id_13["_frozen_importlib"]
    class id_18 core;
    class id_13 external;
    id_11["mmcifio"] --> id_48["PDBIO"]
    class id_11 core;
    class id_48 core;
    id_11["mmcifio"] --> id_35["StructureBuilder"]
    class id_11 core;
    class id_35 core;
    id_11["mmcifio"] -.-> id_46["_frozen_importlib_external"]
    class id_11 core;
    class id_46 external;
    id_11["mmcifio"] -.-> id_13["_frozen_importlib"]
    class id_11 core;
    class id_13 external;
    id_11["mmcifio"] -.-> id_38["collections"]
    class id_11 core;
    class id_38 external;
    id_45["NeighborSearch"] --> id_19["PDBExceptions"]
    class id_45 core;
    class id_19 core;
    id_45["NeighborSearch"] -.-> id_46["_frozen_importlib_external"]
    class id_45 core;
    class id_46 external;
    id_45["NeighborSearch"] -.-> id_13["_frozen_importlib"]
    class id_45 core;
    class id_13 external;
    id_45["NeighborSearch"] --> id_9["Selection"]
    class id_45 core;
    class id_9 core;
    id_48["PDBIO"] -.-> id_6["Bio"]
    class id_48 core;
    class id_6 external;
    id_48["PDBIO"] --> id_19["PDBExceptions"]
    class id_48 core;
    class id_19 core;
    id_48["PDBIO"] --> id_35["StructureBuilder"]
    class id_48 core;
    class id_35 core;
    id_48["PDBIO"] -.-> id_46["_frozen_importlib_external"]
    class id_48 core;
    class id_46 external;
    id_48["PDBIO"] -.-> id_13["_frozen_importlib"]
    class id_48 core;
    class id_13 external;
    id_34["PDBList"] -.-> id_23["urllib"]
    class id_34 core;
    class id_23 external;
    id_34["PDBList"] -.-> id_4["concurrent"]
    class id_34 core;
    class id_4 external;
    id_34["PDBList"] -.-> id_46["_frozen_importlib_external"]
    class id_34 core;
    class id_46 external;
    id_34["PDBList"] -.-> id_13["_frozen_importlib"]
    class id_34 core;
    class id_13 external;
    id_3["PDBMLParser"] -.-> id_26["xml"]
    class id_3 core;
    class id_26 external;
    id_3["PDBMLParser"] --> id_20["Structure"]
    class id_3 core;
    class id_20 core;
    id_3["PDBMLParser"] --> id_35["StructureBuilder"]
    class id_3 core;
    class id_35 core;
    id_3["PDBMLParser"] -.-> id_46["_frozen_importlib_external"]
    class id_3 core;
    class id_46 external;
    id_3["PDBMLParser"] -.-> id_13["_frozen_importlib"]
    class id_3 core;
    class id_13 external;
    id_36["PDBParser"] --> id_19["PDBExceptions"]
    class id_36 core;
    class id_19 core;
    id_36["PDBParser"] --> id_35["StructureBuilder"]
    class id_36 core;
    class id_35 core;
    id_36["PDBParser"] -.-> id_46["_frozen_importlib_external"]
    class id_36 core;
    class id_46 external;
    id_36["PDBParser"] -.-> id_13["_frozen_importlib"]
    class id_36 core;
    class id_13 external;
    id_36["PDBParser"] --> id_40["parse_pdb_header"]
    class id_36 core;
    class id_40 core;
    id_36["PDBParser"] -.-> id_6["Bio"]
    class id_36 core;
    class id_6 external;
    id_27["ResidueDepth"] --> id_31["AbstractPropertyMap"]
    class id_27 core;
    class id_31 core;
    id_27["ResidueDepth"] -.-> id_6["Bio"]
    class id_27 core;
    class id_6 external;
    id_27["ResidueDepth"] --> id_36["PDBParser"]
    class id_27 core;
    class id_36 core;
    id_27["ResidueDepth"] -.-> id_46["_frozen_importlib_external"]
    class id_27 core;
    class id_46 external;
    id_27["ResidueDepth"] -.-> id_13["_frozen_importlib"]
    class id_27 core;
    class id_13 external;
    id_27["ResidueDepth"] --> id_5["Polypeptide"]
    class id_27 core;
    class id_5 core;
    id_30["SASA"] -.-> id_38["collections"]
    class id_30 core;
    class id_38 external;
    id_30["SASA"] -.-> id_46["_frozen_importlib_external"]
    class id_30 core;
    class id_46 external;
    id_30["SASA"] -.-> id_13["_frozen_importlib"]
    class id_30 core;
    class id_13 external;
    id_53["StructureAlignment"] -.-> id_46["_frozen_importlib_external"]
    class id_53 core;
    class id_46 external;
    id_53["StructureAlignment"] -.-> id_13["_frozen_importlib"]
    class id_53 core;
    class id_13 external;
    id_53["StructureAlignment"] --> id_5["Polypeptide"]
    class id_53 core;
    class id_5 core;
    id_42["Superimposer"] --> id_19["PDBExceptions"]
    class id_42 core;
    class id_19 core;
    id_42["Superimposer"] -.-> id_6["Bio"]
    class id_42 core;
    class id_6 external;
    id_42["Superimposer"] -.-> id_46["_frozen_importlib_external"]
    class id_42 core;
    class id_46 external;
    id_42["Superimposer"] -.-> id_13["_frozen_importlib"]
    class id_42 core;
    class id_13 external;
    id_55["vectors"] -.-> id_46["_frozen_importlib_external"]
    class id_55 core;
    class id_46 external;
    id_55["vectors"] -.-> id_13["_frozen_importlib"]
    class id_55 core;
    class id_13 external;
    id_15["Dice"] -.-> id_6["Bio"]
    class id_15 core;
    class id_6 external;
    id_15["Dice"] --> id_48["PDBIO"]
    class id_15 core;
    class id_48 core;
    id_15["Dice"] -.-> id_46["_frozen_importlib_external"]
    class id_15 core;
    class id_46 external;
    id_15["Dice"] -.-> id_13["_frozen_importlib"]
    class id_15 core;
    class id_13 external;
    id_15["Dice"] -.-> id_2["re"]
    class id_15 core;
    class id_2 external;
    id_40["parse_pdb_header"] -.-> id_46["_frozen_importlib_external"]
    class id_40 core;
    class id_46 external;
    id_40["parse_pdb_header"] -.-> id_13["_frozen_importlib"]
    class id_40 core;
    class id_13 external;
    id_40["parse_pdb_header"] -.-> id_38["collections"]
    class id_40 core;
    class id_38 external;
    id_31["AbstractPropertyMap"] -.-> id_46["_frozen_importlib_external"]
    class id_31 core;
    class id_46 external;
    id_31["AbstractPropertyMap"] -.-> id_13["_frozen_importlib"]
    class id_31 core;
    class id_13 external;
    id_49["Atom"] --> id_16["Entity"]
    class id_49 core;
    class id_16 core;
    id_49["Atom"] --> id_19["PDBExceptions"]
    class id_49 core;
    class id_19 core;
    id_49["Atom"] --> id_47["Residue"]
    class id_49 core;
    class id_47 core;
    id_49["Atom"] --> id_55["vectors"]
    class id_49 core;
    class id_55 core;
    id_49["Atom"] -.-> id_46["_frozen_importlib_external"]
    class id_49 core;
    class id_46 external;
    id_49["Atom"] -.-> id_13["_frozen_importlib"]
    class id_49 core;
    class id_13 external;
    id_16["Entity"] -.-> id_6["Bio"]
    class id_16 core;
    class id_6 external;
    id_16["Entity"] --> id_19["PDBExceptions"]
    class id_16 core;
    class id_19 core;
    id_16["Entity"] -.-> id_46["_frozen_importlib_external"]
    class id_16 core;
    class id_46 external;
    id_16["Entity"] -.-> id_13["_frozen_importlib"]
    class id_16 core;
    class id_13 external;
    id_16["Entity"] -.-> id_24["copy"]
    class id_16 core;
    class id_24 external;
    id_16["Entity"] -.-> id_38["collections"]
    class id_16 core;
    class id_38 external;
    id_19["PDBExceptions"] -.-> id_6["Bio"]
    class id_19 core;
    class id_6 external;
    id_19["PDBExceptions"] -.-> id_46["_frozen_importlib_external"]
    class id_19 core;
    class id_46 external;
    id_19["PDBExceptions"] -.-> id_13["_frozen_importlib"]
    class id_19 core;
    class id_13 external;
    id_47["Residue"] --> id_16["Entity"]
    class id_47 core;
    class id_16 core;
    id_47["Residue"] --> id_19["PDBExceptions"]
    class id_47 core;
    class id_19 core;
    id_47["Residue"] -.-> id_46["_frozen_importlib_external"]
    class id_47 core;
    class id_46 external;
    id_47["Residue"] -.-> id_13["_frozen_importlib"]
    class id_47 core;
    class id_13 external;
    id_44["Chain"] --> id_16["Entity"]
    class id_44 core;
    class id_16 core;
    id_44["Chain"] --> id_28["internal_coords"]
    class id_44 core;
    class id_28 core;
    id_44["Chain"] -.-> id_46["_frozen_importlib_external"]
    class id_44 core;
    class id_46 external;
    id_44["Chain"] -.-> id_13["_frozen_importlib"]
    class id_44 core;
    class id_13 external;
    id_28["internal_coords"] --> id_49["Atom"]
    class id_28 core;
    class id_49 core;
    id_28["internal_coords"] -.-> id_14["numpy"]
    class id_28 core;
    class id_14 external;
    id_28["internal_coords"] -.-> id_7["numbers"]
    class id_28 core;
    class id_7 external;
    id_28["internal_coords"] -.-> id_46["_frozen_importlib_external"]
    class id_28 core;
    class id_46 external;
    id_28["internal_coords"] -.-> id_13["_frozen_importlib"]
    class id_28 core;
    class id_13 external;
    id_28["internal_coords"] --> id_55["vectors"]
    class id_28 core;
    class id_55 core;
    id_28["internal_coords"] -.-> id_38["collections"]
    class id_28 core;
    class id_38 external;
    id_43["MMCIF2Dict"] -.-> id_46["_frozen_importlib_external"]
    class id_43 core;
    class id_46 external;
    id_43["MMCIF2Dict"] -.-> id_13["_frozen_importlib"]
    class id_43 core;
    class id_13 external;
    id_43["MMCIF2Dict"] -.-> id_6["Bio"]
    class id_43 core;
    class id_6 external;
    id_35["StructureBuilder"] --> id_49["Atom"]
    class id_35 core;
    class id_49 core;
    id_35["StructureBuilder"] --> id_44["Chain"]
    class id_35 core;
    class id_44 core;
    id_35["StructureBuilder"] --> id_47["Residue"]
    class id_35 core;
    class id_47 core;
    id_35["StructureBuilder"] --> id_21["Model"]
    class id_35 core;
    class id_21 core;
    id_35["StructureBuilder"] --> id_19["PDBExceptions"]
    class id_35 core;
    class id_19 core;
    id_35["StructureBuilder"] --> id_20["Structure"]
    class id_35 core;
    class id_20 core;
    id_35["StructureBuilder"] -.-> id_46["_frozen_importlib_external"]
    class id_35 core;
    class id_46 external;
    id_35["StructureBuilder"] -.-> id_13["_frozen_importlib"]
    class id_35 core;
    class id_13 external;
    id_21["Model"] --> id_16["Entity"]
    class id_21 core;
    class id_16 core;
    id_21["Model"] --> id_28["internal_coords"]
    class id_21 core;
    class id_28 core;
    id_21["Model"] -.-> id_46["_frozen_importlib_external"]
    class id_21 core;
    class id_46 external;
    id_21["Model"] -.-> id_13["_frozen_importlib"]
    class id_21 core;
    class id_13 external;
    id_37["NACCESS"] --> id_31["AbstractPropertyMap"]
    class id_37 core;
    class id_31 core;
    id_37["NACCESS"] --> id_48["PDBIO"]
    class id_37 core;
    class id_48 core;
    id_37["NACCESS"] -.-> id_46["_frozen_importlib_external"]
    class id_37 core;
    class id_46 external;
    id_37["NACCESS"] -.-> id_13["_frozen_importlib"]
    class id_37 core;
    class id_13 external;
    id_9["Selection"] --> id_49["Atom"]
    class id_9 core;
    class id_49 core;
    id_9["Selection"] --> id_16["Entity"]
    class id_9 core;
    class id_16 core;
    id_9["Selection"] --> id_19["PDBExceptions"]
    class id_9 core;
    class id_19 core;
    id_9["Selection"] -.-> id_46["_frozen_importlib_external"]
    class id_9 core;
    class id_46 external;
    id_9["Selection"] -.-> id_13["_frozen_importlib"]
    class id_9 core;
    class id_13 external;
    id_20["Structure"] --> id_16["Entity"]
    class id_20 core;
    class id_16 core;
    id_20["Structure"] -.-> id_46["_frozen_importlib_external"]
    class id_20 core;
    class id_46 external;
    id_20["Structure"] -.-> id_13["_frozen_importlib"]
    class id_20 core;
    class id_13 external;
    id_39["PICIO"] --> id_28["internal_coords"]
    class id_39 core;
    class id_28 core;
    id_39["PICIO"] --> id_19["PDBExceptions"]
    class id_39 core;
    class id_19 core;
    id_39["PICIO"] --> id_47["Residue"]
    class id_39 core;
    class id_47 core;
    id_39["PICIO"] -.-> id_54["_io"]
    class id_39 core;
    class id_54 external;
    id_39["PICIO"] --> id_20["Structure"]
    class id_39 core;
    class id_20 core;
    id_39["PICIO"] --> id_35["StructureBuilder"]
    class id_39 core;
    class id_35 core;
    id_39["PICIO"] -.-> id_46["_frozen_importlib_external"]
    class id_39 core;
    class id_46 external;
    id_39["PICIO"] -.-> id_13["_frozen_importlib"]
    class id_39 core;
    class id_13 external;
    id_39["PICIO"] --> id_40["parse_pdb_header"]
    class id_39 core;
    class id_40 core;
    id_39["PICIO"] -.-> id_6["Bio"]
    class id_39 core;
    class id_6 external;
    id_39["PICIO"] -.-> id_50["datetime"]
    class id_39 core;
    class id_50 external;
    id_8["PSEA"] -.-> id_46["_frozen_importlib_external"]
    class id_8 core;
    class id_46 external;
    id_8["PSEA"] -.-> id_13["_frozen_importlib"]
    class id_8 core;
    class id_13 external;
    id_8["PSEA"] --> id_5["Polypeptide"]
    class id_8 core;
    class id_5 core;
    id_29["SCADIO"] --> id_28["internal_coords"]
    class id_29 core;
    class id_28 core;
    id_29["SCADIO"] --> id_19["PDBExceptions"]
    class id_29 core;
    class id_19 core;
    id_29["SCADIO"] -.-> id_46["_frozen_importlib_external"]
    class id_29 core;
    class id_46 external;
    id_29["SCADIO"] -.-> id_13["_frozen_importlib"]
    class id_29 core;
    class id_13 external;
    id_29["SCADIO"] -.-> id_6["Bio"]
    class id_29 core;
    class id_6 external;
    id_29["SCADIO"] --> id_55["vectors"]
    class id_29 core;
    class id_55 core;
    id_51["_bcif_helper"] -.-> id_46["_frozen_importlib_external"]
    class id_51 core;
    class id_46 external;
    id_51["_bcif_helper"] -.-> id_13["_frozen_importlib"]
    class id_51 core;
    class id_13 external;
    id_17["alphafold_db"] -.-> id_38["collections"]
    class id_17 core;
    class id_38 external;
    id_17["alphafold_db"] --> id_41["MMCIFParser"]
    class id_17 core;
    class id_41 core;
    id_17["alphafold_db"] --> id_20["Structure"]
    class id_17 core;
    class id_20 core;
    id_17["alphafold_db"] -.-> id_46["_frozen_importlib_external"]
    class id_17 core;
    class id_46 external;
    id_17["alphafold_db"] -.-> id_13["_frozen_importlib"]
    class id_17 core;
    class id_13 external;
    id_17["alphafold_db"] -.-> id_23["urllib"]
    class id_17 core;
    class id_23 external;
    id_0["ccealign"] -.-> id_46["_frozen_importlib_external"]
    class id_0 core;
    class id_46 external;
    id_0["ccealign"] -.-> id_13["_frozen_importlib"]
    class id_0 core;
    class id_13 external;
    id_22["qcprot"] --> id_19["PDBExceptions"]
    class id_22 core;
    class id_19 core;
    id_22["qcprot"] -.-> id_46["_frozen_importlib_external"]
    class id_22 core;
    class id_46 external;
    id_22["qcprot"] -.-> id_13["_frozen_importlib"]
    class id_22 core;
    class id_13 external;
    id_1["ic_data"] -.-> id_46["_frozen_importlib_external"]
    class id_1 core;
    class id_46 external;
    id_1["ic_data"] -.-> id_13["_frozen_importlib"]
    class id_1 core;
    class id_13 external;
    id_10["ic_rebuild"] --> id_49["Atom"]
    class id_10 core;
    class id_49 core;
    id_10["ic_rebuild"] --> id_44["Chain"]
    class id_10 core;
    class id_44 core;
    id_10["ic_rebuild"] --> id_47["Residue"]
    class id_10 core;
    class id_47 core;
    id_10["ic_rebuild"] --> id_28["internal_coords"]
    class id_10 core;
    class id_28 core;
    id_10["ic_rebuild"] --> id_21["Model"]
    class id_10 core;
    class id_21 core;
    id_10["ic_rebuild"] --> id_19["PDBExceptions"]
    class id_10 core;
    class id_19 core;
    id_10["ic_rebuild"] --> id_48["PDBIO"]
    class id_10 core;
    class id_48 core;
    id_10["ic_rebuild"] -.-> id_54["_io"]
    class id_10 core;
    class id_54 external;
    id_10["ic_rebuild"] --> id_20["Structure"]
    class id_10 core;
    class id_20 core;
    id_10["ic_rebuild"] -.-> id_46["_frozen_importlib_external"]
    class id_10 core;
    class id_46 external;
    id_10["ic_rebuild"] -.-> id_13["_frozen_importlib"]
    class id_10 core;
    class id_13 external;
    id_10["ic_rebuild"] -.-> id_6["Bio"]
    class id_10 core;
    class id_6 external;
    id_10["ic_rebuild"] --> id_39["PICIO"]
    class id_10 core;
    class id_39 core;
    id_10["ic_rebuild"] -.-> id_25["itertools"]
    class id_10 core;
    class id_25 external;
    id_56["kdtrees"] -.-> id_46["_frozen_importlib_external"]
    class id_56 core;
    class id_46 external;
    id_56["kdtrees"] -.-> id_13["_frozen_importlib"]
    class id_56 core;
    class id_13 external;
```

## 🚀 Global Execution Flow & Extraction Guide
This graph visualizes how data flows between functions across the entire project.
```mermaid
graph TD
    classDef main fill:#f9f,stroke:#333,stroke-width:2px;
    classDef func fill:#fff,stroke:#333,stroke-width:1px;
    f_2["__contains__"] -->|id| f_96["_translate_id"]
    class f_2 func;
    class f_96 func;
    f_6["__getitem__"] -->|key| f_96["_translate_id"]
    class f_6 func;
    class f_96 func;
    f_9["__init__"] --> f_329["upper"]
    class f_9 func;
    class f_329 func;
    f_9["__init__"] -->|element| f_20["_assign_element"]
    class f_9 func;
    class f_20 func;
    f_9["__init__"] --> f_19["_assign_atom_mass"]
    class f_9 func;
    class f_19 func;
    f_7["__gt__"] --> f_173["get"]
    class f_7 func;
    class f_173 func;
    f_5["__ge__"] --> f_173["get"]
    class f_5 func;
    class f_173 func;
    f_12["__lt__"] --> f_173["get"]
    class f_12 func;
    class f_173 func;
    f_11["__le__"] --> f_173["get"]
    class f_11 func;
    class f_173 func;
    f_8["__hash__"] --> f_185["get_full_id"]
    class f_8 func;
    class f_185 func;
    f_20["_assign_element"] --> f_318["strip"]
    class f_20 func;
    class f_318 func;
    f_14["__repr__"] --> f_188["get_id"]
    class f_14 func;
    class f_188 func;
    f_306["set_parent"] --> f_185["get_full_id"]
    class f_306 func;
    class f_185 func;
    f_143["copy"] --> f_148["detach_parent"]
    class f_143 func;
    class f_148 func;
    f_143["copy"] --> f_298["set_coord"]
    class f_143 func;
    class f_298 func;
    f_143["copy"] --> f_184["get_coord"]
    class f_143 func;
    class f_184 func;
    f_10["__iter__"] --> f_153["disordered_get_list"]
    class f_10 func;
    class f_153 func;
    f_136["center_of_mass"] --> f_153["disordered_get_list"]
    class f_136 func;
    class f_153 func;
    f_153["disordered_get_list"] --> f_333["values"]
    class f_153 func;
    class f_333 func;
    f_150["disordered_add"] --> f_169["flag_disorder"]
    class f_150 func;
    class f_169 func;
    f_150["disordered_add"] --> f_195["get_parent"]
    class f_150 func;
    class f_195 func;
    f_150["disordered_add"] -->|residue| f_306["set_parent"]
    class f_150 func;
    class f_306 func;
    f_150["disordered_add"] --> f_177["get_altloc"]
    class f_150 func;
    class f_177 func;
    f_150["disordered_add"] --> f_194["get_occupancy"]
    class f_150 func;
    class f_194 func;
    f_150["disordered_add"] -->|altloc| f_156["disordered_select"]
    class f_150 func;
    class f_156 func;
    f_155["disordered_remove"] --> f_148["detach_parent"]
    class f_155 func;
    class f_148 func;
    f_155["disordered_remove"] --> f_333["values"]
    class f_155 func;
    class f_333 func;
    f_155["disordered_remove"] --> f_156["disordered_select"]
    class f_155 func;
    class f_156 func;
    f_6["__getitem__"] -->|id| f_96["_translate_id"]
    class f_6 func;
    class f_96 func;
    f_4["__delitem__"] -->|id| f_96["_translate_id"]
    class f_4 func;
    class f_96 func;
    f_215["get_unpacked_list"] --> f_191["get_list"]
    class f_215 func;
    class f_191 func;
    f_215["get_unpacked_list"] --> f_233["is_disordered"]
    class f_215 func;
    class f_233 func;
    f_215["get_unpacked_list"] --> f_153["disordered_get_list"]
    class f_215 func;
    class f_153 func;
    f_217["has_id"] -->|id| f_96["_translate_id"]
    class f_217 func;
    class f_96 func;
    f_180["get_atoms"] --> f_199["get_residues"]
    class f_180 func;
    class f_199 func;
    f_335["version"] -->|int| f_243["map"]
    class f_335 func;
    class f_243 func;
    f_335["version"] --> f_314["split"]
    class f_335 func;
    class f_314 func;
    f_163["dssp_dict_from_pdb_file"] -->|dssp_version| f_335["version"]
    class f_163 func;
    class f_335 func;
    f_163["dssp_dict_from_pdb_file"] --> f_335["version"]
    class f_163 func;
    class f_335 func;
    f_163["dssp_dict_from_pdb_file"] --> f_318["strip"]
    class f_163 func;
    class f_318 func;
    f_163["dssp_dict_from_pdb_file"] --> f_61["_make_dssp_dict"]
    class f_163 func;
    class f_61 func;
    f_241["make_dssp_dict"] -->|handle| f_61["_make_dssp_dict"]
    class f_241 func;
    class f_61 func;
    f_61["_make_dssp_dict"] --> f_314["split"]
    class f_61 func;
    class f_314 func;
    f_61["_make_dssp_dict"] --> f_167["find"]
    class f_61 func;
    class f_167 func;
    f_9["__init__"] -->|version_string| f_285["search"]
    class f_9 func;
    class f_285 func;
    f_9["__init__"] -->|in_file<br>dssp<br>dssp_version| f_163["dssp_dict_from_pdb_file"]
    class f_9 func;
    class f_163 func;
    f_9["__init__"] -->|in_file| f_241["make_dssp_dict"]
    class f_9 func;
    class f_241 func;
    f_9["__init__"] -->|dssp_version| f_335["version"]
    class f_9 func;
    class f_335 func;
    f_9["__init__"] --> f_335["version"]
    class f_9 func;
    class f_335 func;
    f_9["__init__"] --> f_233["is_disordered"]
    class f_9 func;
    class f_233 func;
    f_9["__init__"] --> f_152["disordered_get_id_list"]
    class f_9 func;
    class f_152 func;
    f_9["__init__"] --> f_177["get_altloc"]
    class f_9 func;
    class f_177 func;
    f_9["__init__"] --> f_191["get_list"]
    class f_9 func;
    class f_191 func;
    f_9["__init__"] -->|rk| f_156["disordered_select"]
    class f_9 func;
    class f_156 func;
    f_9["__init__"] --> f_156["disordered_select"]
    class f_9 func;
    class f_156 func;
    f_9["__init__"] --> f_215["get_unpacked_list"]
    class f_9 func;
    class f_215 func;
    f_9["__init__"] --> f_200["get_resname"]
    class f_9 func;
    class f_200 func;
    f_9["__init__"] -->|resname| f_173["get"]
    class f_9 func;
    class f_173 func;
    f_106["accept_model"] --> f_188["get_id"]
    class f_106 func;
    class f_188 func;
    f_105["accept_chain"] --> f_188["get_id"]
    class f_105 func;
    class f_188 func;
    f_107["accept_residue"] --> f_188["get_id"]
    class f_107 func;
    class f_188 func;
    f_104["accept_atom"] --> f_188["get_id"]
    class f_104 func;
    class f_188 func;
    f_166["extract"] -->|structure| f_311["set_structure"]
    class f_166 func;
    class f_311 func;
    f_166["extract"] -->|filename<br>sel| f_284["save"]
    class f_166 func;
    class f_284 func;
    f_4["__delitem__"] -->|id| f_147["detach_child"]
    class f_4 func;
    class f_147 func;
    f_82["_reset_full_id"] --> f_38["_generate_full_id"]
    class f_82 func;
    class f_38 func;
    f_38["_generate_full_id"] --> f_188["get_id"]
    class f_38 func;
    class f_188 func;
    f_38["_generate_full_id"] --> f_195["get_parent"]
    class f_38 func;
    class f_195 func;
    f_219["id"] --> f_82["_reset_full_id"]
    class f_219 func;
    class f_82 func;
    f_306["set_parent"] --> f_82["_reset_full_id"]
    class f_306 func;
    class f_82 func;
    f_147["detach_child"] --> f_148["detach_parent"]
    class f_147 func;
    class f_148 func;
    f_108["add"] --> f_188["get_id"]
    class f_108 func;
    class f_188 func;
    f_108["add"] -->|entity_id| f_217["has_id"]
    class f_108 func;
    class f_217 func;
    f_108["add"] -->|self| f_306["set_parent"]
    class f_108 func;
    class f_306 func;
    f_229["insert"] --> f_188["get_id"]
    class f_229 func;
    class f_188 func;
    f_229["insert"] -->|entity_id| f_217["has_id"]
    class f_229 func;
    class f_217 func;
    f_229["insert"] -->|self| f_306["set_parent"]
    class f_229 func;
    class f_306 func;
    f_191["get_list"] --> f_143["copy"]
    class f_191 func;
    class f_143 func;
    f_185["get_full_id"] --> f_38["_generate_full_id"]
    class f_185 func;
    class f_38 func;
    f_322["transform"] --> f_191["get_list"]
    class f_322 func;
    class f_191 func;
    f_136["center_of_mass"] --> f_215["get_unpacked_list"]
    class f_136 func;
    class f_215 func;
    f_143["copy"] --> f_108["add"]
    class f_143 func;
    class f_108 func;
    f_143["copy"] --> f_153["disordered_get_list"]
    class f_143 func;
    class f_153 func;
    f_143["copy"] --> f_150["disordered_add"]
    class f_143 func;
    class f_150 func;
    f_317["strictly_equals"] --> f_188["get_id"]
    class f_317 func;
    class f_188 func;
    f_317["strictly_equals"] --> f_237["keys"]
    class f_317 func;
    class f_237 func;
    f_148["detach_parent"] --> f_153["disordered_get_list"]
    class f_148 func;
    class f_153 func;
    f_306["set_parent"] --> f_153["disordered_get_list"]
    class f_306 func;
    class f_153 func;
    f_78["_read_fragments"] --> f_314["split"]
    class f_78 func;
    class f_314 func;
    f_78["_read_fragments"] -->|coord| f_109["add_residue"]
    class f_78 func;
    class f_109 func;
    f_15["__sub__"] --> f_288["set"]
    class f_15 func;
    class f_288 func;
    f_15["__sub__"] --> f_281["run"]
    class f_15 func;
    class f_281 func;
    f_15["__sub__"] --> f_201["get_rms"]
    class f_15 func;
    class f_201 func;
    f_62["_make_fragment_list"] --> f_200["get_resname"]
    class f_62 func;
    class f_200 func;
    f_62["_make_fragment_list"] --> f_217["has_id"]
    class f_62 func;
    class f_217 func;
    f_62["_make_fragment_list"] --> f_233["is_disordered"]
    class f_62 func;
    class f_233 func;
    f_62["_make_fragment_list"] --> f_184["get_coord"]
    class f_62 func;
    class f_184 func;
    f_62["_make_fragment_list"] -->|resname<br>ca_coord| f_109["add_residue"]
    class f_62 func;
    class f_109 func;
    f_64["_map_fragment_list"] --> f_313["sort"]
    class f_64 func;
    class f_313 func;
    f_9["__init__"] -->|lsize<br>flength<br>fdir| f_78["_read_fragments"]
    class f_9 func;
    class f_78 func;
    f_9["__init__"] --> f_63["_map"]
    class f_9 func;
    class f_63 func;
    f_63["_map"] -->|model| f_131["build_peptides"]
    class f_63 func;
    class f_131 func;
    f_63["_map"] -->|pp| f_62["_make_fragment_list"]
    class f_63 func;
    class f_62 func;
    f_63["_map"] -->|flist| f_64["_map_fragment_list"]
    class f_63 func;
    class f_64 func;
    f_9["__init__"] -->|model| f_131["build_peptides"]
    class f_9 func;
    class f_131 func;
    f_9["__init__"] -->|r1<br>r2<br>r3| f_45["_get_cb"]
    class f_9 func;
    class f_45 func;
    f_9["__init__"] --> f_216["get_vector"]
    class f_9 func;
    class f_216 func;
    f_9["__init__"] -->|ro| f_231["is_aa"]
    class f_9 func;
    class f_231 func;
    f_9["__init__"] --> f_217["has_id"]
    class f_9 func;
    class f_217 func;
    f_9["__init__"] --> f_249["norm"]
    class f_9 func;
    class f_249 func;
    f_9["__init__"] -->|pcb| f_115["angle"]
    class f_9 func;
    class f_115 func;
    f_9["__init__"] --> f_188["get_id"]
    class f_9 func;
    class f_188 func;
    f_9["__init__"] --> f_195["get_parent"]
    class f_9 func;
    class f_195 func;
    f_46["_get_gly_cb_vector"] --> f_216["get_vector"]
    class f_46 func;
    class f_216 func;
    f_46["_get_gly_cb_vector"] -->|c_v| f_277["rotaxis"]
    class f_46 func;
    class f_277 func;
    f_46["_get_gly_cb_vector"] -->|rot| f_238["left_multiply"]
    class f_46 func;
    class f_238 func;
    f_45["_get_cb"] --> f_216["get_vector"]
    class f_45 func;
    class f_216 func;
    f_45["_get_cb"] --> f_250["normalize"]
    class f_45 func;
    class f_250 func;
    f_45["_get_cb"] --> f_217["has_id"]
    class f_45 func;
    class f_217 func;
    f_45["_get_cb"] -->|b| f_115["angle"]
    class f_45 func;
    class f_115 func;
    f_45["_get_cb"] --> f_200["get_resname"]
    class f_45 func;
    class f_200 func;
    f_45["_get_cb"] -->|r2| f_46["_get_gly_cb_vector"]
    class f_45 func;
    class f_46 func;
    f_253["pcb_vectors_pymol"] --> f_336["write"]
    class f_253 func;
    class f_336 func;
    f_253["pcb_vectors_pymol"] --> f_179["get_array"]
    class f_253 func;
    class f_179 func;
    f_9["__init__"] -->|r1| f_231["is_aa"]
    class f_9 func;
    class f_231 func;
    f_9["__init__"] -->|r2| f_231["is_aa"]
    class f_9 func;
    class f_231 func;
    f_9["__init__"] -->|filename| f_121["as_handle"]
    class f_9 func;
    class f_121 func;
    f_9["__init__"] -->|handle| f_95["_tokenize"]
    class f_9 func;
    class f_95 func;
    f_9["__init__"] --> f_316["startswith"]
    class f_9 func;
    class f_316 func;
    f_9["__init__"] --> f_239["lower"]
    class f_9 func;
    class f_239 func;
    f_95["_tokenize"] --> f_316["startswith"]
    class f_95 func;
    class f_316 func;
    f_95["_tokenize"] --> f_280["rstrip"]
    class f_95 func;
    class f_280 func;
    f_95["_tokenize"] -->|token_buffer| f_236["join"]
    class f_95 func;
    class f_236 func;
    f_95["_tokenize"] --> f_93["_splitline"]
    class f_95 func;
    class f_93 func;
    f_95["_tokenize"] --> f_318["strip"]
    class f_95 func;
    class f_318 func;
    f_210["get_structure"] -->|structure_id| f_24["_build_structure"]
    class f_210 func;
    class f_24 func;
    f_210["get_structure"] --> f_301["set_header"]
    class f_210 func;
    class f_301 func;
    f_210["get_structure"] --> f_47["_get_header"]
    class f_210 func;
    class f_47 func;
    f_97["_update_header_entry"] -->|key| f_173["get"]
    class f_97 func;
    class f_173 func;
    f_47["_get_header"] --> f_97["_update_header_entry"]
    class f_47 func;
    class f_97 func;
    f_24["_build_structure"] -->|structure_id| f_228["init_structure"]
    class f_24 func;
    class f_228 func;
    f_24["_build_structure"] --> f_227["init_seg"]
    class f_24 func;
    class f_227 func;
    f_24["_build_structure"] -->|i| f_304["set_line_counter"]
    class f_24 func;
    class f_304 func;
    f_24["_build_structure"] -->|current_model_id<br>current_serial_id| f_225["init_model"]
    class f_24 func;
    class f_225 func;
    f_24["_build_structure"] -->|current_model_id| f_225["init_model"]
    class f_24 func;
    class f_225 func;
    f_24["_build_structure"] -->|current_chain_id| f_223["init_chain"]
    class f_24 func;
    class f_223 func;
    f_24["_build_structure"] -->|resname<br>hetatm_flag<br>int_resseq<br>icode| f_226["init_residue"]
    class f_24 func;
    class f_226 func;
    f_24["_build_structure"] --> f_329["upper"]
    class f_24 func;
    class f_329 func;
    f_24["_build_structure"] -->|name<br>coord<br>tempfactor<br>occupancy<br>altloc<br>name| f_221["init_atom"]
    class f_24 func;
    class f_221 func;
    f_24["_build_structure"] -->|anisou_array| f_295["set_anisou"]
    class f_24 func;
    class f_295 func;
    f_24["_build_structure"] -->|spacegroup<br>cell| f_312["set_symmetry"]
    class f_24 func;
    class f_312 func;
    f_210["get_structure"] -->|filename| f_121["as_handle"]
    class f_210 func;
    class f_121 func;
    f_210["get_structure"] -->|structure_id<br>handle| f_24["_build_structure"]
    class f_210 func;
    class f_24 func;
    f_24["_build_structure"] --> f_316["startswith"]
    class f_24 func;
    class f_316 func;
    f_24["_build_structure"] --> f_318["strip"]
    class f_24 func;
    class f_318 func;
    f_24["_build_structure"] -->|_records| f_243["map"]
    class f_24 func;
    class f_243 func;
    f_24["_build_structure"] -->|_anisors| f_243["map"]
    class f_24 func;
    class f_243 func;
    f_24["_build_structure"] --> f_326["update"]
    class f_24 func;
    class f_326 func;
    f_199["get_residues"] --> f_183["get_chains"]
    class f_199 func;
    class f_183 func;
    f_125["atom_to_internal_coordinates"] --> f_183["get_chains"]
    class f_125 func;
    class f_183 func;
    f_230["internal_to_atom_coordinates"] --> f_183["get_chains"]
    class f_230 func;
    class f_183 func;
    f_282["run_naccess"] -->|handle| f_139["close"]
    class f_282 func;
    class f_139 func;
    f_282["run_naccess"] -->|pdb_file<br>tmp_pdb_file| f_143["copy"]
    class f_282 func;
    class f_143 func;
    f_282["run_naccess"] --> f_311["set_structure"]
    class f_282 func;
    class f_311 func;
    f_282["run_naccess"] --> f_195["get_parent"]
    class f_282 func;
    class f_195 func;
    f_282["run_naccess"] -->|tmp_pdb_file| f_284["save"]
    class f_282 func;
    class f_284 func;
    f_282["run_naccess"] --> f_318["strip"]
    class f_282 func;
    class f_318 func;
    f_282["run_naccess"] --> f_270["readlines"]
    class f_282 func;
    class f_270 func;
    f_261["process_rsa_data"] --> f_316["startswith"]
    class f_261 func;
    class f_316 func;
    f_260["process_asa_data"] --> f_318["strip"]
    class f_260 func;
    class f_318 func;
    f_9["__init__"] -->|model<br>pdb_file| f_282["run_naccess"]
    class f_9 func;
    class f_282 func;
    f_9["__init__"] -->|res_data| f_261["process_rsa_data"]
    class f_9 func;
    class f_261 func;
    f_9["__init__"] -->|atm_data| f_260["process_asa_data"]
    class f_9 func;
    class f_260 func;
    f_1["Main_Script"] --> f_210["get_structure"]
    class f_1 main;
    class f_210 func;
    f_9["__init__"] --> f_184["get_coord"]
    class f_9 func;
    class f_184 func;
    f_54["_get_unique_parent_pairs"] --> f_195["get_parent"]
    class f_54 func;
    class f_195 func;
    f_54["_get_unique_parent_pairs"] -->|parent_pair_list| f_325["uniqueify"]
    class f_54 func;
    class f_325 func;
    f_285["search"] -->|atom_list<br>level| f_324["unfold_entities"]
    class f_285 func;
    class f_324 func;
    f_286["search_all"] -->|next_level_pair_list| f_54["_get_unique_parent_pairs"]
    class f_286 func;
    class f_54 func;
    f_311["set_structure"] --> f_228["init_structure"]
    class f_311 func;
    class f_228 func;
    f_311["set_structure"] --> f_227["init_seg"]
    class f_311 func;
    class f_227 func;
    f_311["set_structure"] --> f_108["add"]
    class f_311 func;
    class f_108 func;
    f_311["set_structure"] --> f_143["copy"]
    class f_311 func;
    class f_143 func;
    f_311["set_structure"] --> f_225["init_model"]
    class f_311 func;
    class f_225 func;
    f_311["set_structure"] -->|chain_id| f_223["init_chain"]
    class f_311 func;
    class f_223 func;
    f_311["set_structure"] --> f_226["init_residue"]
    class f_311 func;
    class f_226 func;
    f_41["_get_atom_line"] --> f_329["upper"]
    class f_41 func;
    class f_329 func;
    f_41["_get_atom_line"] --> f_318["strip"]
    class f_41 func;
    class f_318 func;
    f_85["_revert_write"] --> f_139["close"]
    class f_85 func;
    class f_139 func;
    f_85["_revert_write"] -->|truncate_to| f_287["seek"]
    class f_85 func;
    class f_287 func;
    f_85["_revert_write"] --> f_323["truncate"]
    class f_85 func;
    class f_323 func;
    f_284["save"] --> f_321["tell"]
    class f_284 func;
    class f_321 func;
    f_284["save"] --> f_191["get_list"]
    class f_284 func;
    class f_191 func;
    f_284["save"] -->|model| f_106["accept_model"]
    class f_284 func;
    class f_106 func;
    f_284["save"] --> f_336["write"]
    class f_284 func;
    class f_336 func;
    f_284["save"] -->|chain| f_105["accept_chain"]
    class f_284 func;
    class f_105 func;
    f_284["save"] -->|fhandle| f_85["_revert_write"]
    class f_284 func;
    class f_85 func;
    f_284["save"] --> f_215["get_unpacked_list"]
    class f_284 func;
    class f_215 func;
    f_284["save"] -->|residue| f_107["accept_residue"]
    class f_284 func;
    class f_107 func;
    f_284["save"] -->|atom| f_104["accept_atom"]
    class f_284 func;
    class f_104 func;
    f_284["save"] -->|s| f_336["write"]
    class f_284 func;
    class f_336 func;
    f_284["save"] --> f_139["close"]
    class f_284 func;
    class f_139 func;
    f_9["__init__"] --> f_236["join"]
    class f_9 func;
    class f_236 func;
    f_77["_print_default_format_warning"] --> f_336["write"]
    class f_77 func;
    class f_336 func;
    f_208["get_status_list"] -->|url| f_331["urlopen"]
    class f_208 func;
    class f_331 func;
    f_208["get_status_list"] --> f_318["strip"]
    class f_208 func;
    class f_318 func;
    f_198["get_recent_changes"] --> f_208["get_status_list"]
    class f_198 func;
    class f_208 func;
    f_175["get_all_entries"] -->|url| f_331["urlopen"]
    class f_175 func;
    class f_331 func;
    f_175["get_all_entries"] --> f_270["readlines"]
    class f_175 func;
    class f_270 func;
    f_176["get_all_obsolete"] -->|url| f_331["urlopen"]
    class f_176 func;
    class f_331 func;
    f_176["get_all_obsolete"] --> f_316["startswith"]
    class f_176 func;
    class f_316 func;
    f_176["get_all_obsolete"] --> f_314["split"]
    class f_176 func;
    class f_314 func;
    f_276["retrieve_pdb_file"] -->|file_format| f_77["_print_default_format_warning"]
    class f_276 func;
    class f_77 func;
    f_276["retrieve_pdb_file"] --> f_239["lower"]
    class f_276 func;
    class f_239 func;
    f_276["retrieve_pdb_file"] -->|archive_dict| f_236["join"]
    class f_276 func;
    class f_236 func;
    f_276["retrieve_pdb_file"] -->|path| f_236["join"]
    class f_276 func;
    class f_236 func;
    f_276["retrieve_pdb_file"] -->|path<br>archive| f_236["join"]
    class f_276 func;
    class f_236 func;
    f_276["retrieve_pdb_file"] --> f_330["urlcleanup"]
    class f_276 func;
    class f_330 func;
    f_276["retrieve_pdb_file"] -->|url<br>filename| f_332["urlretrieve"]
    class f_276 func;
    class f_332 func;
    f_276["retrieve_pdb_file"] -->|gz| f_340["writelines"]
    class f_276 func;
    class f_340 func;
    f_328["update_pdb"] -->|file_format| f_77["_print_default_format_warning"]
    class f_328 func;
    class f_77 func;
    f_328["update_pdb"] --> f_198["get_recent_changes"]
    class f_328 func;
    class f_198 func;
    f_328["update_pdb"] -->|pdb_code| f_276["retrieve_pdb_file"]
    class f_328 func;
    class f_276 func;
    f_328["update_pdb"] --> f_174["get_all_assemblies"]
    class f_328 func;
    class f_174 func;
    f_328["update_pdb"] -->|pdb_code<br>assembly_num| f_275["retrieve_assembly_file"]
    class f_328 func;
    class f_275 func;
    f_328["update_pdb"] --> f_236["join"]
    class f_328 func;
    class f_236 func;
    f_328["update_pdb"] -->|new_dir| f_236["join"]
    class f_328 func;
    class f_236 func;
    f_162["download_pdb_files"] -->|file_format| f_77["_print_default_format_warning"]
    class f_162 func;
    class f_77 func;
    f_162["download_pdb_files"] -->|pdb_codes| f_243["map"]
    class f_162 func;
    class f_243 func;
    f_174["get_all_assemblies"] -->|request| f_331["urlopen"]
    class f_174 func;
    class f_331 func;
    f_174["get_all_assemblies"] --> f_267["read"]
    class f_174 func;
    class f_267 func;
    f_322["transform"] --> f_314["split"]
    class f_322 func;
    class f_314 func;
    f_322["transform"] --> f_239["lower"]
    class f_322 func;
    class f_239 func;
    f_174["get_all_assemblies"] -->|transform<br>assemblies| f_243["map"]
    class f_174 func;
    class f_243 func;
    f_275["retrieve_assembly_file"] --> f_239["lower"]
    class f_275 func;
    class f_239 func;
    f_275["retrieve_assembly_file"] -->|file_format| f_77["_print_default_format_warning"]
    class f_275 func;
    class f_77 func;
    f_275["retrieve_assembly_file"] -->|path| f_236["join"]
    class f_275 func;
    class f_236 func;
    f_275["retrieve_assembly_file"] -->|path<br>archive_fn| f_236["join"]
    class f_275 func;
    class f_236 func;
    f_275["retrieve_assembly_file"] --> f_330["urlcleanup"]
    class f_275 func;
    class f_330 func;
    f_275["retrieve_assembly_file"] -->|url<br>assembly_gz_file| f_332["urlretrieve"]
    class f_275 func;
    class f_332 func;
    f_275["retrieve_assembly_file"] -->|gz| f_340["writelines"]
    class f_275 func;
    class f_340 func;
    f_158["download_all_assemblies"] -->|file_format| f_77["_print_default_format_warning"]
    class f_158 func;
    class f_77 func;
    f_158["download_all_assemblies"] --> f_174["get_all_assemblies"]
    class f_158 func;
    class f_174 func;
    f_158["download_all_assemblies"] -->|pdb_code<br>assembly_num| f_320["submit"]
    class f_158 func;
    class f_320 func;
    f_158["download_all_assemblies"] --> f_340["writelines"]
    class f_158 func;
    class f_340 func;
    f_160["download_entire_pdb"] -->|file_format| f_77["_print_default_format_warning"]
    class f_160 func;
    class f_77 func;
    f_160["download_entire_pdb"] --> f_175["get_all_entries"]
    class f_160 func;
    class f_175 func;
    f_160["download_entire_pdb"] -->|entries| f_243["map"]
    class f_160 func;
    class f_243 func;
    f_160["download_entire_pdb"] --> f_340["writelines"]
    class f_160 func;
    class f_340 func;
    f_161["download_obsolete_entries"] -->|file_format| f_77["_print_default_format_warning"]
    class f_161 func;
    class f_77 func;
    f_161["download_obsolete_entries"] --> f_176["get_all_obsolete"]
    class f_161 func;
    class f_176 func;
    f_161["download_obsolete_entries"] -->|entries| f_243["map"]
    class f_161 func;
    class f_243 func;
    f_161["download_obsolete_entries"] --> f_340["writelines"]
    class f_161 func;
    class f_340 func;
    f_204["get_seqres_file"] -->|url<br>savefile| f_332["urlretrieve"]
    class f_204 func;
    class f_332 func;
    f_1["Main_Script"] --> f_328["update_pdb"]
    class f_1 main;
    class f_328 func;
    f_1["Main_Script"] --> f_160["download_entire_pdb"]
    class f_1 main;
    class f_160 func;
    f_1["Main_Script"] --> f_158["download_all_assemblies"]
    class f_1 main;
    class f_158 func;
    f_1["Main_Script"] -->|pdb_path| f_161["download_obsolete_entries"]
    class f_1 main;
    class f_161 func;
    f_1["Main_Script"] -->|pdb_code| f_276["retrieve_pdb_file"]
    class f_1 main;
    class f_276 func;
    f_1["Main_Script"] --> f_174["get_all_assemblies"]
    class f_1 main;
    class f_174 func;
    f_1["Main_Script"] -->|pdb_code<br>assembly_num| f_275["retrieve_assembly_file"]
    class f_1 main;
    class f_275 func;
    f_1["Main_Script"] -->|pdb_id| f_276["retrieve_pdb_file"]
    class f_1 main;
    class f_276 func;
    f_1["Main_Script"] -->|pdb_id<br>assembly_num| f_275["retrieve_assembly_file"]
    class f_1 main;
    class f_275 func;
    f_74["_parse_resolution_from"] -->|candidate<br>namespaces| f_167["find"]
    class f_74 func;
    class f_167 func;
    f_70["_parse_header_from"] -->|namespaces| f_167["find"]
    class f_70 func;
    class f_167 func;
    f_70["_parse_header_from"] -->|tree<br>namespaces| f_74["_parse_resolution_from"]
    class f_70 func;
    class f_74 func;
    f_68["_parse_atom_from"] -->|namespaces| f_167["find"]
    class f_68 func;
    class f_167 func;
    f_73["_parse_residue_id_from"] -->|namespaces| f_167["find"]
    class f_73 func;
    class f_167 func;
    f_210["get_structure"] --> f_287["seek"]
    class f_210 func;
    class f_287 func;
    f_210["get_structure"] -->|tree<br>namespaces| f_70["_parse_header_from"]
    class f_210 func;
    class f_70 func;
    f_210["get_structure"] --> f_228["init_structure"]
    class f_210 func;
    class f_228 func;
    f_210["get_structure"] -->|header| f_301["set_header"]
    class f_210 func;
    class f_301 func;
    f_210["get_structure"] -->|namespaces| f_167["find"]
    class f_210 func;
    class f_167 func;
    f_210["get_structure"] -->|element<br>namespaces| f_73["_parse_residue_id_from"]
    class f_210 func;
    class f_73 func;
    f_210["get_structure"] -->|builder_model_count<br>model_number| f_225["init_model"]
    class f_210 func;
    class f_225 func;
    f_210["get_structure"] -->|chain_id| f_223["init_chain"]
    class f_210 func;
    class f_223 func;
    f_210["get_structure"] -->|component_id| f_226["init_residue"]
    class f_210 func;
    class f_226 func;
    f_210["get_structure"] --> f_221["init_atom"]
    class f_210 func;
    class f_221 func;
    f_210["get_structure"] -->|element<br>namespaces| f_68["_parse_atom_from"]
    class f_210 func;
    class f_68 func;
    f_210["get_structure"] -->|id| f_228["init_structure"]
    class f_210 func;
    class f_228 func;
    f_210["get_structure"] -->|file| f_121["as_handle"]
    class f_210 func;
    class f_121 func;
    f_210["get_structure"] --> f_270["readlines"]
    class f_210 func;
    class f_270 func;
    f_210["get_structure"] -->|lines| f_67["_parse"]
    class f_210 func;
    class f_67 func;
    f_67["_parse"] -->|header_coords_trailer| f_47["_get_header"]
    class f_67 func;
    class f_47 func;
    f_67["_parse"] -->|coords_trailer| f_69["_parse_coordinates"]
    class f_67 func;
    class f_69 func;
    f_47["_get_header"] --> f_304["set_line_counter"]
    class f_47 func;
    class f_304 func;
    f_47["_get_header"] -->|header| f_71["_parse_pdb_header_list"]
    class f_47 func;
    class f_71 func;
    f_69["_parse_coordinates"] --> f_280["rstrip"]
    class f_69 func;
    class f_280 func;
    f_69["_parse_coordinates"] -->|global_line_counter| f_304["set_line_counter"]
    class f_69 func;
    class f_304 func;
    f_69["_parse_coordinates"] --> f_318["strip"]
    class f_69 func;
    class f_318 func;
    f_69["_parse_coordinates"] -->|current_model_id| f_225["init_model"]
    class f_69 func;
    class f_225 func;
    f_69["_parse_coordinates"] --> f_314["split"]
    class f_69 func;
    class f_314 func;
    f_69["_parse_coordinates"] -->|global_line_counter| f_55["_handle_PDB_exception"]
    class f_69 func;
    class f_55 func;
    f_69["_parse_coordinates"] -->|message<br>global_line_counter| f_55["_handle_PDB_exception"]
    class f_69 func;
    class f_55 func;
    f_69["_parse_coordinates"] --> f_329["upper"]
    class f_69 func;
    class f_329 func;
    f_69["_parse_coordinates"] -->|current_segid| f_227["init_seg"]
    class f_69 func;
    class f_227 func;
    f_69["_parse_coordinates"] -->|current_chain_id| f_223["init_chain"]
    class f_69 func;
    class f_223 func;
    f_69["_parse_coordinates"] -->|resname<br>hetero_flag<br>resseq<br>icode| f_226["init_residue"]
    class f_69 func;
    class f_226 func;
    f_69["_parse_coordinates"] -->|name<br>coord<br>bfactor<br>occupancy<br>altloc<br>fullname<br>serial_number<br>element| f_221["init_atom"]
    class f_69 func;
    class f_221 func;
    f_69["_parse_coordinates"] -->|name<br>coord<br>pqr_charge<br>radius<br>altloc<br>fullname<br>serial_number<br>element<br>pqr_charge<br>radius| f_221["init_atom"]
    class f_69 func;
    class f_221 func;
    f_69["_parse_coordinates"] -->|anisou_array| f_295["set_anisou"]
    class f_69 func;
    class f_295 func;
    f_69["_parse_coordinates"] -->|current_model_id<br>serial_num| f_225["init_model"]
    class f_69 func;
    class f_225 func;
    f_69["_parse_coordinates"] -->|siguij_array| f_310["set_siguij"]
    class f_69 func;
    class f_310 func;
    f_69["_parse_coordinates"] -->|sigatm_array| f_309["set_sigatm"]
    class f_69 func;
    class f_309 func;
    f_268["read_PIC"] --> f_71["_parse_pdb_header_list"]
    class f_268 func;
    class f_71 func;
    f_268["read_PIC"] --> f_288["set"]
    class f_268 func;
    class f_288 func;
    f_146["default_hedron"] -->|rhcl| f_236["join"]
    class f_146 func;
    class f_236 func;
    f_145["default_dihedron"] --> f_108["add"]
    class f_145 func;
    class f_108 func;
    f_112["ake_recurse"] -->|ak| f_229["insert"]
    class f_112 func;
    class f_229 func;
    f_111["ak_expand"] --> f_315["split_akl"]
    class f_111 func;
    class f_315 func;
    f_110["ak_add"] -->|ak| f_108["add"]
    class f_110 func;
    class f_108 func;
    f_168["finish_chain"] -->|shl12<br>sha<br>shl23<br>sda<br>bfacs| f_56["_hedraDict2chain"]
    class f_168 func;
    class f_56 func;
    f_268["read_PIC"] -->|file| f_121["as_handle"]
    class f_268 func;
    class f_121 func;
    f_268["read_PIC"] --> f_270["readlines"]
    class f_268 func;
    class f_270 func;
    f_268["read_PIC"] --> f_316["startswith"]
    class f_268 func;
    class f_316 func;
    f_268["read_PIC"] --> f_318["strip"]
    class f_268 func;
    class f_318 func;
    f_268["read_PIC"] -->|header_dict| f_301["set_header"]
    class f_268 func;
    class f_301 func;
    f_268["read_PIC"] --> f_226["init_residue"]
    class f_268 func;
    class f_226 func;
    f_268["read_PIC"] --> f_233["is_disordered"]
    class f_268 func;
    class f_233 func;
    f_268["read_PIC"] --> f_333["values"]
    class f_268 func;
    class f_333 func;
    f_268["read_PIC"] -->|ak| f_108["add"]
    class f_268 func;
    class f_108 func;
    f_268["read_PIC"] -->|coord| f_221["init_atom"]
    class f_268 func;
    class f_221 func;
    f_268["read_PIC"] --> f_210["get_structure"]
    class f_268 func;
    class f_210 func;
    f_269["read_PIC_seq"] --> f_272["replace"]
    class f_269 func;
    class f_272 func;
    f_269["read_PIC_seq"] --> f_314["split"]
    class f_269 func;
    class f_314 func;
    f_269["read_PIC_seq"] --> f_329["upper"]
    class f_269 func;
    class f_329 func;
    f_269["read_PIC_seq"] -->|output| f_336["write"]
    class f_269 func;
    class f_336 func;
    f_269["read_PIC_seq"] --> f_287["seek"]
    class f_269 func;
    class f_287 func;
    f_269["read_PIC_seq"] -->|sp| f_268["read_PIC"]
    class f_269 func;
    class f_268 func;
    f_98["_wpr"] --> f_173["get"]
    class f_98 func;
    class f_173 func;
    f_98["_wpr"] --> f_336["write"]
    class f_98 func;
    class f_336 func;
    f_98["_wpr"] -->|pdbid<br>chainid| f_100["_write_PIC"]
    class f_98 func;
    class f_100 func;
    f_98["_wpr"] -->|entity| f_84["_residue_string"]
    class f_98 func;
    class f_84 func;
    f_34["_enumerate_entity_atoms"] --> f_180["get_atoms"]
    class f_34 func;
    class f_180 func;
    f_34["_enumerate_entity_atoms"] --> f_206["get_serial_number"]
    class f_34 func;
    class f_206 func;
    f_34["_enumerate_entity_atoms"] --> f_199["get_residues"]
    class f_34 func;
    class f_199 func;
    f_34["_enumerate_entity_atoms"] --> f_233["is_disordered"]
    class f_34 func;
    class f_233 func;
    f_34["_enumerate_entity_atoms"] --> f_333["values"]
    class f_34 func;
    class f_333 func;
    f_34["_enumerate_entity_atoms"] --> f_215["get_unpacked_list"]
    class f_34 func;
    class f_215 func;
    f_34["_enumerate_entity_atoms"] -->|anum| f_308["set_serial_number"]
    class f_34 func;
    class f_308 func;
    f_165["enumerate_atoms"] --> f_195["get_parent"]
    class f_165 func;
    class f_195 func;
    f_165["enumerate_atoms"] -->|mdl| f_34["_enumerate_entity_atoms"]
    class f_165 func;
    class f_34 func;
    f_165["enumerate_atoms"] -->|entity| f_34["_enumerate_entity_atoms"]
    class f_165 func;
    class f_34 func;
    f_338["write_PIC"] -->|entity| f_165["enumerate_atoms"]
    class f_338 func;
    class f_165 func;
    f_338["write_PIC"] -->|file| f_121["as_handle"]
    class f_338 func;
    class f_121 func;
    f_338["write_PIC"] --> f_233["is_disordered"]
    class f_338 func;
    class f_233 func;
    f_338["write_PIC"] --> f_333["values"]
    class f_338 func;
    class f_333 func;
    f_338["write_PIC"] -->|r<br>fp<br>pdbid<br>chainid| f_98["_wpr"]
    class f_338 func;
    class f_98 func;
    f_338["write_PIC"] -->|entity<br>fp<br>pdbid<br>chainid| f_98["_wpr"]
    class f_338 func;
    class f_98 func;
    f_338["write_PIC"] --> f_173["get"]
    class f_338 func;
    class f_173 func;
    f_338["write_PIC"] --> f_254["pdb_date"]
    class f_338 func;
    class f_254 func;
    f_338["write_PIC"] --> f_336["write"]
    class f_338 func;
    class f_336 func;
    f_338["write_PIC"] --> f_329["upper"]
    class f_338 func;
    class f_329 func;
    f_283["run_psea"] --> f_314["split"]
    class f_283 func;
    class f_314 func;
    f_283["run_psea"] -->|cmd| f_281["run"]
    class f_283 func;
    class f_281 func;
    f_283["run_psea"] --> f_318["strip"]
    class f_283 func;
    class f_318 func;
    f_263["psea"] -->|pname| f_283["run_psea"]
    class f_263 func;
    class f_283 func;
    f_118["annotate"] --> f_191["get_list"]
    class f_118 func;
    class f_191 func;
    f_118["annotate"] -->|res| f_231["is_aa"]
    class f_118 func;
    class f_231 func;
    f_9["__init__"] -->|filename| f_263["psea"]
    class f_9 func;
    class f_263 func;
    f_9["__init__"] -->|ss_seq| f_264["psea2HEC"]
    class f_9 func;
    class f_264 func;
    f_9["__init__"] -->|model<br>ss_seq| f_118["annotate"]
    class f_9 func;
    class f_118 func;
    f_1["Main_Script"] --> f_235["items"]
    class f_1 main;
    class f_235 func;
    f_231["is_aa"] --> f_200["get_resname"]
    class f_231 func;
    class f_200 func;
    f_231["is_aa"] --> f_329["upper"]
    class f_231 func;
    class f_329 func;
    f_234["is_nucleic"] --> f_200["get_resname"]
    class f_234 func;
    class f_200 func;
    f_234["is_nucleic"] --> f_329["upper"]
    class f_234 func;
    class f_329 func;
    f_196["get_phi_psi_list"] --> f_216["get_vector"]
    class f_196 func;
    class f_216 func;
    f_196["get_phi_psi_list"] -->|cp<br>n<br>ca<br>c| f_134["calc_dihedral"]
    class f_196 func;
    class f_134 func;
    f_196["get_phi_psi_list"] -->|n<br>ca<br>c<br>nn| f_134["calc_dihedral"]
    class f_196 func;
    class f_134 func;
    f_212["get_tau_list"] --> f_182["get_ca_list"]
    class f_212 func;
    class f_182 func;
    f_212["get_tau_list"] --> f_216["get_vector"]
    class f_212 func;
    class f_216 func;
    f_212["get_tau_list"] -->|v1<br>v2<br>v3<br>v4| f_134["calc_dihedral"]
    class f_212 func;
    class f_134 func;
    f_212["get_tau_list"] --> f_195["get_parent"]
    class f_212 func;
    class f_195 func;
    f_213["get_theta_list"] --> f_182["get_ca_list"]
    class f_213 func;
    class f_182 func;
    f_213["get_theta_list"] --> f_216["get_vector"]
    class f_213 func;
    class f_216 func;
    f_213["get_theta_list"] -->|v1<br>v2<br>v3| f_133["calc_angle"]
    class f_213 func;
    class f_133 func;
    f_213["get_theta_list"] --> f_195["get_parent"]
    class f_213 func;
    class f_195 func;
    f_205["get_sequence"] --> f_236["join"]
    class f_205 func;
    class f_236 func;
    f_205["get_sequence"] --> f_173["get"]
    class f_205 func;
    class f_173 func;
    f_205["get_sequence"] --> f_200["get_resname"]
    class f_205 func;
    class f_200 func;
    f_16["_accept"] -->|residue| f_231["is_aa"]
    class f_16 func;
    class f_231 func;
    f_16["_accept"] --> f_200["get_resname"]
    class f_16 func;
    class f_200 func;
    f_131["build_peptides"] --> f_190["get_level"]
    class f_131 func;
    class f_190 func;
    f_131["build_peptides"] --> f_191["get_list"]
    class f_131 func;
    class f_191 func;
    f_59["_is_connected"] --> f_217["has_id"]
    class f_59 func;
    class f_217 func;
    f_59["_is_connected"] --> f_233["is_disordered"]
    class f_59 func;
    class f_233 func;
    f_59["_is_connected"] --> f_153["disordered_get_list"]
    class f_59 func;
    class f_153 func;
    f_59["_is_connected"] --> f_177["get_altloc"]
    class f_59 func;
    class f_177 func;
    f_59["_is_connected"] -->|c_altloc| f_156["disordered_select"]
    class f_59 func;
    class f_156 func;
    f_59["_is_connected"] -->|n_altloc| f_156["disordered_select"]
    class f_59 func;
    class f_156 func;
    f_14["__repr__"] --> f_200["get_resname"]
    class f_14 func;
    class f_200 func;
    f_108["add"] -->|atom_id| f_217["has_id"]
    class f_108 func;
    class f_217 func;
    f_108["add"] --> f_151["disordered_get"]
    class f_108 func;
    class f_151 func;
    f_108["add"] --> f_233["is_disordered"]
    class f_108 func;
    class f_233 func;
    f_108["add"] --> f_200["get_resname"]
    class f_108 func;
    class f_200 func;
    f_313["sort"] --> f_153["disordered_get_list"]
    class f_313 func;
    class f_153 func;
    f_150["disordered_add"] --> f_200["get_resname"]
    class f_150 func;
    class f_200 func;
    f_150["disordered_add"] -->|chain| f_306["set_parent"]
    class f_150 func;
    class f_306 func;
    f_150["disordered_add"] -->|resname| f_154["disordered_has_id"]
    class f_150 func;
    class f_154 func;
    f_150["disordered_add"] -->|resname| f_156["disordered_select"]
    class f_150 func;
    class f_156 func;
    f_155["disordered_remove"] -->|child| f_156["disordered_select"]
    class f_155 func;
    class f_156 func;
    f_42["_get_atom_radius"] --> f_316["startswith"]
    class f_42 func;
    class f_316 func;
    f_42["_get_atom_radius"] --> f_164["endswith"]
    class f_42 func;
    class f_164 func;
    f_79["_read_vertex_array"] --> f_314["split"]
    class f_79 func;
    class f_314 func;
    f_211["get_surface"] -->|model| f_324["unfold_entities"]
    class f_211 func;
    class f_324 func;
    f_211["get_surface"] -->|atom| f_42["_get_atom_radius"]
    class f_211 func;
    class f_42 func;
    f_211["get_surface"] --> f_336["write"]
    class f_211 func;
    class f_336 func;
    f_211["get_surface"] -->|surface_file| f_79["_read_vertex_array"]
    class f_211 func;
    class f_79 func;
    f_274["residue_depth"] --> f_215["get_unpacked_list"]
    class f_274 func;
    class f_215 func;
    f_274["residue_depth"] --> f_184["get_coord"]
    class f_274 func;
    class f_184 func;
    f_274["residue_depth"] -->|coord<br>surface| f_244["min_dist"]
    class f_274 func;
    class f_244 func;
    f_132["ca_depth"] --> f_217["has_id"]
    class f_132 func;
    class f_217 func;
    f_132["ca_depth"] --> f_184["get_coord"]
    class f_132 func;
    class f_184 func;
    f_132["ca_depth"] -->|coord<br>surface| f_244["min_dist"]
    class f_132 func;
    class f_244 func;
    f_9["__init__"] -->|model| f_324["unfold_entities"]
    class f_9 func;
    class f_324 func;
    f_9["__init__"] -->|model| f_211["get_surface"]
    class f_9 func;
    class f_211 func;
    f_9["__init__"] -->|residue| f_231["is_aa"]
    class f_9 func;
    class f_231 func;
    f_9["__init__"] -->|residue<br>surface| f_274["residue_depth"]
    class f_9 func;
    class f_274 func;
    f_9["__init__"] -->|residue<br>surface| f_132["ca_depth"]
    class f_9 func;
    class f_132 func;
    f_1["Main_Script"] --> f_326["update"]
    class f_1 main;
    class f_326 func;
    f_9["__init__"] --> f_143["copy"]
    class f_9 func;
    class f_143 func;
    f_9["__init__"] -->|radii_dict| f_326["update"]
    class f_9 func;
    class f_326 func;
    f_9["__init__"] --> f_30["_compute_sphere"]
    class f_9 func;
    class f_30 func;
    f_141["compute"] --> f_180["get_atoms"]
    class f_141 func;
    class f_180 func;
    f_141["compute"] --> f_288["set"]
    class f_141 func;
    class f_288 func;
    f_141["compute"] --> f_143["copy"]
    class f_141 func;
    class f_143 func;
    f_141["compute"] -->|twice_maxradii| f_285["search"]
    class f_141 func;
    class f_285 func;
    f_141["compute"] --> f_285["search"]
    class f_141 func;
    class f_285 func;
    f_141["compute"] -->|atoms| f_288["set"]
    class f_141 func;
    class f_288 func;
    f_88["_scale_residue"] -->|scaleMtx| f_120["applyMtx"]
    class f_88 func;
    class f_120 func;
    f_339["write_SCAD"] --> f_183["get_chains"]
    class f_339 func;
    class f_183 func;
    f_339["write_SCAD"] --> f_125["atom_to_internal_coordinates"]
    class f_339 func;
    class f_125 func;
    f_339["write_SCAD"] --> f_230["internal_to_atom_coordinates"]
    class f_339 func;
    class f_230 func;
    f_339["write_SCAD"] -->|scale| f_218["homog_scale_mtx"]
    class f_339 func;
    class f_218 func;
    f_339["write_SCAD"] -->|file| f_121["as_handle"]
    class f_339 func;
    class f_121 func;
    f_339["write_SCAD"] -->|peptide_scad| f_336["write"]
    class f_339 func;
    class f_336 func;
    f_339["write_SCAD"] --> f_173["get"]
    class f_339 func;
    class f_173 func;
    f_339["write_SCAD"] --> f_336["write"]
    class f_339 func;
    class f_336 func;
    f_339["write_SCAD"] -->|fp| f_101["_write_SCAD"]
    class f_339 func;
    class f_101 func;
    f_325["uniqueify"] -->|items| f_288["set"]
    class f_325 func;
    class f_288 func;
    f_214["get_unique_parents"] --> f_195["get_parent"]
    class f_214 func;
    class f_195 func;
    f_324["unfold_entities"] --> f_190["get_level"]
    class f_324 func;
    class f_190 func;
    f_324["unfold_entities"] -->|target_level| f_220["index"]
    class f_324 func;
    class f_220 func;
    f_324["unfold_entities"] -->|level| f_220["index"]
    class f_324 func;
    class f_220 func;
    f_324["unfold_entities"] --> f_195["get_parent"]
    class f_324 func;
    class f_195 func;
    f_183["get_chains"] --> f_192["get_models"]
    class f_183 func;
    class f_192 func;
    f_9["__init__"] -->|m1| f_324["unfold_entities"]
    class f_9 func;
    class f_324 func;
    f_9["__init__"] -->|m2| f_324["unfold_entities"]
    class f_9 func;
    class f_324 func;
    f_9["__init__"] -->|r1<br>aa1| f_94["_test_equivalence"]
    class f_9 func;
    class f_94 func;
    f_9["__init__"] -->|r2<br>aa2| f_94["_test_equivalence"]
    class f_9 func;
    class f_94 func;
    f_94["_test_equivalence"] --> f_200["get_resname"]
    class f_94 func;
    class f_200 func;
    f_58["_is_completely_disordered"] --> f_215["get_unpacked_list"]
    class f_58 func;
    class f_215 func;
    f_58["_is_completely_disordered"] --> f_177["get_altloc"]
    class f_58 func;
    class f_177 func;
    f_225["init_model"] --> f_108["add"]
    class f_225 func;
    class f_108 func;
    f_223["init_chain"] -->|chain_id| f_217["has_id"]
    class f_223 func;
    class f_217 func;
    f_223["init_chain"] --> f_108["add"]
    class f_223 func;
    class f_108 func;
    f_226["init_residue"] -->|res_id| f_217["has_id"]
    class f_226 func;
    class f_217 func;
    f_226["init_residue"] --> f_233["is_disordered"]
    class f_226 func;
    class f_233 func;
    f_226["init_residue"] -->|resname| f_154["disordered_has_id"]
    class f_226 func;
    class f_154 func;
    f_226["init_residue"] -->|resname| f_156["disordered_select"]
    class f_226 func;
    class f_156 func;
    f_226["init_residue"] -->|new_residue| f_150["disordered_add"]
    class f_226 func;
    class f_150 func;
    f_226["init_residue"] -->|duplicate_residue| f_58["_is_completely_disordered"]
    class f_226 func;
    class f_58 func;
    f_226["init_residue"] -->|res_id| f_147["detach_child"]
    class f_226 func;
    class f_147 func;
    f_226["init_residue"] -->|disordered_residue| f_108["add"]
    class f_226 func;
    class f_108 func;
    f_226["init_residue"] -->|duplicate_residue| f_150["disordered_add"]
    class f_226 func;
    class f_150 func;
    f_226["init_residue"] --> f_108["add"]
    class f_226 func;
    class f_108 func;
    f_221["init_atom"] -->|name| f_217["has_id"]
    class f_221 func;
    class f_217 func;
    f_221["init_atom"] --> f_186["get_fullname"]
    class f_221 func;
    class f_186 func;
    f_221["init_atom"] --> f_233["is_disordered"]
    class f_221 func;
    class f_233 func;
    f_221["init_atom"] --> f_150["disordered_add"]
    class f_221 func;
    class f_150 func;
    f_221["init_atom"] -->|name| f_147["detach_child"]
    class f_221 func;
    class f_147 func;
    f_221["init_atom"] -->|disordered_atom| f_108["add"]
    class f_221 func;
    class f_108 func;
    f_221["init_atom"] -->|duplicate_atom| f_150["disordered_add"]
    class f_221 func;
    class f_150 func;
    f_221["init_atom"] --> f_170["flag_disordered"]
    class f_221 func;
    class f_170 func;
    f_221["init_atom"] --> f_108["add"]
    class f_221 func;
    class f_108 func;
    f_296["set_atoms"] --> f_184["get_coord"]
    class f_296 func;
    class f_184 func;
    f_296["set_atoms"] -->|fixed_coord<br>moving_coord| f_288["set"]
    class f_296 func;
    class f_288 func;
    f_296["set_atoms"] --> f_281["run"]
    class f_296 func;
    class f_281 func;
    f_296["set_atoms"] --> f_201["get_rms"]
    class f_296 func;
    class f_201 func;
    f_296["set_atoms"] --> f_202["get_rotran"]
    class f_296 func;
    class f_202 func;
    f_119["apply"] -->|rot<br>tran| f_322["transform"]
    class f_119 func;
    class f_322 func;
    f_197["get_predictions"] -->|url| f_331["urlopen"]
    class f_197 func;
    class f_331 func;
    f_197["get_predictions"] --> f_267["read"]
    class f_197 func;
    class f_267 func;
    f_51["_get_mmcif_file_path_for"] --> f_314["split"]
    class f_51 func;
    class f_314 func;
    f_51["_get_mmcif_file_path_for"] -->|directory<br>file_name| f_236["join"]
    class f_51 func;
    class f_236 func;
    f_159["download_cif_for"] -->|prediction<br>directory| f_51["_get_mmcif_file_path_for"]
    class f_159 func;
    class f_51 func;
    f_159["download_cif_for"] -->|cif_url| f_331["urlopen"]
    class f_159 func;
    class f_331 func;
    f_159["download_cif_for"] --> f_267["read"]
    class f_159 func;
    class f_267 func;
    f_159["download_cif_for"] -->|data| f_336["write"]
    class f_159 func;
    class f_336 func;
    f_209["get_structural_models_for"] -->|qualifier| f_197["get_predictions"]
    class f_209 func;
    class f_197 func;
    f_209["get_structural_models_for"] -->|prediction<br>directory| f_51["_get_mmcif_file_path_for"]
    class f_209 func;
    class f_51 func;
    f_209["get_structural_models_for"] -->|prediction<br>directory| f_159["download_cif_for"]
    class f_209 func;
    class f_159 func;
    f_209["get_structural_models_for"] -->|qualifier<br>mmcif_path| f_210["get_structure"]
    class f_209 func;
    class f_210 func;
    f_187["get_guide_coord_from_structure"] --> f_183["get_chains"]
    class f_187 func;
    class f_183 func;
    f_307["set_reference"] -->|structure| f_187["get_guide_coord_from_structure"]
    class f_307 func;
    class f_187 func;
    f_113["align"] -->|structure| f_187["get_guide_coord_from_structure"]
    class f_113 func;
    class f_187 func;
    f_113["align"] -->|coordsA<br>coordsB| f_288["set"]
    class f_113 func;
    class f_288 func;
    f_113["align"] --> f_281["run"]
    class f_113 func;
    class f_281 func;
    f_113["align"] -->|best_alignment| f_66["_optimize"]
    class f_113 func;
    class f_66 func;
    f_113["align"] --> f_183["get_chains"]
    class f_113 func;
    class f_183 func;
    f_113["align"] --> f_215["get_unpacked_list"]
    class f_113 func;
    class f_215 func;
    f_113["align"] -->|rotmtx<br>trvec| f_322["transform"]
    class f_113 func;
    class f_322 func;
    f_66["_optimize"] -->|coordsA<br>coordsB| f_288["set"]
    class f_66 func;
    class f_288 func;
    f_66["_optimize"] --> f_281["run"]
    class f_66 func;
    class f_281 func;
    f_319["structure_rebuild_test"] -->|verbose| f_125["atom_to_internal_coordinates"]
    class f_319 func;
    class f_125 func;
    f_319["structure_rebuild_test"] -->|entity<br>sp| f_338["write_PIC"]
    class f_319 func;
    class f_338 func;
    f_319["structure_rebuild_test"] --> f_287["seek"]
    class f_319 func;
    class f_287 func;
    f_319["structure_rebuild_test"] -->|sp| f_268["read_PIC"]
    class f_319 func;
    class f_268 func;
    f_319["structure_rebuild_test"] --> f_183["get_chains"]
    class f_319 func;
    class f_183 func;
    f_319["structure_rebuild_test"] -->|pdb2| f_273["report_IC"]
    class f_319 func;
    class f_273 func;
    f_319["structure_rebuild_test"] -->|verbose| f_230["internal_to_atom_coordinates"]
    class f_319 func;
    class f_230 func;
    f_319["structure_rebuild_test"] -->|entity<br>pdb2| f_140["compare_residues"]
    class f_319 func;
    class f_140 func;
    f_273["report_IC"] --> f_173["get"]
    class f_273 func;
    class f_173 func;
    f_0["IC_duplicate"] --> f_199["get_residues"]
    class f_0 func;
    class f_199 func;
    f_0["IC_duplicate"] --> f_125["atom_to_internal_coordinates"]
    class f_0 func;
    class f_125 func;
    f_0["IC_duplicate"] -->|entity<br>sp| f_338["write_PIC"]
    class f_0 func;
    class f_338 func;
    f_0["IC_duplicate"] --> f_287["seek"]
    class f_0 func;
    class f_287 func;
    f_0["IC_duplicate"] -->|sp| f_268["read_PIC"]
    class f_0 func;
    class f_268 func;
    f_22["_atmfid_d2h"] --> f_185["get_full_id"]
    class f_22 func;
    class f_185 func;
    f_28["_cmp_atm"] --> f_185["get_full_id"]
    class f_28 func;
    class f_185 func;
    f_28["_cmp_atm"] -->|a0| f_22["_atmfid_d2h"]
    class f_28 func;
    class f_22 func;
    f_28["_cmp_atm"] --> f_184["get_coord"]
    class f_28 func;
    class f_184 func;
    f_29["_cmp_res"] -->|ak| f_173["get"]
    class f_29 func;
    class f_173 func;
    f_29["_cmp_res"] -->|aknd| f_173["get"]
    class f_29 func;
    class f_173 func;
    f_29["_cmp_res"] --> f_233["is_disordered"]
    class f_29 func;
    class f_233 func;
    f_29["_cmp_res"] -->|r0<br>r1<br>a0<br>a1<br>verbose<br>cmpdict| f_28["_cmp_atm"]
    class f_29 func;
    class f_28 func;
    f_29["_cmp_res"] -->|r0<br>r1<br>verbose<br>cmpdict| f_28["_cmp_atm"]
    class f_29 func;
    class f_28 func;
    f_29["_cmp_res"] -->|da0k| f_173["get"]
    class f_29 func;
    class f_173 func;
    f_29["_cmp_res"] --> f_185["get_full_id"]
    class f_29 func;
    class f_185 func;
    f_140["compare_residues"] --> f_185["get_full_id"]
    class f_140 func;
    class f_185 func;
    f_140["compare_residues"] --> f_183["get_chains"]
    class f_140 func;
    class f_183 func;
    f_140["compare_residues"] --> f_199["get_residues"]
    class f_140 func;
    class f_199 func;
    f_140["compare_residues"] --> f_233["is_disordered"]
    class f_140 func;
    class f_233 func;
    f_140["compare_residues"] --> f_333["values"]
    class f_140 func;
    class f_333 func;
    f_140["compare_residues"] -->|dr0<br>dr1<br>verbose<br>cmpdict| f_29["_cmp_res"]
    class f_140 func;
    class f_29 func;
    f_140["compare_residues"] -->|r0<br>r1<br>verbose<br>cmpdict| f_29["_cmp_res"]
    class f_140 func;
    class f_29 func;
    f_337["write_PDB"] -->|entity| f_165["enumerate_atoms"]
    class f_337 func;
    class f_165 func;
    f_337["write_PDB"] -->|file| f_121["as_handle"]
    class f_337 func;
    class f_121 func;
    f_337["write_PDB"] --> f_173["get"]
    class f_337 func;
    class f_173 func;
    f_337["write_PDB"] --> f_254["pdb_date"]
    class f_337 func;
    class f_254 func;
    f_337["write_PDB"] --> f_336["write"]
    class f_337 func;
    class f_336 func;
    f_337["write_PDB"] --> f_329["upper"]
    class f_337 func;
    class f_329 func;
    f_337["write_PDB"] -->|entity| f_311["set_structure"]
    class f_337 func;
    class f_311 func;
    f_337["write_PDB"] -->|fp| f_284["save"]
    class f_337 func;
    class f_284 func;
    f_9["__init__"] -->|verbose| f_91["_set_residues"]
    class f_9 func;
    class f_91 func;
    f_3["__deepcopy__"] --> f_173["get"]
    class f_3 func;
    class f_173 func;
    f_3["__deepcopy__"] --> f_13["__new__"]
    class f_3 func;
    class f_13 func;
    f_3["__deepcopy__"] --> f_143["copy"]
    class f_3 func;
    class f_143 func;
    f_289["setResAtmVws"] --> f_180["get_atoms"]
    class f_289 func;
    class f_180 func;
    f_289["setResAtmVws"] --> f_233["is_disordered"]
    class f_289 func;
    class f_233 func;
    f_289["setResAtmVws"] --> f_333["values"]
    class f_289 func;
    class f_333 func;
    f_3["__deepcopy__"] --> f_333["values"]
    class f_3 func;
    class f_333 func;
    f_76["_peptide_check"] --> f_173["get"]
    class f_76 func;
    class f_173 func;
    f_76["_peptide_check"] --> f_233["is_disordered"]
    class f_76 func;
    class f_233 func;
    f_76["_peptide_check"] -->|Natom<br>pCatom| f_21["_atm_dist_chk"]
    class f_76 func;
    class f_21 func;
    f_76["_peptide_check"] --> f_333["values"]
    class f_76 func;
    class f_333 func;
    f_76["_peptide_check"] -->|n<br>c| f_21["_atm_dist_chk"]
    class f_76 func;
    class f_21 func;
    f_137["clear_ic"] --> f_199["get_residues"]
    class f_137 func;
    class f_199 func;
    f_18["_add_residue"] -->|res| f_76["_peptide_check"]
    class f_18 func;
    class f_76 func;
    f_18["_add_residue"] --> f_259["pretty_str"]
    class f_18 func;
    class f_259 func;
    f_18["_add_residue"] -->|str| f_135["cast"]
    class f_18 func;
    class f_135 func;
    f_18["_add_residue"] --> f_315["split_akl"]
    class f_18 func;
    class f_315 func;
    f_91["_set_residues"] --> f_288["set"]
    class f_91 func;
    class f_288 func;
    f_91["_set_residues"] --> f_199["get_residues"]
    class f_91 func;
    class f_199 func;
    f_91["_set_residues"] --> f_233["is_disordered"]
    class f_91 func;
    class f_233 func;
    f_91["_set_residues"] --> f_333["values"]
    class f_91 func;
    class f_333 func;
    f_91["_set_residues"] -->|r<br>last_res<br>last_ord_res<br>verbose| f_18["_add_residue"]
    class f_91 func;
    class f_18 func;
    f_91["_set_residues"] --> f_326["update"]
    class f_91 func;
    class f_326 func;
    f_91["_set_residues"] -->|res<br>last_res<br>last_ord_res<br>verbose| f_18["_add_residue"]
    class f_91 func;
    class f_18 func;
    f_290["setResAtms"] --> f_180["get_atoms"]
    class f_290 func;
    class f_180 func;
    f_290["setResAtms"] --> f_233["is_disordered"]
    class f_290 func;
    class f_233 func;
    f_290["setResAtms"] --> f_333["values"]
    class f_290 func;
    class f_333 func;
    f_129["build_atomArray"] --> f_23["_build_rak_cache"]
    class f_129 func;
    class f_23 func;
    f_130["build_edraArrays"] --> f_143["copy"]
    class f_130 func;
    class f_143 func;
    f_130["build_edraArrays"] --> f_235["items"]
    class f_130 func;
    class f_235 func;
    f_130["build_edraArrays"] --> f_333["values"]
    class f_130 func;
    class f_333 func;
    f_130["build_edraArrays"] --> f_237["keys"]
    class f_130 func;
    class f_237 func;
    f_56["_hedraDict2chain"] --> f_180["get_atoms"]
    class f_56 func;
    class f_180 func;
    f_56["_hedraDict2chain"] --> f_233["is_disordered"]
    class f_56 func;
    class f_233 func;
    f_56["_hedraDict2chain"] --> f_333["values"]
    class f_56 func;
    class f_333 func;
    f_56["_hedraDict2chain"] --> f_315["split_akl"]
    class f_56 func;
    class f_315 func;
    f_56["_hedraDict2chain"] --> f_60["_link_dihedra"]
    class f_56 func;
    class f_60 func;
    f_56["_hedraDict2chain"] --> f_129["build_atomArray"]
    class f_56 func;
    class f_129 func;
    f_56["_hedraDict2chain"] --> f_235["items"]
    class f_56 func;
    class f_235 func;
    f_56["_hedraDict2chain"] --> f_173["get"]
    class f_56 func;
    class f_173 func;
    f_56["_hedraDict2chain"] -->|atm| f_217["has_id"]
    class f_56 func;
    class f_217 func;
    f_56["_hedraDict2chain"] -->|altloc| f_154["disordered_has_id"]
    class f_56 func;
    class f_154 func;
    f_56["_hedraDict2chain"] -->|newAtom| f_108["add"]
    class f_56 func;
    class f_108 func;
    f_56["_hedraDict2chain"] -->|disordered_atom| f_108["add"]
    class f_56 func;
    class f_108 func;
    f_56["_hedraDict2chain"] -->|newAtom| f_150["disordered_add"]
    class f_56 func;
    class f_150 func;
    f_56["_hedraDict2chain"] --> f_170["flag_disordered"]
    class f_56 func;
    class f_170 func;
    f_56["_hedraDict2chain"] -->|altloc| f_156["disordered_select"]
    class f_56 func;
    class f_156 func;
    f_56["_hedraDict2chain"] -->|bfac| f_297["set_bfactor"]
    class f_56 func;
    class f_297 func;
    f_56["_hedraDict2chain"] -->|occ| f_305["set_occupancy"]
    class f_56 func;
    class f_305 func;
    f_56["_hedraDict2chain"] --> f_206["get_serial_number"]
    class f_56 func;
    class f_206 func;
    f_56["_hedraDict2chain"] --> f_237["keys"]
    class f_56 func;
    class f_237 func;
    f_56["_hedraDict2chain"] --> f_130["build_edraArrays"]
    class f_56 func;
    class f_130 func;
    f_123["assemble_residues"] -->|workSet| f_245["multi_coord_space"]
    class f_123 func;
    class f_245 func;
    f_124["assemble_residues_ser"] --> f_122["assemble"]
    class f_124 func;
    class f_122 func;
    f_124["assemble_residues_ser"] --> f_288["set"]
    class f_124 func;
    class f_288 func;
    f_124["assemble_residues_ser"] --> f_237["keys"]
    class f_124 func;
    class f_237 func;
    f_224["init_edra"] --> f_31["_create_edra"]
    class f_224 func;
    class f_31 func;
    f_224["init_edra"] --> f_129["build_atomArray"]
    class f_224 func;
    class f_129 func;
    f_224["init_edra"] --> f_237["keys"]
    class f_224 func;
    class f_237 func;
    f_224["init_edra"] --> f_130["build_edraArrays"]
    class f_224 func;
    class f_130 func;
    f_222["init_atom_coords"] -->|a4shift| f_108["add"]
    class f_222 func;
    class f_108 func;
    f_222["init_atom_coords"] --> f_247["multi_rot_Z"]
    class f_222 func;
    class f_247 func;
    f_327["update_dCoordSpace"] --> f_123["assemble_residues"]
    class f_327 func;
    class f_123 func;
    f_327["update_dCoordSpace"] -->|workSet| f_245["multi_coord_space"]
    class f_327 func;
    class f_245 func;
    f_262["propagate_changes"] --> f_288["set"]
    class f_262 func;
    class f_288 func;
    f_262["propagate_changes"] -->|pos| f_108["add"]
    class f_262 func;
    class f_108 func;
    f_230["internal_to_atom_coordinates"] --> f_262["propagate_changes"]
    class f_230 func;
    class f_262 func;
    f_230["internal_to_atom_coordinates"] --> f_222["init_atom_coords"]
    class f_230 func;
    class f_222 func;
    f_230["internal_to_atom_coordinates"] --> f_123["assemble_residues"]
    class f_230 func;
    class f_123 func;
    f_230["internal_to_atom_coordinates"] --> f_333["values"]
    class f_230 func;
    class f_333 func;
    f_230["internal_to_atom_coordinates"] --> f_259["pretty_str"]
    class f_230 func;
    class f_259 func;
    f_230["internal_to_atom_coordinates"] --> f_315["split_akl"]
    class f_230 func;
    class f_315 func;
    f_230["internal_to_atom_coordinates"] --> f_124["assemble_residues_ser"]
    class f_230 func;
    class f_124 func;
    f_125["atom_to_internal_coordinates"] --> f_224["init_edra"]
    class f_125 func;
    class f_224 func;
    f_125["atom_to_internal_coordinates"] --> f_249["norm"]
    class f_125 func;
    class f_249 func;
    f_125["atom_to_internal_coordinates"] -->|dha| f_245["multi_coord_space"]
    class f_125 func;
    class f_245 func;
    f_125["atom_to_internal_coordinates"] --> f_92["_spec_glyCB"]
    class f_125 func;
    class f_92 func;
    f_92["_spec_glyCB"] --> f_333["values"]
    class f_92 func;
    class f_333 func;
    f_92["_spec_glyCB"] --> f_266["rak"]
    class f_92 func;
    class f_266 func;
    f_92["_spec_glyCB"] --> f_173["get"]
    class f_92 func;
    class f_173 func;
    f_102["_write_mtx"] --> f_336["write"]
    class f_102 func;
    class f_336 func;
    f_99["_writeSCAD_dihed"] --> f_336["write"]
    class f_99 func;
    class f_336 func;
    f_99["_writeSCAD_dihed"] -->|fp| f_102["_write_mtx"]
    class f_99 func;
    class f_102 func;
    f_101["_write_SCAD"] --> f_336["write"]
    class f_101 func;
    class f_336 func;
    f_101["_write_SCAD"] --> f_235["items"]
    class f_101 func;
    class f_235 func;
    f_101["_write_SCAD"] --> f_288["set"]
    class f_101 func;
    class f_288 func;
    f_101["_write_SCAD"] --> f_138["clear_transforms"]
    class f_101 func;
    class f_138 func;
    f_101["_write_SCAD"] --> f_122["assemble"]
    class f_101 func;
    class f_122 func;
    f_101["_write_SCAD"] --> f_232["is_backbone"]
    class f_101 func;
    class f_232 func;
    f_101["_write_SCAD"] -->|fp<br>d<br>hedraNdx<br>hedraSet| f_99["_writeSCAD_dihed"]
    class f_101 func;
    class f_99 func;
    f_101["_write_SCAD"] --> f_108["add"]
    class f_101 func;
    class f_108 func;
    f_101["_write_SCAD"] --> f_293["set_accuracy_95"]
    class f_101 func;
    class f_293 func;
    f_101["_write_SCAD"] -->|atm| f_173["get"]
    class f_101 func;
    class f_173 func;
    f_101["_write_SCAD"] -->|res| f_173["get"]
    class f_101 func;
    class f_173 func;
    f_101["_write_SCAD"] -->|ak| f_108["add"]
    class f_101 func;
    class f_108 func;
    f_101["_write_SCAD"] -->|atom_str| f_336["write"]
    class f_101 func;
    class f_336 func;
    f_101["_write_SCAD"] -->|atom_done_str| f_336["write"]
    class f_101 func;
    class f_336 func;
    f_101["_write_SCAD"] -->|wstr| f_336["write"]
    class f_101 func;
    class f_336 func;
    f_101["_write_SCAD"] --> f_230["internal_to_atom_coordinates"]
    class f_101 func;
    class f_230 func;
    f_101["_write_SCAD"] --> f_142["coord_space"]
    class f_101 func;
    class f_142 func;
    f_101["_write_SCAD"] -->|fp<br>mtr| f_102["_write_mtx"]
    class f_101 func;
    class f_102 func;
    f_157["distance_plot"] --> f_249["norm"]
    class f_157 func;
    class f_249 func;
    f_242["make_extended"] --> f_294["set_angle"]
    class f_242 func;
    class f_294 func;
    f_9["__init__"] --> f_288["set"]
    class f_9 func;
    class f_288 func;
    f_9["__init__"] --> f_180["get_atoms"]
    class f_9 func;
    class f_180 func;
    f_9["__init__"] --> f_17["_add_atom"]
    class f_9 func;
    class f_17 func;
    f_9["__init__"] --> f_333["values"]
    class f_9 func;
    class f_333 func;
    f_9["__init__"] -->|atm| f_17["_add_atom"]
    class f_9 func;
    class f_17 func;
    f_9["__init__"] -->|atom| f_17["_add_atom"]
    class f_9 func;
    class f_17 func;
    f_9["__init__"] --> f_23["_build_rak_cache"]
    class f_9 func;
    class f_23 func;
    f_3["__deepcopy__"] --> f_326["update"]
    class f_3 func;
    class f_326 func;
    f_23["_build_rak_cache"] -->|atmName| f_173["get"]
    class f_23 func;
    class f_173 func;
    f_17["_add_atom"] -->|atm| f_266["rak"]
    class f_17 func;
    class f_266 func;
    f_17["_add_atom"] -->|ak| f_108["add"]
    class f_17 func;
    class f_108 func;
    f_60["_link_dihedra"] --> f_333["values"]
    class f_60 func;
    class f_333 func;
    f_60["_link_dihedra"] --> f_326["update"]
    class f_60 func;
    class f_326 func;
    f_60["_link_dihedra"] --> f_23["_build_rak_cache"]
    class f_60 func;
    class f_23 func;
    f_60["_link_dihedra"] --> f_315["split_akl"]
    class f_60 func;
    class f_315 func;
    f_299["set_flexible"] --> f_333["values"]
    class f_299 func;
    class f_333 func;
    f_299["set_flexible"] --> f_164["endswith"]
    class f_299 func;
    class f_164 func;
    f_299["set_flexible"] --> f_316["startswith"]
    class f_299 func;
    class f_316 func;
    f_300["set_hbond"] --> f_333["values"]
    class f_300 func;
    class f_333 func;
    f_32["_default_startpos"] -->|akl| f_173["get"]
    class f_32 func;
    class f_173 func;
    f_53["_get_startpos"] --> f_32["_default_startpos"]
    class f_53 func;
    class f_32 func;
    f_138["clear_transforms"] --> f_333["values"]
    class f_138 func;
    class f_333 func;
    f_122["assemble"] --> f_315["split_akl"]
    class f_122 func;
    class f_315 func;
    f_122["assemble"] --> f_266["rak"]
    class f_122 func;
    class f_266 func;
    f_122["assemble"] --> f_32["_default_startpos"]
    class f_122 func;
    class f_32 func;
    f_122["assemble"] --> f_53["_get_startpos"]
    class f_122 func;
    class f_53 func;
    f_122["assemble"] -->|HKT| f_135["cast"]
    class f_122 func;
    class f_135 func;
    f_122["assemble"] --> f_258["pop"]
    class f_122 func;
    class f_258 func;
    f_122["assemble"] -->|h1k| f_173["get"]
    class f_122 func;
    class f_173 func;
    f_122["assemble"] --> f_142["coord_space"]
    class f_122 func;
    class f_142 func;
    f_122["assemble"] -->|a| f_173["get"]
    class f_122 func;
    class f_173 func;
    f_315["split_akl"] --> f_288["set"]
    class f_315 func;
    class f_288 func;
    f_315["split_akl"] -->|ak2| f_114["altloc_match"]
    class f_315 func;
    class f_114 func;
    f_315["split_akl"] -->|altloc| f_108["add"]
    class f_315 func;
    class f_108 func;
    f_37["_gen_edra"] -->|tlst| f_315["split_akl"]
    class f_37 func;
    class f_315 func;
    f_31["_create_edra"] --> f_266["rak"]
    class f_31 func;
    class f_266 func;
    f_31["_create_edra"] --> f_315["split_akl"]
    class f_31 func;
    class f_315 func;
    f_31["_create_edra"] -->|ak| f_108["add"]
    class f_31 func;
    class f_108 func;
    f_31["_create_edra"] -->|ak| f_114["altloc_match"]
    class f_31 func;
    class f_114 func;
    f_31["_create_edra"] -->|rn_ak| f_108["add"]
    class f_31 func;
    class f_108 func;
    f_31["_create_edra"] --> f_37["_gen_edra"]
    class f_31 func;
    class f_37 func;
    f_31["_create_edra"] --> f_173["get"]
    class f_31 func;
    class f_173 func;
    f_31["_create_edra"] -->|nCB| f_108["add"]
    class f_31 func;
    class f_108 func;
    f_31["_create_edra"] -->|atom| f_266["rak"]
    class f_31 func;
    class f_266 func;
    f_31["_create_edra"] -->|r_edra| f_37["_gen_edra"]
    class f_31 func;
    class f_37 func;
    f_31["_create_edra"] --> f_108["add"]
    class f_31 func;
    class f_108 func;
    f_31["_create_edra"] -->|sCB| f_108["add"]
    class f_31 func;
    class f_108 func;
    f_31["_create_edra"] -->|htpl| f_37["_gen_edra"]
    class f_31 func;
    class f_37 func;
    f_31["_create_edra"] -->|dtpl| f_37["_gen_edra"]
    class f_31 func;
    class f_37 func;
    f_31["_create_edra"] --> f_90["_set_hedra"]
    class f_31 func;
    class f_90 func;
    f_31["_create_edra"] -->|verbose| f_60["_link_dihedra"]
    class f_31 func;
    class f_60 func;
    f_31["_create_edra"] --> f_235["items"]
    class f_31 func;
    class f_235 func;
    f_75["_pdb_atom_string"] --> f_233["is_disordered"]
    class f_75 func;
    class f_233 func;
    f_75["_pdb_atom_string"] --> f_333["values"]
    class f_75 func;
    class f_333 func;
    f_255["pdb_residue_string"] --> f_75["_pdb_atom_string"]
    class f_255 func;
    class f_75 func;
    f_84["_residue_string"] --> f_203["get_segid"]
    class f_84 func;
    class f_203 func;
    f_84["_residue_string"] --> f_185["get_full_id"]
    class f_84 func;
    class f_185 func;
    f_1["Main_Script"] --> f_248["namedtuple"]
    class f_1 main;
    class f_248 func;
    f_103["_write_pic_bfac"] -->|atm| f_266["rak"]
    class f_103 func;
    class f_266 func;
    f_103["_write_pic_bfac"] --> f_181["get_bfactor"]
    class f_103 func;
    class f_181 func;
    f_100["_write_PIC"] --> f_84["_residue_string"]
    class f_100 func;
    class f_84 func;
    f_100["_write_PIC"] --> f_256["pick_angle"]
    class f_100 func;
    class f_256 func;
    f_100["_write_PIC"] --> f_75["_pdb_atom_string"]
    class f_100 func;
    class f_75 func;
    f_100["_write_PIC"] --> f_333["values"]
    class f_100 func;
    class f_333 func;
    f_100["_write_PIC"] --> f_126["bits"]
    class f_100 func;
    class f_126 func;
    f_100["_write_PIC"] --> f_180["get_atoms"]
    class f_100 func;
    class f_180 func;
    f_100["_write_PIC"] --> f_233["is_disordered"]
    class f_100 func;
    class f_233 func;
    f_100["_write_PIC"] -->|s<br>col| f_103["_write_pic_bfac"]
    class f_100 func;
    class f_103 func;
    f_100["_write_PIC"] -->|atm<br>s<br>col| f_103["_write_pic_bfac"]
    class f_100 func;
    class f_103 func;
    f_100["_write_PIC"] -->|a<br>s<br>col| f_103["_write_pic_bfac"]
    class f_100 func;
    class f_103 func;
    f_39["_get_ak_tuple"] --> f_314["split"]
    class f_39 func;
    class f_314 func;
    f_39["_get_ak_tuple"] --> f_266["rak"]
    class f_39 func;
    class f_266 func;
    f_39["_get_ak_tuple"] -->|a| f_266["rak"]
    class f_39 func;
    class f_266 func;
    f_40["_get_angle_for_tuple"] --> f_173["get"]
    class f_40 func;
    class f_173 func;
    f_40["_get_angle_for_tuple"] -->|DKT<br>angle_key| f_135["cast"]
    class f_40 func;
    class f_135 func;
    f_40["_get_angle_for_tuple"] -->|HKT<br>angle_key| f_135["cast"]
    class f_40 func;
    class f_135 func;
    f_256["pick_angle"] -->|angle_key| f_40["_get_angle_for_tuple"]
    class f_256 func;
    class f_40 func;
    f_256["pick_angle"] -->|EKT| f_135["cast"]
    class f_256 func;
    class f_135 func;
    f_256["pick_angle"] --> f_39["_get_ak_tuple"]
    class f_256 func;
    class f_39 func;
    f_256["pick_angle"] -->|str<br>angle_key| f_135["cast"]
    class f_256 func;
    class f_135 func;
    f_256["pick_angle"] --> f_266["rak"]
    class f_256 func;
    class f_266 func;
    f_256["pick_angle"] --> f_173["get"]
    class f_256 func;
    class f_173 func;
    f_256["pick_angle"] --> f_316["startswith"]
    class f_256 func;
    class f_316 func;
    f_256["pick_angle"] -->|a| f_266["rak"]
    class f_256 func;
    class f_266 func;
    f_256["pick_angle"] -->|DKT| f_135["cast"]
    class f_256 func;
    class f_135 func;
    f_256["pick_angle"] -->|tklst| f_173["get"]
    class f_256 func;
    class f_173 func;
    f_178["get_angle"] -->|angle_key| f_256["pick_angle"]
    class f_178 func;
    class f_256 func;
    f_294["set_angle"] -->|angle_key| f_256["pick_angle"]
    class f_294 func;
    class f_256 func;
    f_294["set_angle"] -->|v| f_116["angle_dif"]
    class f_294 func;
    class f_116 func;
    f_294["set_angle"] -->|edron<br>delta| f_33["_do_bond_rotate"]
    class f_294 func;
    class f_33 func;
    f_127["bond_rotate"] -->|angle_key| f_256["pick_angle"]
    class f_127 func;
    class f_256 func;
    f_127["bond_rotate"] -->|base<br>delta| f_33["_do_bond_rotate"]
    class f_127 func;
    class f_33 func;
    f_128["bond_set"] -->|angle_key| f_256["pick_angle"]
    class f_128 func;
    class f_256 func;
    f_128["bond_set"] -->|val| f_116["angle_dif"]
    class f_128 func;
    class f_116 func;
    f_128["bond_set"] -->|base<br>delta| f_33["_do_bond_rotate"]
    class f_128 func;
    class f_33 func;
    f_257["pick_length"] -->|BKT| f_135["cast"]
    class f_257 func;
    class f_135 func;
    f_257["pick_length"] -->|ak_spec| f_39["_get_ak_tuple"]
    class f_257 func;
    class f_39 func;
    f_257["pick_length"] --> f_235["items"]
    class f_257 func;
    class f_235 func;
    f_189["get_length"] -->|ak_spec| f_257["pick_length"]
    class f_189 func;
    class f_257 func;
    f_303["set_length"] -->|ak_spec| f_257["pick_length"]
    class f_303 func;
    class f_257 func;
    f_120["applyMtx"] -->|k| f_173["get"]
    class f_120 func;
    class f_173 func;
    f_120["applyMtx"] --> f_237["keys"]
    class f_120 func;
    class f_237 func;
    f_172["gen_tuple"] --> f_314["split"]
    class f_172 func;
    class f_314 func;
    f_9["__init__"] -->|atomkeys| f_171["gen_key"]
    class f_9 func;
    class f_171 func;
    f_9["__init__"] --> f_108["add"]
    class f_9 func;
    class f_108 func;
    f_9["__init__"] --> f_144["cr_class"]
    class f_9 func;
    class f_144 func;
    f_7["__gt__"] -->|other| f_27["_cmp"]
    class f_7 func;
    class f_27 func;
    f_7["__gt__"] -->|rslt| f_135["cast"]
    class f_7 func;
    class f_135 func;
    f_5["__ge__"] -->|other| f_27["_cmp"]
    class f_5 func;
    class f_27 func;
    f_5["__ge__"] -->|rslt| f_135["cast"]
    class f_5 func;
    class f_135 func;
    f_12["__lt__"] -->|other| f_27["_cmp"]
    class f_12 func;
    class f_27 func;
    f_12["__lt__"] -->|rslt| f_135["cast"]
    class f_12 func;
    class f_135 func;
    f_11["__le__"] -->|other| f_27["_cmp"]
    class f_11 func;
    class f_27 func;
    f_11["__le__"] -->|rslt| f_135["cast"]
    class f_11 func;
    class f_135 func;
    f_303["set_length"] --> f_57["_invalidate_atoms"]
    class f_303 func;
    class f_57 func;
    f_9["__init__"] -->|HKT| f_135["cast"]
    class f_9 func;
    class f_135 func;
    f_9["__init__"] --> f_89["_setPrimary"]
    class f_9 func;
    class f_89 func;
    f_48["_get_hedron"] -->|id3| f_173["get"]
    class f_48 func;
    class f_173 func;
    f_90["_set_hedra"] -->|res<br>h1key| f_48["_get_hedron"]
    class f_90 func;
    class f_48 func;
    f_90["_set_hedra"] -->|HKT| f_135["cast"]
    class f_90 func;
    class f_135 func;
    f_90["_set_hedra"] -->|res<br>h2key| f_48["_get_hedron"]
    class f_90 func;
    class f_48 func;
    f_117["angle_pop_sd"] -->|alst<br>avg| f_116["angle_dif"]
    class f_117 func;
    class f_116 func;
    f_149["difference"] --> f_116["angle_dif"]
    class f_149 func;
    class f_116 func;
    f_126["bits"] --> f_173["get"]
    class f_126 func;
    class f_173 func;
    f_9["__init__"] --> f_243["map"]
    class f_9 func;
    class f_243 func;
    f_9["__init__"] -->|k| f_173["get"]
    class f_9 func;
    class f_173 func;
    f_9["__init__"] --> f_173["get"]
    class f_9 func;
    class f_173 func;
    f_27["_cmp"] -->|s| f_173["get"]
    class f_27 func;
    class f_173 func;
    f_27["_cmp"] -->|o| f_173["get"]
    class f_27 func;
    class f_173 func;
    f_284["save"] -->|fp<br>select<br>preserve_atom_numbering| f_87["_save_structure"]
    class f_284 func;
    class f_87 func;
    f_284["save"] -->|fp| f_86["_save_dict"]
    class f_284 func;
    class f_86 func;
    f_86["_save_dict"] -->|key| f_314["split"]
    class f_86 func;
    class f_314 func;
    f_86["_save_dict"] --> f_235["items"]
    class f_86 func;
    class f_235 func;
    f_86["_save_dict"] -->|i| f_220["index"]
    class f_86 func;
    class f_220 func;
    f_86["_save_dict"] --> f_336["write"]
    class f_86 func;
    class f_336 func;
    f_86["_save_dict"] -->|value_no_list| f_36["_format_mmcif_col"]
    class f_86 func;
    class f_36 func;
    f_86["_save_dict"] -->|val| f_81["_requires_quote"]
    class f_86 func;
    class f_81 func;
    f_86["_save_dict"] -->|val| f_80["_requires_newline"]
    class f_86 func;
    class f_80 func;
    f_86["_save_dict"] --> f_36["_format_mmcif_col"]
    class f_86 func;
    class f_36 func;
    f_36["_format_mmcif_col"] -->|val| f_80["_requires_newline"]
    class f_36 func;
    class f_80 func;
    f_36["_format_mmcif_col"] -->|val| f_81["_requires_quote"]
    class f_36 func;
    class f_81 func;
    f_81["_requires_quote"] --> f_316["startswith"]
    class f_81 func;
    class f_316 func;
    f_87["_save_structure"] --> f_191["get_list"]
    class f_87 func;
    class f_191 func;
    f_87["_save_structure"] -->|model| f_106["accept_model"]
    class f_87 func;
    class f_106 func;
    f_87["_save_structure"] -->|chain| f_105["accept_chain"]
    class f_87 func;
    class f_105 func;
    f_87["_save_structure"] --> f_188["get_id"]
    class f_87 func;
    class f_188 func;
    f_87["_save_structure"] --> f_215["get_unpacked_list"]
    class f_87 func;
    class f_215 func;
    f_87["_save_structure"] -->|residue| f_107["accept_residue"]
    class f_87 func;
    class f_107 func;
    f_87["_save_structure"] --> f_200["get_resname"]
    class f_87 func;
    class f_200 func;
    f_87["_save_structure"] -->|entity_id| f_50["_get_label_asym_id"]
    class f_87 func;
    class f_50 func;
    f_87["_save_structure"] -->|atom| f_104["accept_atom"]
    class f_87 func;
    class f_104 func;
    f_87["_save_structure"] --> f_206["get_serial_number"]
    class f_87 func;
    class f_206 func;
    f_87["_save_structure"] --> f_318["strip"]
    class f_87 func;
    class f_318 func;
    f_87["_save_structure"] --> f_193["get_name"]
    class f_87 func;
    class f_193 func;
    f_87["_save_structure"] --> f_177["get_altloc"]
    class f_87 func;
    class f_177 func;
    f_87["_save_structure"] --> f_184["get_coord"]
    class f_87 func;
    class f_184 func;
    f_87["_save_structure"] --> f_194["get_occupancy"]
    class f_87 func;
    class f_194 func;
    f_87["_save_structure"] --> f_181["get_bfactor"]
    class f_87 func;
    class f_181 func;
    f_87["_save_structure"] -->|c| f_272["replace"]
    class f_87 func;
    class f_272 func;
    f_87["_save_structure"] -->|out_file| f_86["_save_dict"]
    class f_87 func;
    class f_86 func;
    f_44["_get_biomoltrans"] --> f_316["startswith"]
    class f_44 func;
    class f_316 func;
    f_44["_get_biomoltrans"] --> f_314["split"]
    class f_44 func;
    class f_314 func;
    f_44["_get_biomoltrans"] --> f_318["strip"]
    class f_44 func;
    class f_318 func;
    f_44["_get_biomoltrans"] --> f_272["replace"]
    class f_44 func;
    class f_272 func;
    f_49["_get_journal"] -->|line| f_285["search"]
    class f_49 func;
    class f_285 func;
    f_49["_get_journal"] --> f_239["lower"]
    class f_49 func;
    class f_239 func;
    f_52["_get_references"] -->|line| f_285["search"]
    class f_52 func;
    class f_285 func;
    f_52["_get_references"] --> f_239["lower"]
    class f_52 func;
    class f_239 func;
    f_35["_format_date"] --> f_220["index"]
    class f_35 func;
    class f_220 func;
    f_65["_nice_case"] --> f_239["lower"]
    class f_65 func;
    class f_239 func;
    f_65["_nice_case"] --> f_329["upper"]
    class f_65 func;
    class f_329 func;
    f_252["parse_pdb_header"] -->|infile| f_121["as_handle"]
    class f_252 func;
    class f_121 func;
    f_252["parse_pdb_header"] -->|header| f_71["_parse_pdb_header_list"]
    class f_252 func;
    class f_71 func;
    f_72["_parse_remark_465"] --> f_314["split"]
    class f_72 func;
    class f_314 func;
    f_71["_parse_pdb_header_list"] -->|header| f_52["_get_references"]
    class f_71 func;
    class f_52 func;
    f_71["_parse_pdb_header_list"] -->|header| f_49["_get_journal"]
    class f_71 func;
    class f_49 func;
    f_71["_parse_pdb_header_list"] -->|header| f_44["_get_biomoltrans"]
    class f_71 func;
    class f_44 func;
    f_71["_parse_pdb_header_list"] --> f_318["strip"]
    class f_71 func;
    class f_318 func;
    f_71["_parse_pdb_header_list"] --> f_239["lower"]
    class f_71 func;
    class f_239 func;
    f_71["_parse_pdb_header_list"] -->|tail| f_25["_chop_end_codes"]
    class f_71 func;
    class f_25 func;
    f_71["_parse_pdb_header_list"] --> f_236["join"]
    class f_71 func;
    class f_236 func;
    f_71["_parse_pdb_header_list"] -->|tail| f_285["search"]
    class f_71 func;
    class f_285 func;
    f_71["_parse_pdb_header_list"] --> f_35["_format_date"]
    class f_71 func;
    class f_35 func;
    f_71["_parse_pdb_header_list"] --> f_65["_nice_case"]
    class f_71 func;
    class f_65 func;
    f_71["_parse_pdb_header_list"] -->|tail| f_26["_chop_end_misc"]
    class f_71 func;
    class f_26 func;
    f_71["_parse_pdb_header_list"] -->|tt| f_285["search"]
    class f_71 func;
    class f_285 func;
    f_71["_parse_pdb_header_list"] --> f_314["split"]
    class f_71 func;
    class f_314 func;
    f_71["_parse_pdb_header_list"] -->|hh| f_285["search"]
    class f_71 func;
    class f_285 func;
    f_71["_parse_pdb_header_list"] --> f_25["_chop_end_codes"]
    class f_71 func;
    class f_25 func;
    f_71["_parse_pdb_header_list"] --> f_316["startswith"]
    class f_71 func;
    class f_316 func;
    f_71["_parse_pdb_header_list"] -->|tail| f_72["_parse_remark_465"]
    class f_71 func;
    class f_72 func;
    f_71["_parse_pdb_header_list"] --> f_272["replace"]
    class f_71 func;
    class f_272 func;
    f_9["__init__"] --> f_83["_reset_properties"]
    class f_9 func;
    class f_83 func;
    f_296["set_atoms"] -->|fix_coord<br>mov_coord| f_288["set"]
    class f_296 func;
    class f_288 func;
    f_288["set"] --> f_83["_reset_properties"]
    class f_288 func;
    class f_83 func;
    f_281["run"] --> f_143["copy"]
    class f_281 func;
    class f_143 func;
    f_281["run"] -->|coords_ref<br>coords| f_265["qcp"]
    class f_281 func;
    class f_265 func;
    f_240["m2rotaxis"] --> f_250["normalize"]
    class f_240 func;
    class f_250 func;
    f_334["vector_to_axis"] --> f_251["normalized"]
    class f_334 func;
    class f_251 func;
    f_334["vector_to_axis"] --> f_249["norm"]
    class f_334 func;
    class f_249 func;
    f_334["vector_to_axis"] -->|point| f_115["angle"]
    class f_334 func;
    class f_115 func;
    f_278["rotaxis2m"] --> f_251["normalized"]
    class f_278 func;
    class f_251 func;
    f_278["rotaxis2m"] --> f_179["get_array"]
    class f_278 func;
    class f_179 func;
    f_271["refmat"] --> f_251["normalized"]
    class f_271 func;
    class f_251 func;
    f_271["refmat"] --> f_249["norm"]
    class f_271 func;
    class f_249 func;
    f_271["refmat"] --> f_250["normalize"]
    class f_271 func;
    class f_250 func;
    f_271["refmat"] --> f_179["get_array"]
    class f_271 func;
    class f_179 func;
    f_279["rotmat"] -->|q| f_271["refmat"]
    class f_279 func;
    class f_271 func;
    f_279["rotmat"] -->|p| f_271["refmat"]
    class f_279 func;
    class f_271 func;
    f_133["calc_angle"] -->|v3| f_115["angle"]
    class f_133 func;
    class f_115 func;
    f_134["calc_dihedral"] -->|v| f_115["angle"]
    class f_134 func;
    class f_115 func;
    f_134["calc_dihedral"] -->|w| f_115["angle"]
    class f_134 func;
    class f_115 func;
    f_250["normalize"] --> f_249["norm"]
    class f_250 func;
    class f_249 func;
    f_251["normalized"] --> f_143["copy"]
    class f_251 func;
    class f_143 func;
    f_251["normalized"] --> f_250["normalize"]
    class f_251 func;
    class f_250 func;
    f_115["angle"] --> f_249["norm"]
    class f_115 func;
    class f_249 func;
    f_207["get_spherical_coordinates"] -->|xyz| f_249["norm"]
    class f_207 func;
    class f_249 func;
    f_207["get_spherical_coordinates"] --> f_43["_get_azimuth"]
    class f_207 func;
    class f_43 func;
    f_142["coord_space"] -->|tm| f_302["set_homog_trans_mtx"]
    class f_142 func;
    class f_302 func;
    f_142["coord_space"] -->|p| f_207["get_spherical_coordinates"]
    class f_142 func;
    class f_207 func;
    f_142["coord_space"] -->|mrz| f_292["set_Z_homog_rot_mtx"]
    class f_142 func;
    class f_292 func;
    f_142["coord_space"] -->|mry| f_291["set_Y_homog_rot_mtx"]
    class f_142 func;
    class f_291 func;
    f_142["coord_space"] --> f_43["_get_azimuth"]
    class f_142 func;
    class f_43 func;
    f_142["coord_space"] -->|mrz2| f_292["set_Z_homog_rot_mtx"]
    class f_142 func;
    class f_292 func;
    f_142["coord_space"] -->|azimuth2<br>mrz2| f_292["set_Z_homog_rot_mtx"]
    class f_142 func;
    class f_292 func;
    f_245["multi_coord_space"] -->|p| f_249["norm"]
    class f_245 func;
    class f_249 func;
    f_245["multi_coord_space"] --> f_247["multi_rot_Z"]
    class f_245 func;
    class f_247 func;
    f_245["multi_coord_space"] --> f_246["multi_rot_Y"]
    class f_245 func;
    class f_246 func;
    f_245["multi_coord_space"] -->|azimuth2| f_247["multi_rot_Z"]
    class f_245 func;
    class f_247 func;
    f_245["multi_coord_space"] -->|polar_angle| f_246["multi_rot_Y"]
    class f_245 func;
    class f_246 func;
    f_245["multi_coord_space"] -->|azimuth| f_247["multi_rot_Z"]
    class f_245 func;
    class f_247 func;
```

### ✂️ Navigator: Snippet Extractor
Want to use a specific function without the whole library? Here is the **Dependency Closure** for **Top 20** key functions.
#### To extract `__init__`:
> You need these **137** components:
`__init__, _add_atom, _add_residue, _assign_atom_mass, _assign_element, _atm_dist_chk, _build_rak_cache, _compute_sphere, _format_mmcif_col, _generate_full_id, _get_atom_radius, _get_cb, _get_gly_cb_vector, _get_label_asym_id, _is_completely_disordered, _make_dssp_dict, _make_fragment_list, _map, _map_fragment_list, _peptide_check, _read_fragments, _read_vertex_array, _requires_newline, _requires_quote, _reset_full_id, _reset_properties, _revert_write, _save_dict, _save_structure, _setPrimary, _set_residues, _splitline, _test_equivalence, _tokenize, _translate_id, accept_atom, accept_chain, accept_model, accept_residue, add, add_residue, altloc_match, angle, annotate, as_handle, build_peptides, ca_depth, cast, close, copy, cr_class, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_id_list, disordered_get_list, disordered_has_id, disordered_select, dssp_dict_from_pdb_file, endswith, find, flag_disorder, gen_key, get, get_altloc, get_atoms, get_bfactor, get_chains, get_coord, get_full_id, get_id, get_level, get_list, get_models, get_name, get_occupancy, get_parent, get_residues, get_resname, get_serial_number, get_surface, get_unpacked_list, get_vector, has_id, index, init_chain, init_model, init_residue, init_seg, init_structure, is_aa, is_disordered, items, join, left_multiply, lower, make_dssp_dict, map, min_dist, norm, normalize, pretty_str, process_asa_data, process_rsa_data, psea, psea2HEC, qcp, rak, readlines, replace, residue_depth, rotaxis, rstrip, run, run_naccess, run_psea, save, search, seek, set, set_coord, set_parent, set_structure, sort, split, split_akl, startswith, strip, tell, truncate, unfold_entities, update, upper, values, version, write`

#### To extract `_hedraDict2chain`:
> You need these **45** components:
`_build_rak_cache, _generate_full_id, _hedraDict2chain, _link_dihedra, _reset_full_id, _reset_properties, _translate_id, add, altloc_match, build_atomArray, build_edraArrays, copy, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, flag_disordered, get, get_altloc, get_atoms, get_chains, get_coord, get_full_id, get_id, get_models, get_occupancy, get_parent, get_residues, get_resname, get_serial_number, has_id, is_disordered, items, keys, set, set_bfactor, set_coord, set_occupancy, set_parent, split_akl, update, values`

#### To extract `_save_structure`:
> You need these **46** components:
`_format_mmcif_col, _generate_full_id, _get_label_asym_id, _requires_newline, _requires_quote, _reset_full_id, _save_dict, _save_structure, _translate_id, accept_atom, accept_chain, accept_model, accept_residue, add, copy, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, get_altloc, get_bfactor, get_coord, get_full_id, get_id, get_list, get_name, get_occupancy, get_parent, get_resname, get_serial_number, get_unpacked_list, has_id, index, is_disordered, items, replace, set_coord, set_parent, split, startswith, strip, values, write`

#### To extract `get_structure`:
> You need these **81** components:
`_build_structure, _chop_end_codes, _chop_end_misc, _format_date, _generate_full_id, _get_biomoltrans, _get_header, _get_journal, _get_references, _handle_PDB_exception, _is_completely_disordered, _nice_case, _parse, _parse_atom_from, _parse_coordinates, _parse_header_from, _parse_pdb_header_list, _parse_remark_465, _parse_residue_id_from, _parse_resolution_from, _reset_full_id, _translate_id, _update_header_entry, add, as_handle, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, find, flag_disorder, flag_disordered, get, get_altloc, get_coord, get_full_id, get_fullname, get_id, get_level, get_list, get_occupancy, get_parent, get_resname, get_structure, get_unpacked_list, has_id, index, init_atom, init_chain, init_model, init_residue, init_seg, init_structure, is_disordered, join, lower, map, readlines, replace, rstrip, search, seek, set_anisou, set_coord, set_header, set_line_counter, set_parent, set_sigatm, set_siguij, set_symmetry, split, startswith, strip, unfold_entities, update, upper, values`

#### To extract `_parse_pdb_header_list`:
> You need these **21** components:
`_chop_end_codes, _chop_end_misc, _format_date, _get_biomoltrans, _get_journal, _get_references, _nice_case, _parse_pdb_header_list, _parse_remark_465, get_level, get_parent, index, join, lower, replace, search, split, startswith, strip, unfold_entities, upper`

#### To extract `_build_structure`:
> You need these **45** components:
`_build_structure, _generate_full_id, _is_completely_disordered, _reset_full_id, _translate_id, add, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, flag_disordered, get_altloc, get_coord, get_full_id, get_fullname, get_id, get_list, get_occupancy, get_parent, get_resname, get_unpacked_list, has_id, init_atom, init_chain, init_model, init_residue, init_seg, init_structure, is_disordered, map, set_anisou, set_coord, set_line_counter, set_parent, set_symmetry, startswith, strip, update, upper, values`

#### To extract `_parse_coordinates`:
> You need these **45** components:
`_generate_full_id, _handle_PDB_exception, _is_completely_disordered, _parse_coordinates, _reset_full_id, _translate_id, add, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, flag_disordered, get_altloc, get_coord, get_full_id, get_fullname, get_id, get_list, get_occupancy, get_parent, get_resname, get_unpacked_list, has_id, init_atom, init_chain, init_model, init_residue, init_seg, is_disordered, rstrip, set_anisou, set_coord, set_line_counter, set_parent, set_sigatm, set_siguij, split, strip, upper, values`

#### To extract `read_PIC`:
> You need these **84** components:
`_build_structure, _chop_end_codes, _chop_end_misc, _format_date, _generate_full_id, _get_biomoltrans, _get_header, _get_journal, _get_references, _handle_PDB_exception, _is_completely_disordered, _nice_case, _parse, _parse_atom_from, _parse_coordinates, _parse_header_from, _parse_pdb_header_list, _parse_remark_465, _parse_residue_id_from, _parse_resolution_from, _reset_full_id, _reset_properties, _translate_id, _update_header_entry, add, as_handle, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, find, flag_disorder, flag_disordered, get, get_altloc, get_coord, get_full_id, get_fullname, get_id, get_level, get_list, get_occupancy, get_parent, get_resname, get_structure, get_unpacked_list, has_id, index, init_atom, init_chain, init_model, init_residue, init_seg, init_structure, is_disordered, join, lower, map, read_PIC, readlines, replace, rstrip, search, seek, set, set_anisou, set_coord, set_header, set_line_counter, set_parent, set_sigatm, set_siguij, set_symmetry, split, startswith, strip, unfold_entities, update, upper, values`

#### To extract `_write_SCAD`:
> You need these **52** components:
`_default_startpos, _generate_full_id, _get_azimuth, _get_startpos, _reset_full_id, _reset_properties, _translate_id, _writeSCAD_dihed, _write_SCAD, _write_mtx, add, altloc_match, assemble, assemble_residues, assemble_residues_ser, cast, clear_transforms, coord_space, disordered_get, disordered_get_list, get, get_chains, get_full_id, get_id, get_models, get_parent, get_resname, get_spherical_coordinates, has_id, init_atom_coords, internal_to_atom_coordinates, is_backbone, is_disordered, items, keys, multi_coord_space, multi_rot_Y, multi_rot_Z, norm, pop, pretty_str, propagate_changes, rak, set, set_Y_homog_rot_mtx, set_Z_homog_rot_mtx, set_accuracy_95, set_homog_trans_mtx, set_parent, split_akl, values, write`

#### To extract `save`:
> You need these **52** components:
`_format_mmcif_col, _generate_full_id, _get_label_asym_id, _requires_newline, _requires_quote, _reset_full_id, _revert_write, _save_dict, _save_structure, _translate_id, accept_atom, accept_chain, accept_model, accept_residue, add, close, copy, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, get_altloc, get_bfactor, get_coord, get_full_id, get_id, get_list, get_name, get_occupancy, get_parent, get_resname, get_serial_number, get_unpacked_list, has_id, index, is_disordered, items, replace, save, seek, set_coord, set_parent, split, startswith, strip, tell, truncate, values, write`

#### To extract `Main_Script`:
> You need these **103** components:
`Main_Script, _build_structure, _chop_end_codes, _chop_end_misc, _format_date, _generate_full_id, _get_biomoltrans, _get_header, _get_journal, _get_references, _handle_PDB_exception, _is_completely_disordered, _nice_case, _parse, _parse_atom_from, _parse_coordinates, _parse_header_from, _parse_pdb_header_list, _parse_remark_465, _parse_residue_id_from, _parse_resolution_from, _print_default_format_warning, _reset_full_id, _translate_id, _update_header_entry, add, as_handle, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, download_all_assemblies, download_entire_pdb, download_obsolete_entries, find, flag_disorder, flag_disordered, get, get_all_assemblies, get_all_entries, get_all_obsolete, get_altloc, get_coord, get_full_id, get_fullname, get_id, get_level, get_list, get_occupancy, get_parent, get_recent_changes, get_resname, get_status_list, get_structure, get_unpacked_list, has_id, index, init_atom, init_chain, init_model, init_residue, init_seg, init_structure, is_disordered, items, join, lower, map, namedtuple, read, readlines, replace, retrieve_assembly_file, retrieve_pdb_file, rstrip, search, seek, set_anisou, set_coord, set_header, set_line_counter, set_parent, set_sigatm, set_siguij, set_symmetry, split, startswith, strip, submit, unfold_entities, update, update_pdb, upper, urlcleanup, urlopen, urlretrieve, values, write, writelines`

#### To extract `write_PIC`:
> You need these **55** components:
`_enumerate_entity_atoms, _generate_full_id, _get_ak_tuple, _get_angle_for_tuple, _pdb_atom_string, _reset_full_id, _residue_string, _translate_id, _wpr, _write_PIC, _write_pic_bfac, add, as_handle, bits, cast, copy, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, enumerate_atoms, flag_disorder, get, get_altloc, get_atoms, get_bfactor, get_chains, get_coord, get_full_id, get_id, get_list, get_models, get_occupancy, get_parent, get_residues, get_resname, get_segid, get_serial_number, get_unpacked_list, has_id, is_disordered, pdb_date, pick_angle, rak, set_coord, set_parent, set_serial_number, split, startswith, upper, values, write, write_PIC`

#### To extract `internal_to_atom_coordinates`:
> You need these **44** components:
`_default_startpos, _generate_full_id, _get_azimuth, _get_startpos, _reset_full_id, _reset_properties, _translate_id, add, altloc_match, assemble, assemble_residues, assemble_residues_ser, cast, coord_space, disordered_get, disordered_get_list, get, get_chains, get_full_id, get_id, get_models, get_parent, get_resname, get_spherical_coordinates, has_id, init_atom_coords, internal_to_atom_coordinates, is_disordered, keys, multi_coord_space, multi_rot_Y, multi_rot_Z, norm, pop, pretty_str, propagate_changes, rak, set, set_Y_homog_rot_mtx, set_Z_homog_rot_mtx, set_homog_trans_mtx, set_parent, split_akl, values`

#### To extract `_create_edra`:
> You need these **29** components:
`_build_rak_cache, _create_edra, _gen_edra, _generate_full_id, _get_hedron, _link_dihedra, _reset_full_id, _reset_properties, _set_hedra, _translate_id, add, altloc_match, cast, disordered_get, disordered_get_list, get, get_full_id, get_id, get_parent, get_resname, has_id, is_disordered, items, rak, set, set_parent, split_akl, update, values`

#### To extract `disordered_add`:
> You need these **15** components:
`_generate_full_id, _reset_full_id, disordered_add, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, get_altloc, get_full_id, get_id, get_occupancy, get_parent, get_resname, set_parent, values`

#### To extract `write_SCAD`:
> You need these **77** components:
`_build_rak_cache, _create_edra, _default_startpos, _gen_edra, _generate_full_id, _get_azimuth, _get_hedron, _get_startpos, _link_dihedra, _reset_full_id, _reset_properties, _set_hedra, _spec_glyCB, _translate_id, _writeSCAD_dihed, _write_SCAD, _write_mtx, add, altloc_match, as_handle, assemble, assemble_residues, assemble_residues_ser, atom_to_internal_coordinates, build_atomArray, build_edraArrays, cast, clear_transforms, coord_space, copy, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, get, get_altloc, get_chains, get_coord, get_full_id, get_id, get_models, get_occupancy, get_parent, get_resname, get_spherical_coordinates, has_id, homog_scale_mtx, init_atom_coords, init_edra, internal_to_atom_coordinates, is_backbone, is_disordered, items, keys, multi_coord_space, multi_rot_Y, multi_rot_Z, norm, pop, pretty_str, propagate_changes, rak, set, set_Y_homog_rot_mtx, set_Z_homog_rot_mtx, set_accuracy_95, set_coord, set_homog_trans_mtx, set_parent, split_akl, update, values, write, write_SCAD`

#### To extract `init_residue`:
> You need these **29** components:
`_generate_full_id, _is_completely_disordered, _reset_full_id, _translate_id, add, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, flag_disorder, get_altloc, get_coord, get_full_id, get_id, get_list, get_occupancy, get_parent, get_resname, get_unpacked_list, has_id, init_residue, is_disordered, set_coord, set_parent, values`

#### To extract `structure_rebuild_test`:
> You need these **149** components:
`_atmfid_d2h, _build_rak_cache, _build_structure, _chop_end_codes, _chop_end_misc, _cmp_atm, _cmp_res, _create_edra, _default_startpos, _enumerate_entity_atoms, _format_date, _gen_edra, _generate_full_id, _get_ak_tuple, _get_angle_for_tuple, _get_azimuth, _get_biomoltrans, _get_header, _get_hedron, _get_journal, _get_references, _get_startpos, _handle_PDB_exception, _is_completely_disordered, _link_dihedra, _nice_case, _parse, _parse_atom_from, _parse_coordinates, _parse_header_from, _parse_pdb_header_list, _parse_remark_465, _parse_residue_id_from, _parse_resolution_from, _pdb_atom_string, _reset_full_id, _reset_properties, _residue_string, _set_hedra, _spec_glyCB, _translate_id, _update_header_entry, _wpr, _write_PIC, _write_pic_bfac, add, altloc_match, as_handle, assemble, assemble_residues, assemble_residues_ser, atom_to_internal_coordinates, bits, build_atomArray, build_edraArrays, cast, compare_residues, coord_space, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, enumerate_atoms, find, flag_disorder, flag_disordered, get, get_altloc, get_atoms, get_bfactor, get_chains, get_coord, get_full_id, get_fullname, get_id, get_level, get_list, get_models, get_occupancy, get_parent, get_residues, get_resname, get_segid, get_serial_number, get_spherical_coordinates, get_structure, get_unpacked_list, has_id, index, init_atom, init_atom_coords, init_chain, init_edra, init_model, init_residue, init_seg, init_structure, internal_to_atom_coordinates, is_disordered, items, join, keys, lower, map, multi_coord_space, multi_rot_Y, multi_rot_Z, norm, pdb_date, pick_angle, pop, pretty_str, propagate_changes, rak, read_PIC, readlines, replace, report_IC, rstrip, search, seek, set, set_Y_homog_rot_mtx, set_Z_homog_rot_mtx, set_anisou, set_coord, set_header, set_homog_trans_mtx, set_line_counter, set_parent, set_serial_number, set_sigatm, set_siguij, set_symmetry, split, split_akl, startswith, strip, structure_rebuild_test, unfold_entities, update, upper, values, write, write_PIC`

#### To extract `write_PDB`:
> You need these **72** components:
`_enumerate_entity_atoms, _format_mmcif_col, _generate_full_id, _get_label_asym_id, _is_completely_disordered, _requires_newline, _requires_quote, _reset_full_id, _revert_write, _save_dict, _save_structure, _translate_id, accept_atom, accept_chain, accept_model, accept_residue, add, as_handle, close, copy, detach_child, detach_parent, disordered_add, disordered_get, disordered_get_list, disordered_has_id, disordered_select, enumerate_atoms, flag_disorder, get, get_altloc, get_atoms, get_bfactor, get_chains, get_coord, get_full_id, get_id, get_list, get_models, get_name, get_occupancy, get_parent, get_residues, get_resname, get_serial_number, get_unpacked_list, has_id, index, init_chain, init_model, init_residue, init_seg, init_structure, is_disordered, items, pdb_date, replace, save, seek, set_coord, set_parent, set_serial_number, set_structure, split, startswith, strip, tell, truncate, upper, values, write, write_PDB`

#### To extract `assemble`:
> You need these **32** components:
`_default_startpos, _generate_full_id, _get_azimuth, _get_startpos, _reset_full_id, _reset_properties, _translate_id, add, altloc_match, assemble, cast, coord_space, disordered_get, disordered_get_list, get, get_full_id, get_id, get_parent, get_resname, get_spherical_coordinates, has_id, is_disordered, norm, pop, rak, set, set_Y_homog_rot_mtx, set_Z_homog_rot_mtx, set_homog_trans_mtx, set_parent, split_akl, values`

## 📑 Top-Level API Contents & Logic Flow
### 🔧 Functions
#### `calc_angle(v1, v2, v3)`
> Calculate angle method.
<details><summary>Full Docstring</summary>

```text
Calculate angle method.

Calculate the angle between 3 vectors
representing 3 connected points.

:param v1, v2, v3: the tree points that define the angle
:type v1, v2, v3: L{Vector}

:return: angle
:rtype: float
```
</details>

#### `calc_dihedral(v1, v2, v3, v4)`
> Calculate dihedral angle method.
<details><summary>Full Docstring</summary>

```text
Calculate dihedral angle method.

Calculate the dihedral angle between 4 vectors
representing 4 connected points. The angle is in
]-pi, pi].

:param v1, v2, v3, v4: the four points that define the dihedral angle
:type v1, v2, v3, v4: L{Vector}
```
</details>

#### `extract(structure, chain_id, start, end, filename)`
> Write out selected portion to filename.
<details><summary>Full Docstring</summary>

```text
Write out selected portion to filename.
```
</details>

#### `get_surface(model, MSMS='msms')`
> Represent molecular surface as a vertex list array.
<details><summary>Full Docstring</summary>

```text
Represent molecular surface as a vertex list array.

Return a NumPy array that represents the vertex list of the
molecular surface.

Arguments:
 - model - BioPython PDB model object (used to get atoms for input model)
 - MSMS - msms executable (used as argument to subprocess.call)
```
</details>

#### `is_aa(residue, standard=False)`
> Return True if residue object/string is an amino acid.
<details><summary>Full Docstring</summary>

```text
Return True if residue object/string is an amino acid.

:param residue: a L{Residue} object OR a three letter amino acid code
:type residue: L{Residue} or string

:param standard: flag to check for the 20 AA (default false)
:type standard: boolean

>>> is_aa('ALA')
True

Known three letter codes for modified amino acids are supported,

>>> is_aa('FME')
True
>>> is_aa('FME', standard=True)
False
```
</details>

#### `is_nucleic(residue, standard=False)`
> Return True if residue object/string is a nucleic acid.
<details><summary>Full Docstring</summary>

```text
Return True if residue object/string is a nucleic acid.

:param residue: a L{Residue} object OR a three letter code
:type residue: L{Residue} or string

:param standard: flag to check for the 8 (DNA + RNA) canonical bases.
    Default is False.
:type standard: boolean

>>> is_nucleic('DA ')
True

>>> is_nucleic('A  ')
True

Known three letter codes for modified nucleotides are supported,

>>> is_nucleic('A2L')
True
>>> is_nucleic('A2L', standard=True)
False
```
</details>

#### `m2rotaxis(m)`
> Return angles, axis pair that corresponds to rotation matrix m.
<details><summary>Full Docstring</summary>

```text
Return angles, axis pair that corresponds to rotation matrix m.

The case where ``m`` is the identity matrix corresponds to a singularity
where any rotation axis is valid. In that case, ``Vector([1, 0, 0])``,
is returned.
```
</details>

#### `make_dssp_dict(filename)`
> DSSP dictionary mapping identifiers to properties.
<details><summary>Full Docstring</summary>

```text
DSSP dictionary mapping identifiers to properties.

Return a DSSP dictionary that maps (chainid, resid) to
aa, ss and accessibility, from a DSSP file.

Parameters
----------
filename : string
    the DSSP output file
```
</details>

#### `parse_pdb_header(infile)`
> Return the header lines of a pdb file as a dictionary.
<details><summary>Full Docstring</summary>

```text
Return the header lines of a pdb file as a dictionary.

Dictionary keys are: head, deposition_date, release_date, structure_method,
resolution, structure_reference, journal_reference, author and
compound.
```
</details>

#### `refmat(p, q)`
> Return a (left multiplying) matrix that mirrors p onto q.
<details><summary>Full Docstring</summary>

```text
Return a (left multiplying) matrix that mirrors p onto q.

:type p,q: L{Vector}
:return: The mirror operation, a 3x3 NumPy array.

Examples
--------
>>> from Bio.PDB.vectors import refmat
>>> p, q = Vector(1, 2, 3), Vector(2, 3, 5)
>>> mirror = refmat(p, q)
>>> qq = p.left_multiply(mirror)
>>> print(q)
<Vector 2.00, 3.00, 5.00>
>>> print(qq)
<Vector 1.21, 1.82, 3.03>
```
</details>

#### `rotaxis(theta, vector)`
> Calculate left multiplying rotation matrix.
<details><summary>Full Docstring</summary>

```text
Calculate left multiplying rotation matrix.

Calculate a left multiplying rotation matrix that rotates
theta rad around vector.

:type theta: float
:param theta: the rotation angle

:type vector: L{Vector}
:param vector: the rotation axis

:return: The rotation matrix, a 3x3 NumPy array.

Examples
--------
>>> from numpy import pi
>>> from Bio.PDB.vectors import rotaxis2m
>>> from Bio.PDB.vectors import Vector
>>> m = rotaxis2m(pi, Vector(1, 0, 0))
>>> Vector(1, 2, 3).left_multiply(m)
<Vector 1.00, -2.00, -3.00>
```
</details>

#### `rotaxis2m(theta, vector)`
> Calculate left multiplying rotation matrix.
<details><summary>Full Docstring</summary>

```text
Calculate left multiplying rotation matrix.

Calculate a left multiplying rotation matrix that rotates
theta rad around vector.

:type theta: float
:param theta: the rotation angle

:type vector: L{Vector}
:param vector: the rotation axis

:return: The rotation matrix, a 3x3 NumPy array.

Examples
--------
>>> from numpy import pi
>>> from Bio.PDB.vectors import rotaxis2m
>>> from Bio.PDB.vectors import Vector
>>> m = rotaxis2m(pi, Vector(1, 0, 0))
>>> Vector(1, 2, 3).left_multiply(m)
<Vector 1.00, -2.00, -3.00>
```
</details>

#### `rotmat(p, q)`
> Return a (left multiplying) matrix that rotates p onto q.
<details><summary>Full Docstring</summary>

```text
Return a (left multiplying) matrix that rotates p onto q.

:param p: moving vector
:type p: L{Vector}

:param q: fixed vector
:type q: L{Vector}

:return: rotation matrix that rotates p onto q
:rtype: 3x3 NumPy array

Examples
--------
>>> from Bio.PDB.vectors import rotmat
>>> p, q = Vector(1, 2, 3), Vector(2, 3, 5)
>>> r = rotmat(p, q)
>>> print(q)
<Vector 2.00, 3.00, 5.00>
>>> print(p)
<Vector 1.00, 2.00, 3.00>
>>> p.left_multiply(r)
<Vector 1.21, 1.82, 3.03>
```
</details>

#### `vector_to_axis(line, point)`
> Vector to axis method.
<details><summary>Full Docstring</summary>

```text
Vector to axis method.

Return the vector between a point and
the closest point on a line (ie. the perpendicular
projection of the point on the line).

:type line: L{Vector}
:param line: vector defining a line

:type point: L{Vector}
:param point: vector defining the point
```
</details>

### 📦 Classes
#### `class CEAligner(window_size=8, max_gap=30)`
Protein Structure Alignment by Combinatorial Extension.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, window_size=8, max_gap=30)` | Superimpose one set of atoms onto another using structural data. |
| **align** | `(self, structure, transform=True, *, final_optimization=True)` | Align the input structure onto the reference structure. |
| **get_guide_coord_from_structure** | `(self, structure)` | Return the coordinates of guide atoms in the structure. |
| **set_reference** | `(self, structure)` | Define a reference structure onto which all others will be aligned. |


#### `class CaPPBuilder(radius=4.3)`
Use CA--CA distance to find polypeptides.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, radius=4.3)` | Initialize the class. |
| **build_peptides** | `(self, entity, aa_only=1)` | Build and return a list of Polypeptide objects. |


#### `class DSSP(model, in_file, dssp='dssp', acc_array='Sander', file_type='')`
Run DSSP and parse secondary structure and accessibility.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, model, in_file, dssp='dssp', acc_array='Sander', file_type='')` | Create a DSSP object. |
| **keys** | `(self)` | Return the list of residues. |


#### `class ExposureCN(model, radius=12.0, offset=0)`
Residue exposure as number of CA atoms around its CA atom.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, model, radius=12.0, offset=0)` | Initialize class. |
| **keys** | `(self)` | Return the list of residues. |


#### `class FastMMCIFParser(structure_builder=None, auth_chains=True, auth_residues=True, QUIET=False)`
Parse an MMCIF file and return a Structure object.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, structure_builder=None, auth_chains=True, auth_residues=True, QUIET=False)` | Create a FastMMCIFParser object. |
| **get_structure** | `(self, structure_id, filename)` | Return the structure. |


#### `class FragmentMapper(model, lsize=20, flength=5, fdir='.')`
Map polypeptides in a model to lists of representative fragments.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, model, lsize=20, flength=5, fdir='.')` | Create instance of FragmentMapper. |


#### `class HSExposureCA(model, radius=12, offset=0)`
Class to calculate HSE based on the approximate CA-CB vectors.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, model, radius=12, offset=0)` | Initialize class. |
| **keys** | `(self)` | Return the list of residues. |
| **pcb_vectors_pymol** | `(self, filename='hs_exp.py')` | Write PyMol script for visualization. |


#### `class HSExposureCB(model, radius=12, offset=0)`
Class to calculate HSE based on the real CA-CB vectors.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, model, radius=12, offset=0)` | Initialize class. |
| **keys** | `(self)` | Return the list of residues. |


#### `class MMCIFIO()`
Write a Structure object or a mmCIF dictionary as a mmCIF file.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self)` | Initialise. |
| **save** | `(self, filepath, select=<Select all>, preserve_atom_numbering=False)` | Save the structure to a file. |
| **set_dict** | `(self, dic)` | Set the mmCIF dictionary to be written out. |
| **set_structure** | `(self, pdb_object)` | Check what the user is providing and build a structure. |


#### `class MMCIFParser(structure_builder=None, auth_chains=True, auth_residues=True, QUIET=False)`
Parse a mmCIF file and return a Structure object.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, structure_builder=None, auth_chains=True, auth_residues=True, QUIET=False)` | Create a PDBParser object. |
| **get_structure** | `(self, structure_id, filename)` | Return the structure. |


#### `class NeighborSearch(atom_list, bucket_size=10)`
Class for neighbor searching.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, atom_list, bucket_size=10)` | Create the object. |
| **search** | `(self, center, radius, level='A')` | Neighbor search. |
| **search_all** | `(self, radius, level='A')` | All neighbor search. |


#### `class PDBIO(use_model_flag=0, is_pqr=False)`
Write a Structure object (or a subset of a Structure object) as a PDB or PQR file.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, use_model_flag=0, is_pqr=False)` | Create the PDBIO object. |
| **save** | `(self, file, select=<Select all>, write_end=True, preserve_atom_numbering=False)` | Save structure to a file. |
| **set_structure** | `(self, pdb_object)` | Check what the user is providing and build a structure. |


#### `class PDBList(server='https://files.wwpdb.org', pdb=None, obsolete_pdb=None, verbose=True)`
Quick access to the structure lists on the PDB or its mirrors.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, server='https://files.wwpdb.org', pdb=None, obsolete_pdb=None, verbose=True)` | Initialize the class with the default server or a custom one. |
| **download_all_assemblies** | `(self, listfile: Optional[str] = None, file_format: Optional[str] = None, max_num_threads: Optional[int] = None)` | Retrieve all biological assemblies not in the local PDB copy. |
| **download_entire_pdb** | `(self, listfile: Optional[str] = None, file_format: Optional[str] = None, max_num_threads: Optional[int] = None)` | Retrieve all PDB entries not present in the local PDB copy. |
| **download_obsolete_entries** | `(self, listfile: Optional[str] = None, file_format: Optional[str] = None, max_num_threads: Optional[int] = None)` | Retrieve all obsolete PDB entries not present in local obsolete PDB copy. |
| **download_pdb_files** | `(self, pdb_codes: list[str], obsolete: bool = False, pdir: Optional[str] = None, file_format: Optional[str] = None, overwrite: bool = False, max_num_threads: Optional[int] = None)` | Fetch set of PDB structure files from the PDB server and store them locally. |
| **get_all_assemblies** | `(self, file_format: str = '') -> list[tuple[str, str]]` | Retrieve the list of PDB entries with an associated bio assembly. |
| **get_all_entries** | `(self)` | Retrieve the big file containing all the PDB entries and some annotation. |
| **get_all_obsolete** | `(self)` | Return a list of all obsolete entries ever in the PDB. |
| **get_recent_changes** | `(self)` | Return three lists of the newest weekly files (added,mod,obsolete). |
| **get_seqres_file** | `(self, savefile='pdb_seqres.txt')` | Retrieve and save a (big) file containing all the sequences of PDB entries. |
| **get_status_list** | `(url)` | Retrieve a list of pdb codes in the weekly pdb status file from given URL. |
| **retrieve_assembly_file** | `(self, pdb_code, assembly_num, pdir=None, file_format=None, overwrite=False)` | Fetch one or more assembly structures associated with a PDB entry. |
| **retrieve_pdb_file** | `(self, pdb_code, obsolete=False, pdir=None, file_format=None, overwrite=False)` | Fetch PDB structure file from PDB server, and store it locally. |
| **update_pdb** | `(self, file_format=None, with_assemblies=False)` | Update your local copy of the PDB files. |


#### `class PDBMLParser()`
A parser for PDBML (PDB XML) files. See https://pdbml.wwpdb.org/.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self)` | Initialize a PDBML parser. |
| **get_structure** | `(self, source: Union[int, str, bytes, os.PathLike, TextIO]) -> Bio.PDB.Structure.Structure` | Parse and return the PDB structure from XML source. |


#### `class PDBParser(PERMISSIVE=True, get_header=False, structure_builder=None, QUIET=False, is_pqr=False)`
Parse a PDB file and return a Structure object.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, PERMISSIVE=True, get_header=False, structure_builder=None, QUIET=False, is_pqr=False)` | Create a PDBParser object. |
| **get_header** | `(self)` | Return the header. |
| **get_structure** | `(self, id, file)` | Return the structure. |
| **get_trailer** | `(self)` | Return the trailer. |


#### `class PPBuilder(radius=1.8)`
Use C--N distance to find polypeptides.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, radius=1.8)` | Initialize the class. |
| **build_peptides** | `(self, entity, aa_only=1)` | Build and return a list of Polypeptide objects. |


#### `class ResidueDepth(model, msms_exec=None)`
Calculate residue and CA depth for all residues.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, model, msms_exec=None)` | Initialize the class. |
| **keys** | `(self)` | Return the list of residues. |


#### `class Select()`
Select everything for PDB output (for use as a base class).

| Method | Signature | Description |
| :--- | :--- | :--- |
| **accept_atom** | `(self, atom)` | Overload this to reject atoms for output. |
| **accept_chain** | `(self, chain)` | Overload this to reject chains for output. |
| **accept_model** | `(self, model)` | Overload this to reject models for output. |
| **accept_residue** | `(self, residue)` | Overload this to reject residues for output. |


#### `class ShrakeRupley(probe_radius=1.4, n_points=100, radii_dict=None)`
Calculates SASAs using the Shrake-Rupley algorithm.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, probe_radius=1.4, n_points=100, radii_dict=None)` | Initialize the class. |
| **compute** | `(self, entity, level='A')` | Calculate surface accessibility surface area for an entity. |


#### `class StructureAlignment(fasta_align, m1, m2, si=0, sj=1)`
Class to align two structures based on an alignment of their sequences.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, fasta_align, m1, m2, si=0, sj=1)` | Initialize. |
| **get_iterator** | `(self)` | Create an iterator over all residue pairs. |
| **get_maps** | `(self)` | Map residues between the structures. |


#### `class Superimposer()`
Rotate/translate one set of atoms on top of another to minimize RMSD.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self)` | Initialize the class. |
| **apply** | `(self, atom_list)` | Rotate/translate a list of atoms. |
| **set_atoms** | `(self, fixed, moving)` | Prepare translation/rotation to minimize RMSD between atoms. |


#### `class Vector(x, y=None, z=None)`
3D vector.

| Method | Signature | Description |
| :--- | :--- | :--- |
| **__init__** | `(self, x, y=None, z=None)` | Initialize the class. |
| **angle** | `(self, other)` | Return angle between two vectors. |
| **copy** | `(self)` | Return a deep copy of the Vector. |
| **get_array** | `(self)` | Return (a copy of) the array of coordinates. |
| **left_multiply** | `(self, matrix)` | Return Vector=Matrix x Vector. |
| **norm** | `(self)` | Return vector norm. |
| **normalize** | `(self)` | Normalize the Vector object. |
| **normalized** | `(self)` | Return a normalized copy of the Vector. |
| **normsq** | `(self)` | Return square of vector norm. |
| **right_multiply** | `(self, matrix)` | Return Vector=Vector x Matrix. |

