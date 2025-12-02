# Documentation for `alphaPredict`
**File Path:** `/home/nicai_zht/.local/share/mamba/envs/idr_llm/lib/python3.13/site-packages/alphaPredict/__init__.py`

## Module Docstring
```text
short description of alphaPredict.

Predict confidence scores of alphaFold2
```

## 🚦 Navigator: How to Drive
This section helps you understand how to run this library from the command line or entry points.
- ℹ️ **No Direct Entry Point**: This module seems to be a library intended for import, not direct execution.

### 🐍 Python API Usage (Inferred)
Since no CLI entry point was found, here are the likely **Python API entry points** for your script:

#### 🚀 Top 20 Recommended Entry Points
| Type | API | Description |
| :--- | :--- | :--- |
| `ƒ` | **alphaPredict.predict**(**sequence**) | Function to return confidence scores from |
| `ƒ` | **alphaPredict.graph**(**sequence**, title, confidence_threshold, shaded_regions, shaded_region_color, confidence_line_color, threshold_line_color, DPI, output_file) | No description. |

> **Note:** Bold parameters are required. Others are optional.

#### 🧩 Code Snippets (Auto-Generated)
```python
import alphaPredict

# --- Top 20 Ranked Functions ---
# 1. predict
result_1 = alphaPredict.predict(sequence=...)

# 2. graph
result_2 = alphaPredict.graph(sequence=...)

```

_No explicit `argparse` configuration detected in the main module._


## 📊 Network & Architecture Analysis
### 🌍 Top 20 External Dependencies
| Library | Usage Count |
| :--- | :--- |
| **_frozen_importlib_external** | 4 |
| **_frozen_importlib** | 4 |


### 🕸️ Network Metrics (Advanced)
#### 👑 Top 20 Modules by PageRank (Authority)
| Rank | Module | Score | Type | Role |
| :--- | :--- | :--- | :--- | :--- |
| 1 | `_frozen_importlib_external` | 0.2067 | External | External Lib |
| 2 | `_frozen_importlib` | 0.2067 | External | External Lib |
| 3 | `alpha` | 0.1085 | Internal | Utility / Core |
| 4 | `alpha_exceptions` | 0.1030 | Internal | Utility / Core |
| 5 | `backend.parrot_alpha` | 0.1030 | Internal | External Lib |
| 6 | `backend.alpha_graph` | 0.1030 | Internal | External Lib |
| 7 | `alphaPredict` | 0.0846 | Internal | Utility / Core |
| 8 | `_version` | 0.0846 | Internal | Utility / Core |


### 🗺️ Dependency & Architecture Map
```mermaid
graph TD
    classDef core fill:#f96,stroke:#333,stroke-width:2px;
    classDef external fill:#9cf,stroke:#333,stroke-width:1px;
    id_6["alphaPredict"] -.-> id_1["_frozen_importlib_external"]
    class id_6 core;
    class id_1 external;
    id_6["alphaPredict"] -.-> id_4["_frozen_importlib"]
    class id_6 core;
    class id_4 external;
    id_6["alphaPredict"] --> id_7["alpha"]
    class id_6 core;
    class id_7 core;
    id_7["alpha"] --> id_10["alpha_exceptions"]
    class id_7 core;
    class id_10 core;
    id_7["alpha"] -.-> id_1["_frozen_importlib_external"]
    class id_7 core;
    class id_1 external;
    id_7["alpha"] -.-> id_4["_frozen_importlib"]
    class id_7 core;
    class id_4 external;
    id_7["alpha"] --> id_5["parrot_alpha"]
    class id_7 core;
    class id_5 core;
    id_7["alpha"] --> id_9["alpha_graph"]
    class id_7 core;
    class id_9 core;
    id_2["_version"] -.-> id_1["_frozen_importlib_external"]
    class id_2 core;
    class id_1 external;
    id_2["_version"] -.-> id_4["_frozen_importlib"]
    class id_2 core;
    class id_4 external;
    id_10["alpha_exceptions"] -.-> id_1["_frozen_importlib_external"]
    class id_10 core;
    class id_1 external;
    id_10["alpha_exceptions"] -.-> id_4["_frozen_importlib"]
    class id_10 core;
    class id_4 external;
    id_8["AlphaError"] ==> id_0["Exception"]
    class id_8 core;
    class id_0 external;
    id_8["AlphaError"] ==> id_0["Exception"]
    class id_8 core;
    class id_0 external;
    id_3["DomainError"] ==> id_0["Exception"]
    class id_3 core;
    class id_0 external;
```

## 🚀 Global Execution Flow & Extraction Guide
This graph visualizes how data flows between functions across the entire project.
```mermaid
graph TD
    classDef main fill:#f9f,stroke:#333,stroke-width:2px;
    classDef func fill:#fff,stroke:#333,stroke-width:1px;
    f_3["predict"] -->|sequence| f_0["_alpha_predict"]
    class f_3 func;
    class f_0 func;
    f_2["graph"] -->|sequence<br>title<br>confidence_threshold<br>shaded_regions<br>shaded_region_color<br>confidence_line_color<br>threshold_line_color<br>DPI<br>output_file| f_1["_graph"]
    class f_2 func;
    class f_1 func;
```

### ✂️ Navigator: Snippet Extractor
Want to use a specific function without the whole library? Here is the **Dependency Closure** for **Top 20** key functions.
#### To extract `predict`:
> You need these **2** components:
`_alpha_predict, predict`

#### To extract `graph`:
> You need these **2** components:
`_graph, graph`

## 📑 Top-Level API Contents & Logic Flow
### 🔧 Functions
#### `graph(sequence, title='Predicted Confidence Score', confidence_threshold=50, shaded_regions=None, shaded_region_color='red', confidence_line_color='blue', threshold_line_color='black', DPI=150, output_file=None)`
> No documentation available.
<details><summary>Full Docstring</summary>

```text
No documentation available.
```
</details>

#### `predict(sequence)`
> Function to return confidence scores from
<details><summary>Full Docstring</summary>

```text
Function to return confidence scores from
Alpha fold 2 of a single input sequence. Returns the
predicted values as a float.

Parameters
------------

sequence : str 
    Input amino acid sequence (as string) to be predicted.

Returns
--------

Float
    Returns a float of the confidence score value (predicted)
```
</details>
