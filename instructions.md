# mzXplorer Instructions

*mzXplorer* is an interactive Shiny application for exploring high‑resolution MS data using\
**mass‑defect visualisation**, and an advanced **RT‑ and CCS‑aware homologue‑series detection**.

It is designed for LC/GC‑HRMS, DI/FIA, DIA, and **IM‑MS** workflows.

------------------------------------------------------------------------

# 1. Input Data

## 1.1 Required columns

mzXplorer accepts a single **CSV (.csv)** file with at least:

| Column      | Description                 |
|-------------|-----------------------------|
| `mz`        | Mass‑to‑charge ratio        |
| `rt`        | Retention time (min or sec) |
| `intensity` | Peak height or area         |

-   `rt` may be constant (e.g. direct infusion).
-   Column names are automatically lower‑cased and trimmed.

------------------------------------------------------------------------

## 1.2 Optional columns

Any additional columns are allowed and become available for:

-   plotting\
-   filtering\
-   homologue summary\
-   export

Examples:

-   `ccs` *(collision cross section — fully supported)*\
-   `class`, `group`\
-   `chemical_formula`\
-   `sample_id`, `feature_id`

Names must follow standard R variable rules (letters, digits, `_`, `.`).

------------------------------------------------------------------------

# 2. Mass Defect Units & Derived Variables

mzXplorer supports up to **three mass‑defect bases**, entered as:

| Input      | Meaning                      |
|------------|------------------------------|
| `CH2`      | CH₂ mass                     |
| `Cl-H`     | replacement of H with Cl     |
| `CH2,O`    | first = CH₂, second = O      |
| `CH2,Cl-H` | CH₂ and (Cl−H)               |
| `CH2/10`   | fractional mass with base 10 |

After clicking **Plot**, mzXplorer computes:

-   **OMD** – original mass defect\
-   **RMD** – relative mass defect\
-   **MD1 (MD2**, **MD3)** – user input mass defects

All appear in axis selectors.

------------------------------------------------------------------------

# 3. Interactive Plots

The app shows **two synchronized scatter plots** using Plotly.

## 3.1 Axis controls

Each plot has independent:

-   X variable\
-   Y variable

Any numeric column is selectable.

## 3.2 Intensity scaling

Enabling **Show intensity as size** scales point size by the selected variable.

## 3.3 Homologue colouring

Selecting homologue series:

-   overlays coloured line+marker traces\
-   dims background points\
-   enhances visibility of structural patterns

## 3.4 Lasso/box selection

Plot selection:

-   highlights points across both plots\
-   updates the selected‑data table\
-   updates the barplot

Selections are driven by **plotly::event_data()**.

------------------------------------------------------------------------

# 4. Filters

Three automatically‑generated filters:

-   intensity range\
-   m/z range\
-   retention time (rt) range

Filtering affects:

-   plots\
-   homologue detection\
-   summary table\
-   exported data

Click **Plot** to apply.

------------------------------------------------------------------------

# 5. Selected‑Data Table & Export

## 5.1 Selected‑data table

Shows:

-   only selected points (from plots or series), or\
-   the full dataset if no selection exists.

Supports sorting and column filtering.

## 5.2 Export

**Export Data** downloads:

-   selected rows if any exist\
-   otherwise the full dataset

The export includes all MD variables, homologue IDs, and CCS values if present.

------------------------------------------------------------------------

# 6. Homologue‑Series Detection (Updated)

Enable detection via:

> **Detect homologues in filtered data**

mzXplorer implements a graph‑based homologue finder with **RT** and **CCS** monotonicity.

## 6.1 Parameters

| Parameter | Description |
|----|----|
| Repeating unit | formula (e.g. CH2 or CF2) |
| ppm tolerance | mass error in ppm |
| Minimum length | least number of homologues to be considered a series |
| RT tolerance | max rt difference per step |
| RT trend | increasing / decreasing / any |
| Allow gaps | allows k = 2 jumps (allow skipping one repeating unit in the middle of a series) |
| Spline R² | RT smoothness filter |
| Enable CCS rules | activates CCS settings |
| CCS mode | RT only / CCS only / Both |
| CCS tolerance | max CCS difference per step |

------------------------------------------------------------------------

# 7. Homologue Detection Pipeline

The algorithm performs:

1.  **Candidate edges** via mass defect + optional RT restrictions (`build_edges()`).
2.  **Graph construction** (`build_graph()`).
3.  **Connected components** → provisional series.
4.  **Monotonicity filtering** (`strict_rt_filter()`):
    -   RT monotonicity\
    -   CCS monotonicity\
    -   CCS tolerance\
5.  **Minimum length** check.\
6.  **Chromatographic smoothness** (`apply_shiny_splines()`).
7.  Final renumbering of surviving series.

------------------------------------------------------------------------

# 8. CCS‑Aware Monotonicity

When CCS support is enabled, the user may choose:

### **RT only**

-   classical RT monotonicity
-   CCS ignored

### **CCS only**

-   CCS monotonicity enforced\
-   RT ignored\
-   optional CCS tolerance filter

### **Both RT + CCS**

A point must satisfy:

-   RT monotonicity *and*\
-   CCS monotonicity *and*, if enabled\
-   CCS tolerance per step

Ideal for LC‑IM‑HRMS workflows.

------------------------------------------------------------------------

# 9. Homologue Table (Updated)

Lists all homologue series passing all filters:

Includes:

-   `series_id`, `n`\
-   `mz_min`, `mz_max`\
-   `rt_min`, `rt_max`\
-   `int_sum`

### **New CCS summary columns (if available)**

-   `ccs_min`\
-   `ccs_max`\
-   `ccs_range`

### Selection behaviour

Clicking one or more rows:

-   colours those series in both plots\
-   dims other points\
-   sets *series selection* to **highest priority**\
-   updates **selected‑data table & barplot to show series members**

Use **Clear selection** to return to plot‑based selection.

------------------------------------------------------------------------

# 10. Selection Logic (New)

mzXplorer now uses **unified selection priority**:

## Priority:

1.  **Homologue‑table selection** (highest)
    -   drives table & barplot\
    -   drives colour overlays
2.  **Plot selection**
    -   active only when no series is selected

This ensures intuitive control:\
**series selection overrides plot selection**.

------------------------------------------------------------------------

# 11. Barplot

Barplot shows:

-   x = m/z\
-   y = relative (% of max) of selected variable

Only selected points appear.

Updates when:

-   selection changes (plot ↔ homologue table)\
-   intensity variable changes

------------------------------------------------------------------------

# 12. Known Issues

-   To be added

------------------------------------------------------------------------

# 13. Tips for Effective Use

-   LC/GC‑HRMS:
    -   RT trend = **increasing**\
    -   Enable **Both (RT + CCS)**\
    -   Use spline R² \> 0.95 for stricter monotonic trend
-   IM‑MS:
    -   Activate CCS mode\
    -   Tune CCS tolerance
-   PFAS:
    -   Use `CF2` or `CF2,SO2`
-   Lipids:
    -   Use `CH2`\
    -   Allow gaps for missing chain lengths
-   Direct infusion:
    -   RT trend = **any**\
    -   Spline R² = 0

------------------------------------------------------------------------

# 14. Mass Defect Equations

### OMD

$$
\text{OMD} = \text{round}(m) - m
$$

### RMD

$$
\text{RMD} = \frac{\text{round}(m) - m}{m} \times 10^6
$$

### MD1

$$
m_1 = m \times \frac{\text{round}(u_1)}{u_1}
$$ $$
\text{MD1} = \text{round}(m_1) - m_1
$$

### MD2

$$
m_2 = \frac{\text{MD1}(m)}{\text{MD1}(u_2)}
$$ $$
\text{MD2} = \text{round}(m_2) - m_2
$$

### MD3

$$
m_3 = \frac{\text{MD2}(m)}{\text{MD2}(u_3)}
$$ $$
\text{MD3} = \text{round}(m_3) - m_3
$$

Examples of chemical formula in the MD formula insert box:

$CH2$: corresponds to one MD base unit, in this case a methylene unit, corresponding to the exact mass of 14.01565.

$Cl-H$: calculates the addition of one Cl atom and subtraction of one H atom, corresponding to the exact mass of 33.96103. Use the minus sign “-“ to separate the base units. This app only support two different unit.

$CH2,O$: specifies that $CH_2$ is the first-order MD unit and $O$ is the second-order MD unit, use comma without blank space to separate them. This app only support at most three-order mass defect units.

$CH2,Cl-H$: specifies that $CH_2$ is the first-order MD unit and $Cl-H$ is the second-order MD unit.

Use comma without blank space to separate the units.

**Fractional base unit** can also be used (but not together with subtraction) as follows:

$CH2/X$: where X is the divisor. If X = 10 then the input formula will be: $CH2/10$.
