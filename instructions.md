# mzXplorer Instructions

mzXplorer is an interactive Shiny application for exploring mass spectrometric feature tables through **mass defect visualisation**, **Kendrick-style projections**, and an advanced **homologue‑series detection engine**. It enables rapid interpretation of complex HRMS data from LC/GC-HRMS, FIA, direct infusion, or DIA workflows.

------------------------------------------------------------------------

# 1. Input Data

## 1.1 File format

mzXplorer accepts a single input file in **CSV format (.csv)**.

Your file must contain **three mandatory columns** (with exact names):

| Column | Description |
|----------------------------|--------------------------------------------|
| `mz` | Mass‑to‑charge ratio of the detected feature |
| `rt` | Retention time (in minutes or seconds, consistent across dataset) |
| `intensity` | Peak area or intensity |

If your experiment has **no retention time** (e.g., direct infusion), simply set the `rt` column to a constant such as `1`.\
Column names must not contain spaces.

You may include **additional user‑defined columns**, which will appear in:

-   the data table\
-   variable selectors\
-   exported output

See 1.2 for details.

------------------------------------------------------------------------

## 1.2 Additional user-defined variables

Your CSV may include any number of extra columns. These must follow R naming rules:

-   start with a letter\
-   contain only letters, digits, `_` (underscore), or `.` (dot)\
-   no spaces

Useful examples:

-   `chemical_formula`\
-   `chemical_group` (PFAS, PBDEs, PAHs, lipids, steroids, …)\
-   `component` (DIA deconvolution groups)\
-   `sample_id`\
-   `feature_id`

Short names keep the table readable.

------------------------------------------------------------------------

# 2. Mass Defect Bases and Derived Variables

After uploading your data, enter one or more **mass defect (MD) base units**.\
mzXplorer supports **up to three‑order MD calculations**.

## 2.1 Valid MD base formats

| Input      | Meaning                               |
|------------|---------------------------------------|
| `CH2`      | First‑order MD base (14.01565 Da)     |
| `Cl-H`     | Mass of Cl minus mass of H            |
| `CH2,O`    | CH₂ as first order, O as second order |
| `CH2,Cl-H` | CH₂ (1st) + (Cl − H) (2nd)            |
| `CH2/10`   | Fractional unit (CH₂ / 10)            |

Rules:

-   Separate units by commas (no spaces).\
-   Fractional units cannot be combined with subtraction units.\
-   Maximum three units (MD1, MD2, MD3).

After clicking **Plot**, mzXplorer computes:

-   OMD (Original Mass Defect)\
-   RMD (Relative Mass Defect)\
-   MD1, MD2, MD3 (if applicable)

All these become selectable for plotting.

------------------------------------------------------------------------

# 3. Interactive Plots

mzXplorer displays **two synchronized Plotly scatter plots**.

## 3.1 Axis controls

You can choose any numeric variable as: - X‑axis\
- Y‑axis

for each plot independently.

## 3.2 Intensity scaling (optional)

Checkbox: **Show intensity as size**

-   When enabled → point size is proportional to selected intensity variable\
-   When disabled → constant point size\
-   Applies to both scatter plots and the barplot

## 3.3 Colour behaviour

-   Background points are **always black**\
-   Selected homologue series are drawn with **distinct colours**\
-   When one or more series are selected, background points automatically dim (opacity lowered from 0.35 → 0.10)

## 3.4 Crosstalk selection

Selecting points in one plot:

-   highlights them in the other plot\
-   updates the selected‑data table\
-   updates the barplot

Plots support zooming, panning, and lasso/box selection.

------------------------------------------------------------------------

# 4. Filters

mzXplorer generates three dynamic filters:

-   Intensity range\
-   m/z range\
-   Retention time range

Filters apply to:

-   both plots\
-   homologue detection\
-   the data table\
-   exported output

Click **Plot** to apply changes.

------------------------------------------------------------------------

# 5. Data Table and Export

## 5.1 Selected-data table

Located below the plots:

-   If points are selected → table shows **only selected rows**\
-   If no selection → table shows **the full dataset**

The table supports filtering, sorting, and column reordering.

## 5.2 Export

**Export Data** downloads:

-   selected rows (if any selected)\
-   full dataset otherwise

Export includes all derived MD variables and homologue-series IDs.

------------------------------------------------------------------------

# 6. Homologue Series Detection

mzXplorer includes a fast, graph‑based algorithm for detecting homologue series based on repeating mass differences.

Enable via:

> **Detect homologues in filtered data**

## 6.1 Parameters

-   **Repeating unit** (formula or exact mass)\
-   **ppm tolerance**\
-   **Minimum series length**\
-   **RT tolerance** between successive members\
-   **RT trend** (increasing / decreasing / any)\
-   **Allow gaps** (permit k=2 spacing)\
-   **Spline R² filter** (optional chromatographic smoothness test)

------------------------------------------------------------------------

# 7. How Homologue Detection Works (Function-Level Explanation)

The homologue detection system is powered by several core functions that work together:

------------------------------------------------------------------------

## 7.1 `build_edges()` — *Finding possible homologous relationships*

This function identifies all peak‑to‑peak transitions that could represent **one repeating‑unit increment** (or k = 2 when gaps are allowed).

### What it does:

1.  For each peak *i*, compute the expected m/z of the next peak:

$$
mz_i + k \times \text{unit mass}
$$

where $k = 1$ (normal step) or $k = 2$ (allowed gap).

2.  Convert ppm tolerance into exact m/z limits:

$$
[m_{\text{min}}, m_{\text{max}}]
$$

3.  Use fast binary search (`findInterval`) to find all peaks whose m/z falls inside this window.

4.  Filter edges according to:

-   RT trend (increasing / decreasing / any)\
-   RT tolerance\
-   allowed gap setting (k = 1 or k = 2)

### Output:

A table of candidate edges:

| from | to  |
|------|:----|
| 12   | 18  |
| 18   | 26  |
| 26   | 34  |

------------------------------------------------------------------------

## 7.2 `build_graph()` — *Constructing the homologous network*

This function takes the edges from `build_edges()` and constructs a **directed graph** using igraph.

**What it does:**

1.  Each peak becomes a **node**.\
2.  Each valid pair (from → to) becomes a **directed edge**.\
3.  Duplicate edges and self‑loops are removed.\
4.  Peaks filtered out earlier are excluded.

**Why a graph?**\
A homologue series is simply a **connected component** in this graph where masses increase stepwise.

------------------------------------------------------------------------

## 7.3 Component extraction (`igraph::decompose()`)

Once the graph is constructed:

-   Each connected subgraph is interpreted as **one homologue series**
-   Series are assigned unique `series_id` values

------------------------------------------------------------------------

## 7.4 Additional filtering functions

### **`strict_rt_filter()`**

Ensures that peaks within each candidate homologue series follow the selected retention‑time trend (“increasing” or “decreasing”). Any peak that violates the monotonic trend is removed from the series.

------------------------------------------------------------------------

### **`apply_shiny_splines()` — enforcing smooth chromatographic behaviour (optional)**

This function evaluates how well each detected homologue series follows a **smooth retention‑time progression** across increasing m/z.\
This is important in LC/GC–HRMS data, where homologous compounds (PFAS, PBDEs, lipids, etc.) typically elute in a predictable RT order.

#### **Why this matters**

In real datasets, false-positive series may pass mass-based matching but show: - retention-time jumps, - unstable elution patterns, - or noisy peak shapes. `apply_shiny_splines()` provides a way to **remove such non‑chromatographic series**.

------------------------------------------------------------------------

### **How it works**

For each homologue series:

1.  The peaks in that series are sorted by m/z.

2.  A smooth spline is fitted: $$
    \text{RT} = f(\text{m/z})
    $$ using `stats::smooth.spline()` with a moderate smoothing parameter (`spar = 0.45` in this implementation).

3.  The spline is used to generate predicted RT values.

4.  The correlation between observed RT and spline‑predicted RT is computed: $$
    R^2 = \left(\text{cor}(\text{RT}, \text{RT}_{\text{pred}})\right)^2
    $$

5.  If:

    -   the spline fit is poor\
    -   or $R^2$ is below the user-defined threshold (e.g. **0.98**)

    → the **entire series is discarded**.

------------------------------------------------------------------------

### **Interpretation**

-   **High R² (\~1.0):**\
    Series behaves like a real chromatographic family\
    (e.g., PFAS homologues, fatty acids, alkylated PAHs).
-   **Low R²:**\
    Series is likely noise, co-elution artefacts, or accidental mass matches.

**This step dramatically improves the quality of detected series** in LC/GC–HRMS applications.

------------------------------------------------------------------------

### **When to use it**

Enable spline filtering when:

-   you expect clean monotonic chromatographic behaviour,\
-   the dataset contains many false matches,\
-   retention time increases consistently with m/z.

Disable it (set R² threshold = 0) when:

-   analysing direct-infusion or FIA data (RT is non-informative),\
-   chromatographic behaviour is irregular (short gradients, HILIC plateaus, etc.).

------------------------------------------------------------------------

# 8. Visualisation of Homologue Series

Selecting one or more series in the homologue table:

-   overlays coloured line+marker traces in both scatter plots\
-   displays series members at full opacity\
-   dims background points (opacity = 0.10)\
-   updates table and barplot

A **Clear selection** button resets the view.

------------------------------------------------------------------------

# 9. Barplot of Selected Points

The barplot displays:

-   m/z (x-axis)\
-   relative intensity (y-axis, %)

The variable used is selected above the barplot.

It updates whenever the variable selection changes.

------------------------------------------------------------------------

# 10. Known Issues

## 10.1 Crosstalk table editing

Editing values directly in the filtered table may break selection sync between table and plots.\
Workaround:

-   Export the filtered table\
-   Re-import it for fresh analysis

## 10.2 Duplicate m/z–rt values

If duplicate pairs exist, homologue detection may assign ambiguous series IDs.\
Include a unique identifier column if needed.

------------------------------------------------------------------------

# 11. Mass Defect Equations

### Original Mass Defect (OMD)

$$
\text{OMD} = \text{round}(m) - m
$$

### Relative Mass Defect (RMD)

$$
\text{RMD} = \frac{\text{round}(m) - m}{m} \times 10^6
$$

### First‑order MD

$$
m_1 = m \times \frac{\text{round}(u_1)}{u_1}
$$ $$
\text{MD1} = \text{round}(m_1) - m_1
$$

### Second‑order MD

$$
m_2 = \frac{\text{MD1}(m)}{\text{MD1}(u_2)}
$$ $$
\text{MD2} = \text{round}(m_2) - m_2
$$

### Third‑order MD

$$
m_3 = \frac{\text{MD2}(m)}{\text{MD2}(u_3)}
$$ $$
\text{MD3} = \text{round}(m_3) - m_3
$$

------------------------------------------------------------------------

# 12. Tips for Effective Use

-   Use ppm = 1–5 for high‑resolution instruments.\
-   For LC/GC‑HRMS, set RT trend to **increasing**.\
-   Enable **Allow gaps** for incomplete series.\
-   Typical MD bases:
    -   CH2 → lipids\
    -   CF2 → PFAS\
    -   Cl-H → chlorinated compounds\
-   Experiment with multiple MD bases to reveal structural classes.
