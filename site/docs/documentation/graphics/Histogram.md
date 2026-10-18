# Histogram

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Histogram[data, opts...]`**

Bins the numeric values in data and draws a frequency histogram. Bin count defaults to Sturges' rule: ceil(Log2\[n\]) + 1.

**`Histogram[data, k, opts...]`**

k equal-width bins.

**`Histogram[data, {step}, opts...]`**

Bins of width step.

**`Histogram[data, {min, max, step}, opts...]`**

Explicit range and width. Options: ChartStyle   color/style list cycling through bins BarSpacing   gap fraction of bin width (default 0.2) Standard Graphics options pass through.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Histogram[Table[RandomReal[], {200}]]
Out[1]= -Graphics-

In[2]:= data = Table[RandomReal[], {200}]
Out[2]= {0.178661, 0.430161, 0.812782, 0.00885941, 0.117009, 0.374991, 0.870921, 0.0176291, 0.867194, 0.0433318, 0.458816, 0.372966, 0.726014, 0.685522, 0.515293, 0.74346, 0.978886, 0.493972, 0.482616, 0.697654, 0.441669, 0.625915, 0.864002, 0.534648, 0.706739, 0.571163, 0.716441, 0.694264, 0.586844, 0.579504, 0.723413, 0.385733, 0.842941, 0.0999773, 0.253992, 0.972484, 0.208103, 0.954419, 0.568394, 0.364205, 0.649658, 0.556531, 0.581293, 0.272865, 0.880713, 0.460802, 0.579412, 0.694165, 0.509727, 0.986593, 0.889292, 0.72912, 0.18196, 0.278404, 0.467919, 0.0573136, 0.523487, 0.662998, 0.446285, 0.293101, 0.457043, 0.880796, 0.0443126, 0.28829, 0.484451, 0.0935614, 0.627013, 0.793831, 0.197005, 0.203763, 0.0579358, 0.450136, 0.626751, 0.377681, 0.304491, 0.0557184, 0.341145, 0.0620615, 0.127988, 0.0808755, 0.216694, 0.770363, 0.329797, 0.110209, 0.997586, 0.782126, 0.644404, 0.850379, 0.452234, 0.145283, 0.601646, 0.376486, 0.474005, 0.407595, 0.857994, 0.240346, 0.679083, 0.894269, 0.694642, 0.253708, 0.265091, 0.739935, 0.575176, 0.0416655, 0.833131, 0.681541, 0.341771, 0.866631, 0.302888, 0.0950284, 0.34083, 0.377937, 0.556943, 0.68804, 0.737725, 0.818262, 0.658377, 0.621306, 0.535673, 0.691901, 0.355946, 0.574331, 0.684926, 0.34357, 0.124224, 0.559872, 0.0227967, 0.32526, 0.664334, 0.354363, 0.245496, 0.666638, 0.551862, 0.228036, 0.304041, 0.701588, 0.0178955, 0.82289, 0.54641, 0.592203, 0.661386, 0.109424, 0.911152, 0.979182, 0.622362, 0.83585, 0.220201, 0.278937, 0.400055, 0.876589, 0.00859966, 0.798611, 0.546399, 0.734084, 0.687611, 0.965393, 0.599333, 0.289716, 0.257491, 0.319344, 0.593106, 0.313117, 0.788006, 0.513178, 0.824461, 0.463544, 0.343397, 0.782552, 0.644138, 0.644015, 0.762101, 0.833061, 0.988912, 0.236647, 0.306903, 0.316154, 0.699251, 0.381895, 0.608896, 0.282224, 0.75569, 0.327074, 0.20534, 0.570847, 0.443693, 0.590506, 0.532888, 0.996561, 0.877443, 0.459852, 0.32801, 0.411989, 0.810046, 0.370046, 0.0948037, 0.182889, 0.754547, 0.446398, 0.204962, 0.568322}

In[3]:= Histogram[data, 20]
Out[3]= -Graphics-

In[4]:= Histogram[data, {0.1}]
Out[4]= -Graphics-

In[5]:= Histogram[data, {0, 1, 0.05}]
Out[5]= -Graphics-
```

### Applications (5)

```mathematica
In[6]:= Histogram[{1, 2, 2, 3, 3, 3, 4, 4, 4, 4}]
Out[6]= -Graphics-

In[7]:= Histogram[Range[100], 10]
Out[7]= -Graphics-

In[8]:= Histogram[Range[100], {25}]
Out[8]= -Graphics-

In[9]:= Histogram[Range[100], {0, 100, 25}]
Out[9]= -Graphics-

In[10]:= Length[Cases[Histogram[Range[100], 10], _Rectangle, Infinity]]
Out[10]= 10
```

## Algorithm

barchart.c — BarChart[data, opts...] and Histogram[data, opts...]

BarChart renders a vertical bar chart from explicit heights. Histogram bins numeric data and renders a frequency histogram. Both return Graphics[...] objects auto-displayed by the REPL.

BarChart[{v1,...,vn}, opts...]

```text
  n bars at x = 1..n with heights v1..vn, coloured via ChartStyle or
  the default palette.
```

BarChart[{{v1,...}, {w1,...}, ...}, opts...]

```text
  Multiple grouped datasets, each dataset in a distinct palette colour.
```

Histogram[data, opts...]

```text
Histogram[data, k, opts...]          k equal-width bins
Histogram[data, {step}, opts...]     bins of given width
```

Histogram[data, {min,max,step}, opts...] explicit range + width

Options (both):

```text
  ChartStyle   — color/style list cycling through bars
  ChartLabels  — label list for x-axis ticks
  BarSpacing   — gap fraction of bar width (default 0.2)
  PlotLabel, standard Graphics options pass through 
```

## Implementation notes

**Algorithm.** `builtin_histogram` lives alongside `BarChart` in `barchart.c` and
shares its `split_chart_options` and `emit_bar` machinery. It is a plain
`Protected` builtin (**not** `HoldAll`): the data is evaluated and coerced to
reals. `histogram_parse_bins` reads the bin specification —
`Histogram[data]` uses **Sturges' rule** `n_bins = ceil(log2 n) + 1` clamped to
`[2, 50]`; `Histogram[data, k]` gives `k` bins; `Histogram[data, {step}]` gives
`ceil((dmax - dmin)/step)` bins over the data range; `Histogram[data, {min, max,
step}]` fixes the range and gives `ceil((max - min)/step)` bins. Counting: with
`bin_w = (bin_max - bin_min)/n_bins`, each value lands in bin
`(int)((v - bin_min)/bin_w)` clamped to `[0, n_bins-1]`, accumulated into an
`int` count array. Empty bins are skipped (no zero-height rectangle). Each
populated bin is drawn by the shared `emit_bar` as a 5-primitive bar over
`[x0, x0 + bin_w] x [0, count]`, coloured by `bar_color` indexed by the bin.
The result is an inert `Graphics[prims, opts...]` with `PlotRange` x = the bin
range and y = `[0, max_count*1.06 + 0.5]`.

**Data structures.** `ChartOpts`; a `calloc`'d `int counts[n_bins]`; `Expr**
prims` sized by the populated-bin count.

**Complexity / limits.** `O(n_vals + n_bins)`. Uniform-width bins only; Sturges'
count is hard-capped at 50; frequency counts only (no density/PDF
normalisation). Defaults `Axes -> True`, `AspectRatio -> 0.618` (`BarSpacing` is
parsed but unused in the histogram layout).

**Attributes:** `Protected`.

## References

**See also:** [ChartStyle](../../other-advanced/ChartStyle/), [BarSpacing](../../other-advanced/BarSpacing/)

- H. A. Sturges, *The Choice of a Class Interval*, J. Amer. Statist. Assoc. **21** (1926) 65-66 — the default bin count.
- Source: [`src/graphics/barchart.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/barchart.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)

## Notes & additional examples

### Notes

`Histogram` lives alongside `BarChart` and shares its bar-drawing machinery; it is
`Protected` but **not** `HoldAll`. The bin specification is flexible:
`Histogram[data]` uses **Sturges' rule** (`ceil(log2 n) + 1` bins, clamped to
`[2, 50]`); `Histogram[data, k]` gives `k` bins; `Histogram[data, {step}]` fixes the
bin width; `Histogram[data, {min, max, step}]` fixes the range and width.

Each value is dropped into its bin and the bin counts become bar heights; empty
bins are skipped (no zero-height rectangle). Only equal-width frequency bins are
supported — there is no density/PDF normalisation. Defaults are `Axes -> True`,
`AspectRatio -> 0.618`.
