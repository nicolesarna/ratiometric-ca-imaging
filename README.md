# Ratiometric Calcium Imaging - Single Timepoint, 2D

MATLAB image processing pipeline for ratiometric intracellular calcium analysis of single-timepoint
confocal images acquired with the Fluo-4 / Fura-Red dye pair.

The script segments individual cells in each image, assigns a `CellID` per object,
performs background-corrected fluorescence quantification on the **raw** channel data,
and exports a table with one row per cell.

## Author

**Nicole Suzanne Sarna** - Department of Bioengineering, Rice University
(Laboratory of Dr. Michael R. King)

Written April 30, 2026; last updated May 20, 2026.

## Citation

If you use this software, or data produced with it, please cite:

> Sarna, N. S. (2026). *RatiometricCaImagingSingleTimepoint: MATLAB analysis program
> for ratiometric Fluo-4/Fura-Red calcium imaging* (Version 1.0.1) [Computer software].
> Zenodo. https://doi.org/10.5281/zenodo.24380711

Suggested methods-section language:

> Ratiometric calcium measurements were obtained using a custom MATLAB image analysis program written by
> N. S. Sarna (King Laboratory, Department of Bioengineering, Rice University), which performs
> automated cell segmentation and background-corrected quantification of Fluo-4 and Fura-Red
> fluorescence intensity (doi:10.5281/zenodo.24380711).

## Provenance

This script is adapted from [TCellAnalysis](https://github.com/nicolesarna/TCellAnalysis)
by the same author. The underlying segmentation and fluorescence quantification approach is described in:

> Sarna, N. S., Desai, S. H., Kaufman, B. G., Curry, N. M., Hanna, A. M., & King, M. R. (2024).
> Enhanced and sustained T cell activation in response to fluid shear stress.
> *iScience*, 27(6), 109999. https://doi.org/10.1016/j.isci.2024.109999

The present version extends that pipeline to two-channel ratiometric calcium measurement
with the Fluo-4 / Fura-Red dye pair, and retunes the segmentation parameters for
single-timepoint confocal acquisition.

Work reporting results produced with this script should cite both the software record
above and the iScience article.

## What the script does

1. Prompts for one or more single-timepoint `.czi` files (two channels: ch1 = Fluo-4, ch2 = Fura-Red).
2. Reads pixel calibration from the OME metadata, falling back to a hardcoded value if absent.
3. Builds a binary cell mask by median filtering, contrast adjustment, thresholding, and
   morphological opening of the combined channels, then separates touching cells by
   distance-transform watershed.
4. Clears border objects and removes objects below 1000 px^2.
5. Quantifies the **raw, unmodified** channel data within each object: integrated density,
   corrected total cell fluorescence (CTCF), and background-corrected mean (BCmean).
6. Computes the per-cell Fluo-4 / Fura-Red ratio from the background-corrected means.
7. Optionally writes QC segmentation overlays and a three-panel channel montage with a scale bar.
8. Optionally exports the per-cell table to Excel.

Enhancement and filtering are used **only** to generate segmentation masks. All reported
fluorescence values are measured on the raw channels.

## Outputs

All outputs are written to the **same folder as the selected `.czi` files**. Nothing is
written anywhere else, and no files are overwritten except as noted below.

| File | Written when | Contents |
|---|---|---|
| `<imagename>_segmentation.png` | "Save QC segmentation images?" = Yes | One per input image. QC overlay. |
| `<imagename>_montage.png` | "Save QC segmentation images?" = Yes | One per input image. Display figure. |
| `SingleTimepointCalciumAnalysis.xls` | "Export data to Excel?" = Yes | One workbook for the whole batch, one row per cell. |

### `<imagename>_segmentation.png` - segmentation QC

A composite of the two filtered channels, overlaid with magenta outlines around every
object that passed segmentation, and the `CellID` number printed in red at each object
centroid. The title gives the image name and the number of cells found.

This is the file to check before trusting any row of the table. It lets you confirm that
cells were split correctly, that debris and partial cells at the image border were
excluded, and that any particular row in the spreadsheet can be traced back to the cell it
came from — `CellID` in the table matches the number drawn on the image.

### `<imagename>_montage.png` - three-panel display figure

A publication-style figure with three panels side by side on a white background -
**Fura-Red** (red), **Fluo-4** (green), and **Composite** (both) - each labelled above and
carrying a 20 µm scale bar in the lower right corner.

This figure is for display only. The panels are built from the contrast-adjusted,
median-filtered channels (linear adjustments for visualization purposes). 
**Do not measure intensities from these images** - all quantitative values come from the 
raw channel data and appear in the spreadsheet.

Regarding the scale bar: its length is computed from the pixel size read
from the image metadata, so it is correct for each image individually; but the `20 µm`
text label is commented out in the source, so the bar is drawn unlabelled. Either state
the scale in the figure caption, or uncomment the `text(...)` call in the montage block.

### `SingleTimepointCalciumAnalysis.xls` - per-cell data

A single workbook pooling every image in the batch, one row per cell, sorted by image and
then by `CellID`. Columns are described below.

**This filename is fixed.** Running the script a second time on images in the same folder
overwrites the previous spreadsheet without warning. Rename or move the file between runs
if you are processing several conditions out of one directory. The PNG filenames are
derived from the image names, so those do not collide across images, but re-running on the
same images does regenerate them.

Figures are created invisibly and closed automatically, so nothing appears on screen during
a batch — the only sign of progress is the console, which prints a warning for any image
where the pixel size could not be read, or where no cells or no background region were
found (such images are skipped and contribute no rows).

## Quantities reported

| Column | Meaning |
|---|---|
| `FileName`, `File`, `CellID` | source image and per-image object index |
| `Area_px`, `Area_um2` | object area in pixels and calibrated units |
| `Diameter_px`, `Diameter_um` | equivalent circular diameter |
| `Circularity`, `Eccentricity`, `Solidity` | shape descriptors |
| `RedCTCF`, `GreenCTCF` | corrected total cell fluorescence, `IntDen - (area_px x mean_background)` |
| `RedBCmean`, `GreenBCmean` | `CTCF / area_px`, equivalently `mean - mean_background` |
| `RatioBCmean_Fluo4overFuraRed` | `GreenBCmean / RedBCmean` — the ratiometric calcium readout |

## Requirements

This code was written and run in **MATLAB R2025a**. 
It should run without modification on **R2020a or later**

| Requirement | Notes |
|---|---|
| MATLAB R2020a or later | Developed and tested on R2025a |
| Image Processing Toolbox | **Required.** Nearly every segmentation, morphology, and display call depends on it (`imbinarize`, `bwareafilt`, `watershed`, `bwdist`, `regionprops`, `imoverlay`, and others). |
| Bio-Formats MATLAB toolbox (`bfmatlab`) | **Required**, for `bfGetReader` and `bfGetPlane`, which read the `.czi` files and their OME metadata. Must be on the MATLAB path. |

### Installing Bio-Formats

Bio-Formats is developed by the Open Microscopy Environment and is not redistributed here.
Download `bfmatlab` from the [Bio-Formats downloads page](https://www.openmicroscopy.org/bio-formats/downloads/),
unzip it, and add the folder to your MATLAB path:

```matlab
addpath('/path/to/bfmatlab')
savepath
```

Use a Bio-Formats release compatible with your MATLAB version - each release supports a
specific range of MATLAB and Java versions, documented in the
[Bio-Formats documentation](https://bio-formats.readthedocs.io/). Mismatched versions
typically fail with a Java error on the first call to `bfGetReader`.

To confirm the installation:

```matlab
which bfGetReader   % should return a path inside bfmatlab
bfCheckJavaPath     % should return 1
```

## Usage

Add `bfmatlab` to the MATLAB path, then run:

```matlab
RatiometricCaImagingSingleTimepoint
```

Select one or more `.czi` files from a single directory (shift-click for multiple), then
answer the two prompts for QC images and Excel export. Outputs are written to the same
directory as the selected images.

### Parameters you may need to retune

The threshold and morphology values were set on a representative image and are not
universal. For a new imaging configuration, check the QC overlays and adjust:

- `imbinarize` thresholds for the two channels (currently 0.2 and 0.1)
- `imextendedmin` depth for watershed seeding (currently 7)
- the minimum object size in `bwareafilt` (currently 1000 px)
- the `strel` disk radii

Montage cosmetics are collected at the top of the figure-saving block and flagged with
`EDIT ... here` comments: `labelFontSz` (channel label size, 40), `sbarFontSz` (scale bar
label size, 30), `sbarLineW` (scale bar thickness, 10), and `pad` (white space around and
between panels, 40 px). The scale bar length is set by the `20` in
`scaleBarPx = round(20 / pxSize_um)`; change that value to draw a different bar length.

## Use and attribution

You are welcome to use and adapt this code for research. Please keep the authorship header
in the source file intact, and credit the author in any thesis, presentation, or publication
reporting results produced with it.
