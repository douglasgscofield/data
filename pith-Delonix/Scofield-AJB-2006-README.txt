Pith cell sizes in twigs of *Delonix regia*
===========================================

Data from:

Scofield, DG.  2006.  Medial pith cells per meter in twigs as a proxy for
mitotic growth rate (Φ/m) in the apical meristem.  American Journal of Botany
93(12): 1740-1747.  http://dx.doi.org/doi:10.3732/ajb.93.12.1740

Columns:

tree: the individual tree from which the twig was collected

branch: the individual branch on the tree from which the twig was collected,
    up to 4 within a tree

segment: the segment within the twig from which sections were taken and
    cell size measurements made, up to 3 within a branch

file: location of the measured cell file within the pith (marginal or medial)

len: calculated length of an individual pith cell in micrometers (see below
    for more details and caveats re: variance)

dbh: diameter in meters of the tree at 1.4 m height, NA if not measured

Thus each cell size measurement is nested within segment, which is nested
within branch, which is nested within tree.  These data were subject to a
mixed-effects linear model analysis in the above paper.

See below for more details on data gathering and analysis, and see the paper
for full details and relevant references.  Raw data files and the R script
used for analysis and figure production are available at my github repository 
for the dataset:

    http://github.com/douglasgscofield/data/tree/master/pith-Delonix

The material below is a cut-down and slightly edited version of the Materials
and Methods section of the above paper.

* * *

This work was part of my doctoral research with Stewart T. Schultz and Carol
C. Horvitz at the University of Miami and was funded in part by U.S. National
Science Foundation Doctoral Dissertation Improvement Grant DEB-0309253 and
grants from Sigma Xi and the University of Miami Department of Biology.

* * *


### Sample collection and preparation

*Delonix regia* (Caesalpinioideae) is a fast-growing, familiar tree planted as
an ornamental throughout the world’s tropical and subtropical regions,
originally native to the lowland seasonal tropics of Madagascar but now very
rare there.

Trees used in this study were cultivated specimens growing on the University
of Miami campus, Coral Gables, Florida, USA (25º 43′ N, 80º 16′ W). In January
2003, I collected 4 apical twigs approximately 1 m in length and no more than
~2 cm in diameter at their base from each of 25 separate *D. regia*. Each twig
was collected from a different primary branch attached to the main trunk, and
twig-bearing branches were chosen from several locations varying in exposure
throughout the canopy.  At the time of twig collection, tree size was measured
as trunk diameter at 1.4 m height (dbh).  All subject trees had a single
unbranched trunk to at least 1.4 m.

Three 1 cm segments were cut from internodal locations selected at random
between 30 cm and 60 cm from the growing tip; segments within this range had
diameters of ~1 cm; this region was chosen in an attempt to avoid acropetal
regions of the twig that might still be undergoing internode elongation. Each
segment was stripped of its bark and cut in half longitudinally, and one half
of each segment was selected at random. The three segments of each twig were
then impregnated with Paraplast+ paraffin and four radial or nearly-radial
longitudinal sections of thickness 7 to 8 μm were taken from each segment on a
sliding microtome.   Sections were deparaffined and stained in aqueous 0.05%
toluidine blue O and 0.1% safranin O.


### Pith cell size measurements

Cells were organized into clear longitudinal files, with cells within medial
files clearly larger than cells within marginal files.  The length of medial
cells was measured along the longitudinal axis, using a Reichert microscope
with ocular micrometer.  The ocular micrometer was calibrated with a 2 mm ×
0.01 mm scale stage micrometer prior to the first measurement.  Lengths were
estimated to the nearest 1/2 unit, with the length of a cell extending between
the middle of the cell walls shared with neighboring cells in the same file
along the longitudinal axis. To reduce measurement error, the aggregate length
of ten contiguous cells within the same file was measured.  Mean individual
cell length was thus calculated as 0.1× the length measured. The variance of
individual cell lengths reported here is thus reduced, and should be
multiplied by 10 (√10 ≈ 3.2 for standard deviations and standard errors) to
derive appropriate variances for lenghts of single cells.  Whenever possible
for each segment, 10 files of 10 medial pith cells were chosen haphazardly for
measurement from among the available files of medial cells, along with 4 files
of 10 marginal pith cells.

There were various sources of data loss during sample preparation and
measurement, and the dataset was imbalanced.  A balanced subset was subjected
to most analyses in the paper.  The attached data include all data,
including data for test slides (prefixed with C) used while refining the
methods.  Note in particular that only the lengths of medial files were
reported in the above paper; this dataset also includes lengths of marginal
files.

