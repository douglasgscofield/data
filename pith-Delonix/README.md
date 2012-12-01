Pith cell sizes in twigs of *Delonix regia*
===========================================

Complete data from [Scofield 2006](http://dx.doi.org/10.3732/ajb.93.12.1740), along with the R script used for data analysis and figure generation.

Scofield, DG.  2006.  Medial pith cells per meter in twigs as a proxy for mitotic growth rate (&Phi;/m) in the apical meristem.  *American Journal of Botany* **93**(12): 1740-1747.  <http://dx.doi.org/doi:10.3732/ajb.93.12.1740>

When using these data, please cite the above paper as well as the DataDryad entry for the processed data:

Scofield, DG.  2006.  Data from: Medial pith cells per meter in twigs as a proxy for mitotic growth rate (&Phi;/m) in the apical meristem.  Dryad Digital Repository.  <http:http://dx.doi.org/10.5061/dryad.b1t2b>


I have waived copyright and related or neighboring rights to these data in all territories under the Creative Commons 0 license   <a rel="license" href="http://creativecommons.org/publicdomain/zero/1.0/"><img src="http://i.creativecommons.org/p/zero/1.0/80x15.png" style="border-style: none;" alt="CC0" /></a>.

You can everything here bundled into a single zip file: [pith-Delonix.zip](https://github.com/downloads/douglasgscofield/data/pith-Delonix.zip)

* * *

This work was part of my doctoral research with Stewart T. Schultz and Carol C. Horvitz at the University of Miami and was funded in part by U.S. National Science Foundation Doctoral Dissertation Improvement Grant DEB-0309253 and grants from Sigma Xi and the University of Miami Department of Biology.

* * *

### Contents

Columns in `Scofield-AJB-2006-Delonix-regia-pith-measurements.txt` (the
processed data file deposited in DataDryad) are:

* tree: the individual tree from which the twig was collected
* branch: the individual branch on the tree from which the twig was collected, up to 4 
  within a tree
* segment: the segment within the twig from which sections were taken and cell size 
  measurements made, up to 3 within a branch
* file: location of the measured cell file within the pith (marginal or medial)
* len: calculated length of an individual pith cell in micrometers, see below for more 
  details and caveats re: variance
* dbh: diameter of the tree at 1.4 m height, NA if not measured

Thus every cell length measurement is nested within segment, which is nested
within branch, which is nested within tree. Hence the random-effects analysis
presented in the paper.

Files included here:

* `README.DataDryad.txt`: the README submitted to DataDryad 
* `Scofield-AJB-2006-Delonix-regia-pith-measurements.txt`: the processed data as deposited 
  in DataDryad
* `data_cellsizes.txt`: raw measurements of cell sizes; each row is a collection of 
  measurements for a single slide, with values being the length of cell files in number of 
  ocular micrometer gradations
* `data_trees.txt`: raw measurements of tree circumference (m) and calculated dbh (m)
* `Pith.R`: [R](http://www.r-project.org) source file for processing raw data, running 
  analyses and producing figures for the above paper.
* `Delonixregia_pith_medialfile.png`: Figure 1A below
* `Delonixregia_twig_crosssection.png`: Figure 1B below
* `README.md`: this README

To regenerate the processed data, plus a couple of transformed columns, and to
create Figure 2 from the paper, put all files in the same directory, fire up R,
and enter:

````R
source("Pith.R")
dat = pith.readdata(include.dbh=TRUE)
pith.figure2(dat)
````

See below for more details on data gathering and analysis, and see the paper
for full details and relevant references.  The material below is a cut-down and
slightly edited version of the Materials and Methods section and Figure 1 of
the above paper.

If you have any further questions or something here isn't clear, feel free to
[email me](mailto:douglasgscofield@gmail.com).


### Sample collection and preparation

[*Delonix regia*](http://en.wikipedia.org/wiki/Delonix_regia)
(Caesalpinioideae) is a fast-growing, familiar tree planted as an ornamental
throughout the world’s tropical and subtropical regions, originally native to
the lowland seasonal tropics of Madagascar but now very rare there.

Trees used in this study were cultivated specimens growing on the University of
Miami campus, Coral Gables, Florida, USA (25º 43′ N, 80º 16′ W). In January
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


#### Figure 1A: *Delonix regia*, twig cross-section

![*Delonix regia* twig cross-section showing pith and xylem](https://github.com/douglasgscofield/data/blob/master/pith-Delonix/Delonixregia_twig_crosssection.png?raw=true)

#### Figure 1B: *Delonix regia*, file of 10 medial pith cells (length 386 μm as measured)

![*Delonix regia* file of 10 medial pith cells as measured in the attached data](https://github.com/douglasgscofield/data/blob/master/pith-Delonix/Delonixregia_pith_medialfile.png?raw=true)



### Pith cell size measurements

Cells were organized into clear longitudinal files, with cells within medial
files clearly larger than cells within marginal files (see Figure 1A above).
The length of medial cells was measured along the longitudinal axis, using a
Reichert microscope with ocular micrometer.  The ocular micrometer was
calibrated with a 2 mm × 0.01 mm scale stage micrometer prior to the first
measurement.  Lengths were estimated to the nearest ½ unit, with the length of
a cell extending between the middle of the cell walls shared with neighboring
cells in the same file along the longitudinal axis. To reduce measurement
error, the aggregate length of ten contiguous cells within the same file was
measured (see Figure 1B above).  Mean individual cell length was thus
calculated as 0.1× the length measured. The variance of individual cell lengths
reported here is thus reduced, and should be multiplied by 10 (√10 ≈ 3.2 for
standard deviations and standard errors) to derive appropriate variances for
lenghts of single cells. Whenever possible for each segment, 10 files of 10
medial pith cells were chosen haphazardly for measurement from among the
available files of medial cells, along with 4 files of 10 marginal pith cells.

There were various sources of data loss during sample preparation and
measurement, and the dataset was imbalanced.  A balanced subset was subjected
to most analyses in the paper.  **The attached data include all data, including
data for test slides (prefixed with C) used while refining the methods.  Note
in particular that only the lengths of medial files were reported in the above
paper; this dataset also includes lengths of marginal files.**

### Licensing

When using these data, please cite the above paper and DataDryad data repository.

<p xmlns:dct="http://purl.org/dc/terms/">
  <a rel="license"
     href="http://creativecommons.org/publicdomain/zero/1.0/">
    <img src="http://i.creativecommons.org/p/zero/1.0/80x15.png" style="border-style: none;" alt="CC0" />
  </a>
  <br />
  To the extent possible under law,
  <a rel="dct:publisher"
     href="https://github.com/douglasgscofield/data/tree/master/pith-Delonix">
    <span property="dct:title">Douglas G. Scofield</span></a>
  has waived all copyright and related or neighboring rights to
  <span property="dct:title">Pith cell sizes in twigs of Delonix regia</span>.
</p>

