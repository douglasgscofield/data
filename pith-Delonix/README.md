Pith cell sizes in twigs of *Delonix regia*
===========================================

Complete data from [Scofield 2006](http://dx.doi.org/10.3732/ajb.93.12.1740), along with the R script used for data analysis and figure generation.

Scofield, DG.  2006.  Medial pith cells per meter in twigs as a proxy for mitotic growth rate (&Phi;/m) in the apical meristem.  *American Journal of Botany* **93**(12): 1740-1747.  <http://dx.doi.org/doi:10.3732/ajb.93.12.1740>

All data here are released from copyright under the Creative Commons 0 license.  When using these data, please cite the above paper as well as the DataDryad entry for the processed data:

### Contents

Columns in `Scofield-AJB-2006-Delonix-regia-pith-measurements.txt` (the processed data file deposited in DataDryad) are:

* tree: the individual tree from which the twig was collected
* branch: the individual branch on the tree from which the twig was collected, up to 4 within a tree
* segment: the segment within the twig from which sections were taken and cell size measurements made, up to 3 within a branch
* file: location of the measured cell file within the pith (marginal or medial)
* len: calculated length of an individual pith cell in micrometers, see below for more details and caveats re: variance
* dbh: diameter of the tree at 1.4 m height, NA if not measured

Thus every cell length measurement is nested within segment, which is nested within branch, which is nested within tree.


Files included here:

* `README.DataDryad.txt`: the README submitted to DataDryad 
* `Scofield-AJB-2006-Delonix-regia-pith-measurements.txt`: the processed data as deposited in [DataDryad]()
* `data_cellsizes.txt`: raw measurements of cell sizes; each row is a collection of measurements for a single slide, with values being the length of cell files in number of micrometer gradations
* `data_trees.txt`: raw measurements of tree circumference (m) and calculated dbh (m)
* `Pith.R`: [R](http://www.r-project.org) source file for processing raw data, running analyses and producing figures for the above paper.  For example, to regenerate the processed data along with Figure 2 from the paper, put all files in the same directory, start R, and enter:
````R
source("Pith.R")
dat = pith.readdata(include.dbh=TRUE)
pith.figure2(dat)
````
* `Delonixregia_pith_medialfile.png`: Figure 1A below
* `Delonixregia_twig_crosssection.png`: Figure 1B below
* `README.md`: this README

See below for more details on data gathering and analysis, and see the paper for full details and relevant references.  The material below is a cut-down and slightly edited version of the Materials and Methods section and Figure 1 of the above paper.

If you have any further questions or something here isn't clear, feel free to [email me](mailto:douglasgscofield@gmail.com).


### Sample collection and preparation


Trees used in this study were cultivated specimens growing on the University of Miami campus, Coral Gables, Florida, USA (25º 43′ N, 80º 16′ W). In January 2003, I collected 4 apical twigs approximately 1 m in length and no more than ~2 cm in diameter at their base from each of 25 separate *D. regia*. Each twig was collected from a different primary branch attached to the main trunk, and twig-bearing branches were chosen from several locations varying in exposure throughout the canopy.  At the time of twig collection, tree size was measured as trunk diameter at 1.4 m height (dbh).  All subject trees had a single unbranched trunk to at least 1.4 m.



#### Figure 1A: *Delonix regia*, twig cross-section

![*Delonix regia* twig cross-section showing pith and xylem](Delonixregia_twig_crosssection.png)

#### Figure 1B: *Delonix regia*, file of 10 medial pith cells (length 386 μm as measured)

![*Delonix regia* file of 10 medial pith cells as measured in the attached data](Delonixregia_pith_medialfile.png)



