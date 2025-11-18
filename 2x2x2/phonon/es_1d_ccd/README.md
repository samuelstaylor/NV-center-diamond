This directory: performs ES calculation (DeltaSCF or TDDFT) on each interpolated geometry between GS geo 3A_2 and ES geo 3E

Steps: set up the input files:

gs_coord.in: the coordinated of relaxed gs

es coords-- two options: 
`es_cdft_010_coord.in`
`es_tddft_010_coord.in`

I will name choose the second one since that is what used to calculate the es coords

`prefix.in` will be appended with the run settings to the coordinates for each image. 
Since we want to do a ES calculation, we can choose either TDDFT or CDFT
- lets have it with the same parameters as our tddft calculation