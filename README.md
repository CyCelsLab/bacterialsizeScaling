# Bacterial size scaling
Code relating to the scaling of bacterial cell size

## Instructions for use (FM464 data)
1. Download all contents of folder Code
2. Open file _SegmentationTest.m_ in MATLAB and run it. 
3. A figure window will open with cell overlays detected using the pipeline described in the paper. Select contours that are not correct. 
4. Output will be an overlay of detections and length and width of selected cell contours in pixel units. 


## Instructions for use (Phase Contrast)
1. Download all contents of folder PhaseContrast
2. Add the additional folders to the matlab path (bresenham/interparc/interpclosed)
2. Open file _PhaseSementation_v3.m_ in MATLAB and run it. 
3. A figure window will open with cell overlays detected using the pipeline described in the paper.
4. Output will be an overlay of detections and length and width of selected cell contours in pixel units.
5. A final figure of a L & W plot and histogram of Lengths and widhts will show up.


### List of toolboxes required
1. Image Processing Toolbox
2. Curve Fitting Toolbox
3. Partial Differential Equation Toolbox

### REFERENCE:
In case this code is used in any settings please refer to our paper that corresponds to this
##### Tanvi Kale, Dhruv Khatri and Chaitanya A Athale (2023) Allometry of Escherichia coli surface area with volume: effect of size variability, filamentation and division dynamics. Phys. Biol. 20 046007. DOI:10.1088/1478-3975/acdcda
