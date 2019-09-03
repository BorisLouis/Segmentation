# Segmentation
Image processing for Polymer networks microscopy images (Segmentation, pore properties calculations, 3D rendering,...) in MATLAB.
The repository was used for the following Vandaele et al. Structural characterization of fibrous hydrogels using fluorescence microscopy. (INSERT LINK)

# User guide - WorkFlow

## Before using the code
The code is only taking .tif file, this was chosen because most software can convert to tif and was to avoid to have to plan for every file type that potential user would have. The code will make sure that the pixel size in XY and in Z is the same (by extrapolation if it is not the case), so if possible, it is better to acquire with the same Z step as the XY pixel size.

All the analysis detailed below can be run for testing on the testFile.tif that is in the repository.

## 1. Segmentation
To perform segmentation on the data, one need to use **mainImSegmentation.m**
Within the code there is a section **User Input** within which the parameter for the segmentation (sensitivity, cleanup,...) can be tuned although I would recommend to start with default. Each parameter is briefly explained in comment in the code. Once you have set the parameter, simply press run, it will open a dialog where you can indicate a **folder** in which your data is stored. In this folder, it will search and segment all the .tif files contained and stored the resulting segmentation in a folder called **SegmentedStacks** inside the user-indicated folder. The analysis will be run with adaptive threshold and single-threshold. It is important that you check which one works best on your data so you can choose after which one needs to be used for the property determination. **!!! The segmented image are inverted compare to the image so what is bright in the image will appear dark and vice versa because the rest of the code analyze only the bright part which therefore needs to be the pores**. if needed, Inversion of the data can be perform with ~ before the binary variable containing the data or using imcomplement (in matlab).

To check how good/bad is the segmentation, you can use **mainCheckSegmentation.m**.
Pressing run will aslo open a dialog where the user can select the **folder where segmented images are stored as .tif, SegmentedStacks folder if our routine was used for segmentation**.

I would recommend to first run and check the segmentation on a single file (to avoid long processing time) and then, once you are satisfied with your parameters, put all the files acquired in similar condition in one folder to segment them all in a single "run button" press.

## 2. Pore properties determination
To calculate the properties of the pores, one need to use **mainPorePropsCalc.m**.
The main information that need to be provided here is the pixel size in XY and in Z and whether you want the analysis to be run in 2D, 3D or both. Once you press run, a dialog where the user can select **folder where segmented images are stored as .tif, SegmentedStacks folder if our routine was used for segmentation**. The results will be stored in a folder called PoreSize-Results.

## 3. Data Visualization
To obtain 3D model of the network one need to use **main3DModel.m**.
The main information that needs to be provided here is the pixel size in XY and in Z for scaled representation., color of the model and color of the pores. Once again, pressing run will open a dialog where the user need to select **A binary z-stack folder **. In the indicated folder it will search for a PoreSize-Results folder and check if there is a .mat file in there. If yes, it will also output model of the pores connectivity and the polymer model together with the pores, if not it will only output a model of the network.

To obtain distribution plot of different properties one need to use **plotPoreprops.m**
Once you press run, a dialog where the user can select **folder where PoreSize results are stored as .mat, PoreSize-Results folder.**
The main information that needs to be provided is some details about the plotting (e.g. stepSize, number of bins, log scale or not), if you want to plot 2D or 3D determined properties and which propertie you want to plot (volume, internal radius, throat, diameter,...).
The code will save a figure for each .mat file present in the selected folder. Therefore, if you want to compare/look at different conditions one can put the results output from different conditions in a new folder (give meaningful name to the .mat file as they will be reuse in the figure names) and then run plotPoreProps which will output and save (in a folder called figure) a figure for each file (thus each conditions). 

**Watchout that the results from different z-stack stored in the same folder (segmented together) will be made into one single distribution. This is okay if they were acquired in the same condition (e.g. in our paper, same polymer length, same concentration, same imaging conditions) but may lead to unwanted results if one mix data that should not be analyzed together.**



