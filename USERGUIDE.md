# HyperViewer User Guide

User Guide for [HyperViewer software](www.github.com/springlabnu/hyperviewer). To get started, checkout the [README.md](www.github.com/springlabnu/hyperviewer/README.md) file.

## Contents

* [Introduction](#introduction)
* [Workflow](#workflow)
* [GUI Layout](#gui-layout)
* [Menu Bar](#menu-bar)

## Introduction

The HyperViewer software was written to enhance visualization and analysis of hyperspectral images and videos through spectral unmixing.

## Workflow


### Basis Spectra Library

The Basis Spectra Library (BSL) contains all the information needed for spectral unmixing, along with identifying information used by HyperViewer to display the spectrum in the GUI. The library, contained in ```BSLpub.mat```, is a ```struct``` with the following fields:

* **name**
  * ```char``` Name of the spectral source, like 'AF633' or 'DAPI'
* **color**
  * ```char``` Default color when displayed as a spectrum or basis map. Must be interpretable by ```rgb.m```.
* **target**
  * ```char``` *Optional*. Name of the object associated with the spectral source, like 'Apple' or 'EGFR'.
* **tspec**
  * n-by-3 ```double``` *Optional*. Theoretical spectral characteristics, like those downloaded from the web. First column contains the wavelength in nm, second column contains the theoretical absorption/excitation spectrum, third column contains the theoretical emission spectrum. If one is unknown, populate that column with zeros.
* **mspec**
  * n-by-2 ```double``` Spectral information as measured by the relevant hyperspectral imager. First column is the wavelength in nm, second column is the spectral emission. This data is used for unmixing and should be at the same spectral resolution as the image (that is, the number of wavelength channels should be equal). If not, the spectrum in ```mspec```is interpolated to match the resolution of the hyperspectral image.
* **mpespec**
  * n-by-3 ```double``` *Optional*. Spectral characteristics for multi-photon excitation. Reserved for future use.
* **compIdx**
  * ```logical``` If true, this spectrum will be included in the unmixed composite image. If false, this spectrum will not be included in the composite by default. Some spectra, like autofluorescence spectra, are best left excluded from the final image sets.   
* **abs_max**
  * ```double``` *Optional*. Wavelength of peak absorption
* **ems_max**
  * ```double``` *Optional*. Wavelength of peak emission.
* **filetype**
  * ```char``` Hyperspectral image filetype associated with the spectrum.

The BSL is loaded automatically on startup. The user can move spectra from the BSL to the A matrix using the left/right arrows between the two lists.

There are a few ways to edit the BSL:

##### Variables editor

In MATLAB, navigate to the HyperViewer working directory and load the BSL into the workspace.
```
>> load(['lib' filesep 'BSLpub.mat'])
```
Double click the ```BSL``` variable in the workspace to open it in the Variables editor. Manually add or edit entries by typing information into the matrix. To edit strings, surround the text with single quotes. To add matrices, type the variable name into the appropriate box. When finished, save the new BSL.
```
>> save(['lib' filesep 'BSLpub.mat'], 'BSL')
```

##### Within the GUI

* [Edit Colors](#edit-colors) allows the user to edit some of the properties of the BSL within the GUI. Changes are retained through the session but are not saved permanently.
* [Region of Interest -> Add ROI to BSL](#region-of-interest) adds the spectra defined in the current ROI into the BSL. These changes are saved to the BSL file permanently.

### Hyperspectral Images

An 'experiment' is an image or group of images meant to be unmixed and analyzed with the same BSL. Images should all be saved in the same directory. For filetypes that are saved as tiff series, each image set should be saved in its own subdirectory. When an experiment is loaded, HyperViewer looks for all images within the directory and stages them for analysis.

To avoid extensive loading times, only the first image is loaded and unmixed. Subsequent images are loaded and unmixed as they are selected in the image directory. The 'Unmix All' button will load and unmix all remaining images.




## GUI Layout




## Menu Bar

### File

#### New -> *Ctrl N*
Load a new Experiment.

#### Open -> *Ctrl O*

NOT Functional. Reserved for future versions.

#### Reload

Load the current Experiment fresh using the selected directory. This resets all analyses.

#### Save Analysis-> *Ctrl S*

Save the unmixing analysis of the currently selected image.

#### Save All Analyses

Save the unmixing analyses of all images in the Experiment.

#### Save New Image

Save the current image as new hyperspectral image.

#### Export As

Save the current analysis in a different configuration.

* **Movie**
  * Combines images and analyses and the Experiment into a ```.avi``` video file. Default frame rate is 15 fps, but can be changed in the code.
* **New Basis Spectra**
  * Exports the average spectrum within the most recently created ROI to a text file.
* **Tiles**
  * Creates a sequence of images stitched together with a buffer between images. A rectangle appears for the user to select the region to use in the tile. The composite is placed first, followed by the basis maps in the order of the A matrix library. This feature is useful for creating unmixing grids for visual comparison of spatial patterns.
  * Options to customize the output tile:
    * ```direction```: direction of the stitching. ```'horizontal'``` or ```'vertical'```
    * ```sizex```, ```sizey```: Length and width (in pixels) of the default rectangle. A non-negative integer.
    * ```pxlbuffer```: Width (in pixels) of the buffer between stitched images. A non-negative integer.
    * ```buffercolor```: Color of the pixel buffer. String of any color accepted by ```rgb.m```.



### Tools

#### Edit Colors -> *Ctrl B*

Opens the BSL in a table format to edit BSL entries. The properties that can be edited are: ```name```, ```color```, ```target```, and ```compIdx```. After the desired changes are made, press 'Done', and the new BSL will update in the GUI.

#### SpectraViewer GUI

Not Functional; reserved for future versions. Use [Edit Colors](#edit-colors) to edit the BSL.

#### Region of Interest

The default ROI loaded for every image is the entire image. ROIs or single points can be subsequently added or replaced using the options below creating a stack of ROIs on the image. *Adding* ROIs place additional regions onto the image, while *New* ROIs replace most recent ROI in the stack. The analysis axes always displays information from newest ROI in the stack. ROIs are removed in the reverse order from which the were created, the oldest ROI is never removed.

* **New ROI** -> *Ctrl R*
  * Replace the most recent ROI or point in the stack with a new ROI.
* **Add ROI**
  * Add a new ROI into the stack.
* **New Point** -> *Ctrl A*
  * Replace the most recent ROI or point in the stack with a new point.
* **Add Point**
  * Add a new point into the stack.
* **Add ROI to BSL**
  * Place the newest ROI in the stack into the BSL. Information about the spectrum is given through a popup box.
* **Clear Last ROI** -> *Ctrl Z*
  * Remove the newest ROI from the ROI stack.

#### Average Experiment

Average images together to reduce noise. This only applies to multiple images of the same location, like a long exposure. A dialogue box appears for the user to select the images to average, the default selection is the entire experiment but can be reduced to a subset of *sequential* images from the experiment.

Averaging the raw images (before unmixing) or the basis maps (after unmixing) may produce different results and both may be desirable for different reasons. The user can try both options:

* **Pre-Unmixing**
  * Average only the raw hyperspectral images. Then unmix the result.
* **Post-Unmixing**
  * Average the subset of unmixed abundance maps into a new basis set.


### Options

Options are treated as logical (true/false) selections.

#### Include All Spectra In Composite

Default: False

Some spectra, like autofluorescence, are excluded from the composite by default. If true, the composite is recalculated using *all* spectra from the A matrix list.  

#### Include Saturated Pixels

Default: False

For some image types, saturated pixels are tracked and excluded from the analysis. If true, this will not mask the image and include saturated pixels.

#### Gaussian Blur

Default: True

Adds a gaussian blur of 4 pixels to the image. The width of the gaussian can be changed in the GUI parameters.

#### Global ROI

Default: False

By default each image in the experiment keeps track of it's own ROIs. When Global ROI is on, the ROI stack of the current image is applied to the entire experiment.

#### Apply Fiber Image Filter

Default: False

For SMIRC data only. In some data sets, a binary mask of the fiber core can be applied.

#### Sampling Correction

For SMIRC data only. In some data sets, images are over sampled. This option applies a convolution filter to condense the image down to the correct aspect ratio.

### Demo

#### Apples

Analyze a scene of real and plastic apples taken using the [Specim IQ](https://www.specim.fi/iq/) camera. Expected unmixing time is ~20 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Supplementary Figure 6.

#### Wide-field Fluorescence Imaging of Ovarian Cancer

A mouse model of advanced-stage disseminated ovarian cancer was treated with an anti-EGFR antibody with a fluorescent label and imaged using the CRi Maestro small animal imager. Expected unmixing time is ~20 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Supplementary Figure 5.

#### SMIRC - *In vivo*

A mouse model of advanced-stage disseminated ovarian cancer was imaged via 5-color SMIRC. Expected unmixing time is ~15 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Figure 1d.

#### SMIRC - Dyes in Solution

A set of Alexa Fluor Dyes were dissolved in solution and imaged via SMIRC. Expected unmixing time is ~5-10 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Figure 1a using the image tiling feature.
