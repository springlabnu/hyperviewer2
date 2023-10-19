# HyperViewer

Visualization and analysis tool for hyperspectral images and videos as described in:

[1] Kercher et al. _"Video-rate hyperspectral unmixing for multiplexed microscopy and microendoscopy"_. Sci Rep. (Accepted). 

[2] Harman et al. _"Denoising multiplexed microscopy images in n-dimensional spectral space"_. Biomed Opt Exp 2022; 13(8), 4298–4309.

* [Overview](#overview)
* [System Requirements](#system-requirements)
* [Get Started](#get-started)
* [Supported filetypes](#supported-filetypes)
* [Authors](#authors)
* [License](#license)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)

## Overview

HyperViewer is a tool for visualization, analysis, and batch processing of hyperspectral images and videos. HyperViewer integrates intuitive visualization of hyperspectral images with quantitative unmixing via non-negative least squares fitting routine using a predefined basis spectra library. A number of NNLS configurations are employed for comparison including the seminal active-set method ([*lsqnonneg*](https://www.mathworks.com/help/matlab/ref/lsqnonneg.html) in MATLAB), the improved [fast-NNLS method](https://doi.org/10.1002/(SICI)1099-128X(199709/10)11:5%3C393::AID-CEM483%3E3.0.CO;2-L), and the GPU-FNNLS method developed by the authors. To use GPU-FNNLS, you must have a CUDA-enabled NVIDIA GPU in a supported environment. If a compatible GPU is not available, CPU based unmixing is automatically implemented but will be quite slow. GPU-acceleration is highly recommended for large data sets.

While a number of basic App Designer features have been added (and more are in development), it is by no means complete. HyperViewer is best used within the MATLAB programming environment such that custom processing routines and support for more filetypes may be added at the users discretion. Therefore, we recommend users have a basic working knowledge of MATLAB. We encourage users to submit bug reports or feature requests for continued development of a more general hyperspectral app.

## System Requirements

### HyperViewer GUI:

This code requires MATLAB. Therefore, we recommend using systems that meet MATLAB's systems requirements. In general this means:
* OS that supports MATLAB ([Windows](https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2018b-windows.pdf), [macOS](https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2018b-macintosh.pdf), [Linux](https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2018b-linux.pdf))
* 64-bit Intel or AMD processor
* 4-6 GB Hard drive space (for new MATLAB installs)
* ≥4 GB RAM

This software was developed and tested using MATLAB 2020b (ver 9.9). The following MATLAB Toolboxes are required. Version numbers indicate releases under which this code was developed and tested.

* Image Processing Toolbox (ver 11.2)
* Statistics and Machine Learning Toolbox (ver 12.0)

The following toolboxes are recommended, but not required, for standard app operation:

* [MATLAB Parallel Computing Toolbox](https://www.mathworks.com/support/requirements/parallel-computing-toolbox.html) (ver 7.3)

### GPU-accelerated NNLS:
The following is required to run GPU-FNNLS using MATLAB:

* [MATLAB Parallel Computing Toolbox](https://www.mathworks.com/support/requirements/parallel-computing-toolbox.html) (ver 7.3)
* CUDA-enabled NVIDIA GPU with compute capability 3.0 or higher.
* Download the latest NVIDIA [graphics driver](https://www.nvidia.com/Download/index.aspx)

## Get Started
This code requires MATLAB, download the latest version [here](https://www.mathworks.com/products/matlab.html). New users of MATLAB will need a [MathWorks account](https://www.mathworks.com/mwaccount/) and to associate that account with a MATLAB license. Install times may take up to 1-2 hours depending on internet speed, but are usually faster.

Once installed, confirm MATLAB and toolbox versions are up to date in the MATLAB command window with:
```
>> ver
```

To check that a compatible GPU is installed, run:
```
>> gpuDevice
```
If a compatible device is found and drivers are correctly installed, information about your GPU will output to the command window. Otherwise, this will throw an error.

To use HyperViewer, download or clone this repository to your machine. In MATLAB, navigate to the directory containing ```unmixUI.m``` and add the functions library to the MATLAB path.
```
>> addpath(genpath([cd filesep 'lib']));
>> savepath;
```
To test that GPU-FNNLS executes without error, run:
```
>> nnlstest
```
Specific driver versions and GPU benchmark results are described in the manuscript.

To start HyperViewer, run
```
>> unmixUI
```
in the command window, or open ```unmixUI.m``` in the editor and click "Run" in the menu bar. To view or edit the layout in GUIDE, run
```
>> guide unmixUI
```

A detailed walkthrough of the HyperViewer App functions can be found in the [User Guide](www.github.com/hyperviewer_release/USERGUIDE.md).

## Demo

To see a few examples of hyperspectral unmixing, download the demo package (1.4 GB) from the [Spring Lab website](www.springlabnu.com/software). Unzip the file and move into the main HyperViewer working directory. Run the app, then click 'Demo' in the menu bar and select one of the demos.

#### Apples

Analyze a scene of real and plastic apples taken using the [Specim IQ](https://www.specim.fi/iq/) camera. Expected unmixing time is ~20 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Supplementary Figure 6.

#### Wide-field Fluorescence Imaging of Ovarian Cancer

A mouse model of advanced-stage disseminated ovarian cancer was treated with an anti-EGFR antibody with a fluorescent label and imaged using the CRi Maestro small animal imager. Expected unmixing time is ~20 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Supplementary Figure 5.

#### SMIRC - *In vivo*

A mouse model of advanced-stage disseminated ovarian cancer was imaged via 5-color SMIRC. Expected unmixing time is ~15 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Figure 1d.

#### SMIRC - Dyes in Solution

A set of Alexa Fluor Dyes were dissolved in solution and imaged via SMIRC. Expected unmixing time is ~5-10 seconds per image using a typical CPU. Default settings will reproduce data in Kercher EM et al. Figure 1a using the image tiling feature.

## Supported Filetypes
The following hyperspectral image filetypes are supported
* ENVI
  * Requires ```.dat``` and ```.hdr``` files
* SMIRC tiff stack
  * Automatically recognized by looking for laser power 'mW' string in filename
* CRi Maestro tiff stack
  * Automatically recognized by looking for exposure time 'ms' string in filename
* More to come

## Authors
* **Eric Kercher** - Primary author
* **Becca Harman** - Secondary author
* **Ryan Lang** - Secondary author
* **Paige Leven** - Secondary author
* **Ji Tae Park** - Secondary author
* **Liam Price** - Secondary author
* **Kai Zhang** - Secondary author

## License
This software is distributed under the GNU General Public License version 3 (GPLv3). If you choose to utilize this software in your work, we kindly request that you acknowledge and cite the following in your publication. Your citation is greatly appreciated and helps support ongoing development. 

[1] Kercher et al. _"Video-rate hyperspectral unmixing for multiplexed microscopy and microendoscopy"_. Sci Rep. (Accepted). 

[2] Harman et al. _"Denoising multiplexed microscopy images in n-dimensional spectral space"_. Biomed Opt Exp 2022; 13(8), 4298–4309.

## Acknowledgments
This code includes a number of third-party files from the MATLAB Central File Exchange.
* [**Export Fig**](https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig) - *Yair Altman*
* [**ENVI File reader/writer**](https://www.mathworks.com/matlabcentral/fileexchange/27172-envi-file-reader-writer) - *Felix Totir*
* [**Natural-Order Filename Sort**](https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort) - *Stephen Cobeldick*
* [**Graphical data selction tool**](https://www.mathworks.com/matlabcentral/fileexchange/13857-graphical-data-selection-tool) - *John D'Errico*
* [**nestedSortStruct**](https://www.mathworks.com/matlabcentral/fileexchange/28573-nestedsortstruct) - *Jake Hughey*
* [**Spectral and XYZ Color Functions**](https://www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions) - *Jeff Mather*
* [**RGB triple of color name, version 2**](https://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2) - *Kristjan Jonasson*

