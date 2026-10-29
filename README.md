<img width="21654" height="10454" alt="fig5" src="https://github.com/user-attachments/assets/4d3868b1-e855-40a5-a41e-01048c82b3a5" /># LSM-MGM
The code for the fast implement of the paper "An Improved Neighborhood-based Mismatch Removal Algorithm for Consumer Device". Please cite the following paper if you use the code in this repository as part of a published research project.
Successfully test on Matlab R2019b.

## An Improved Neighborhood-Based Mismatch Removal Algorithm for Consumer Device
<br>Haotian Zhang , Wenlong Hou, Ning Jia , and Weidong Zhao<br>
[DOI 10.1109/TCE.2024.3409774](https://ieeexplore.ieee.org/document/10549964/)

## Usage
We have provided an example with tuned initialization and parameters for LSM-MGM in example.m. Directly running example.m would give the result.

## The comparative test codes
<br>LPM: https://github.com/jiayi-ma/LPM. and https://github.com/jiayi-ma/LPM-GC.
<br>GLOF: https://github.com/gwang-cv/GLOF.
<br>RFM-SCAN: https://github.com/StaRainJ/RFM-SCAN.
<br>EAS-A: https://github.com/AoxiangFan/EifficientDeterministicSearch.
<br>LOGO: https://github.com/XiaYifan1999/LOGO.
<br>GBCPD: https://github.com/ohirose/bcpd.
<br>LCT: https://github.com/gwang-cv/LCT.
<br>ConvMatch: https://github.com/SuhZhang/ConvMatch/tree/convmatch_plus.
<br>CHCANet: https://github.com/gwang-cv/CHCANet.

### Other comparison results with algorithm LSM-MGM
<p align="center">
  <img width="2929" height="3086" alt="test image1" src="data/figureImage/testImage1.jpg">
</p>
<br>The rows in the figure, from top to bottom, correspond to the image sets DAISY, DTU, VGG, HEINLY, and MegaDepth1500. The columns, from left to right, show the cumulative distributions of the inlier ratio, recall, precision, and F1 score, respectively.
<br>
<br>
<p align="center">
  <img width="2929" height="2285" alt="test image2" src="data/figureImage/testImage2.jpg">
</p>
<br>The rows in the figure, from top to bottom, correspond to the image sets RS, Print, and FIRE. The columns, from left to right, show the cumulative distributions of the inlier ratio, recall, precision, and F1 score, respectively.

## Datasets
<br>DAISY: https://www.epfl.ch/labs/cvlab/software/descriptors-and-keypoints/daisy/
<br>DTU: https://roboimagedata.compute.dtu.dk/
<br>VGG: https://www.robots.ox.ac.uk/~vgg/research/affine/
<br>HEINLY: https://www.cs.unc.edu/~jheinly/binary_descriptors.html
<br>MegaDepth1500: https://www.cs.cornell.edu/projects/megadepth/ and https://cvg-data.inf.ethz.ch/megadepth/
<br>RS: https://github.com/jiayi-ma/GLPM/tree/master
<br>FIRE: https://www.ics.forth.gr/cvrl/fire/
<br>Print: The images obtained in "Print" dataset  are sourced from a private engineering project.
<br>SMVS: https://exhibits.stanford.edu/data/catalog/rb470rw0983
