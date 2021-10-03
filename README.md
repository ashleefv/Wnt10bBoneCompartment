# Wnt10bBoneCompartment
Mathematical modeling of the effects of Wnt-10b on bone metabolism<br/>
[Preprint](https://biorxiv.org/cgi/content/short/2021.06.12.448204v1)



## Overview
This ODE model describes the effects of Wnt-10b on bone metablolism through mutilple remodeling cycles. The model consists of differential equations that track the active cell populations of osteocytes, pre-osteoblasts, osteoblasts, and osteoclasts. The model also uses a diffential equation to track bone volume chagnes. Wnt-10b is connected to the sytem through a hill function relationship. The model allows for the exploration of the dynamics through various images created by the code.

## Authors
Carley V. Cook<sup>a,b</sup>,  Mohammad Aminul Islam<sup>a,b</sup>,  Brenda J. Smith<sup>c</sup>,  Ashlee N. Ford Versypt<sup>a,b,d</sup><br/>
<sup>a</sup>School of Chemical Engineering, Oklahoma State University, Stillwater, OK, USA<br/>
<sup>b</sup>Department of Chemical and Biological Engineering, University at Buffalo, The State University of New York, Buffalo, NY,USA<br/>
<sup>c</sup>Department of Nutritional Sciences, Oklahoma State University, Stillwater, OK, USA<br/>
<sup>d</sup>Institute for Computational and Data Sciences, University at Buffalo, The State University of New York, Buffalo, NY, USA<br/>

## Scripts

* BoneCompartmentUp.m This file runs the model and utilizes lsqurve fit to estimate the parameters related to Wnt-10b. 
* GraphsforPaper.m This file runs the model and produces graphs that shows the various dynamics involved in remodeling. The code also has the option to automatically save the images. Note that the images will be overwritten with each new run with the save option on.
* WntBoneRemodeling.sbml This file provides an sbml version of the model. This version does not support multiple remodeling cycles. Currently the Wnt-10b fold change is set to 50, but this can be easily changed by the user.

## Recommended Supplementary Packages
[Filled area plot](https://www.mathworks.com/matlabcentral/fileexchange/69652-filled-area-plot) is used to create the shaded region for the validation graph.

## Acknowledgements
Research reported in this publication was supported by the National Institute of General Medical Sciences of the National Institutes of Health under award number R35GM133763. The content is solely the responsibility of the authors and does not necessarily represent the offcial views of the National Institutes of Health.
