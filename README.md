SpaceIntersection
==========

## Description
This repository is an implementation of space intersection, which can get the position of an object point with given information.
The necessary information include:  
+ focal length [mm]  
+ image coordinates [mm]  
+ exterior orientation parameters of more than two photos with overlap area [m] [deg]  


## Usage
First, create an input file. The format of an input file must like:
```
<f>
<Photo-Name> <x1> <y1> <XL1> <YL1> <ZL1> <Omega1> <Phi1> <Kappa1> <XL1-Err> <YL1-Err> <ZL1-Err> <OmegaL1-Err> <PhiL1-Err> <KappaL1-Err>
<Photo-Name> <x2> <y2> <XL2> <YL2> <ZL2> <Omega2> <Phi2> <Kappa2> <XL2-Err> <YL2-Err> <ZL2-Err> <OmegaL2-Err> <PhiL2-Err> <KappaL2-Err>
<Photo-Name> <x3> <y3> <XL3> <YL3> <ZL3> <Omega3> <Phi3> <Kappa3> <XL3-Err> <YL3-Err> <ZL3-Err> <OmegaL3-Err> <PhiL3-Err> <KappaL3-Err>
```
`<f>` is focal length of camera.  
`<x*> <y*>` stands for the image coordinates of target point.  
`<XL*> <YL*> <ZL*> <Omega*> <Phi*> <Kappa*>` are the exterior orientation parameters.  
Since the exterior orientation parameters are treated as observables with uncertainty, so they must have errors `<XL*-Err> <YL*-Err> <ZL*-Err> <Omega*-Err> <Phi*-Err> <Kappa*-Err>.`

## Requirements

### Python
[Python 3 ](https://www.python.org) with the following modules to be installed.

-[Numpy](http://www.numpy.org)  
-[Sympy](http://www.sympy.org/en/index.html)  
-[Pandas](http://pandas.pydata.org/)  
