Photogrammetry Space Intersection
==========

## _Description_ 
This repository implements photogrammetry space intersection, getting the position of an object point with given information.
The necessary information includes:  
+ focal length f [mm]  
+ image coordinates x y [mm]  
+ exterior orientation parameters (XL, YL, ZL, Omega, Phi, Kappa) of two or more overlapped photos with overlap area [meter] [deg]  


## _Usage_
First, you'll need to create an input file. The format of an input file must like this:
```
<f>
<Photo-Name> <x1> <y1> <XL1> <YL1> <ZL1> <Omega1> <Phi1> <Kappa1> <XL1-Err> <YL1-Err> <ZL1-Err> <OmegaL1-Err> <PhiL1-Err> <KappaL1-Err>
<Photo-Name> <x2> <y2> <XL2> <YL2> <ZL2> <Omega2> <Phi2> <Kappa2> <XL2-Err> <YL2-Err> <ZL2-Err> <OmegaL2-Err> <PhiL2-Err> <KappaL2-Err>
<Photo-Name> <x3> <y3> <XL3> <YL3> <ZL3> <Omega3> <Phi3> <Kappa3> <XL3-Err> <YL3-Err> <ZL3-Err> <OmegaL3-Err> <PhiL3-Err> <KappaL3-Err>
```
_Note_: if the exterior orientation parameters do not include any errors, you can use values near zero for example, 0.0001

`<f>` is the focal length of the camera.  
`<x*> <y*>` stands for the image coordinates of the target point.  
`<XL*> <YL*> <ZL*> <Omega*> <Phi*> <Kappa*>` are the exterior orientation parameters.  
Since the exterior orientation parameters are treated as observables with uncertainty, they must have errors `<XL*-Err> <YL*-Err> <ZL*-Err> <Omega*-Err> <Phi*-Err> <Kappa*-Err>.`

## _Requirements_

### _Python_
[Python 3 ](https://www.python.org) with the following modules to be installed.

-[Numpy]
-[Sympy]
-[Pandas]
