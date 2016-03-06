SpaceIntersection
==========

##Description
This repository is an implementation of space intersection, which can get the position of an object point with given information.
The necessary information include:  
+ focal length [mm]  
+ image coordinates [mm]  
+ exterior orientation parameters of two photos with overlap area [m] [deg]  


##Usage
First, create an input file. The format of an input file must like:
```
<f>
<PhoL-Name> <XL> <YL> <ZL> <OmegaL> <PhiL> <KappaL> <XL-Err> <YL-Err> <ZL-Err> <OmegaL-Err> <PhiL-Err> <KappaL-Err>
<PhoR-Name> <XR> <YR> <ZR> <OmegaR> <PhiR> <KappaR> <XR-Err> <YR-Err> <ZR-Err> <OmegaR-Err> <PhiR-Err> <KappaR-Err>
<P-Name> <xL> <yL> <xR> <yR>
```
`<f>` is focal length of camera.  
`<x*> <y*>` stands for the image coordinates of target point in left or right photo.  
`<XR> <YR> <ZR> <OmegaR> <PhiR> <KappaR>` are the exterior orientation parameters of left photo or right photo.  
Since the exterior orientation parameters are treated as observables with uncertainty, so they must have errors `<XR-Err> <YR-Err> <ZR-Err> <OmegaR-Err> <PhiR-Err> <KappaR-Err>.`

Then you can just call `./spaceIntersection.py -i <input file>` to start the computation.  
You can also type `./spaceIntersection.py -h` for more information about this repository.  
There are already two input files serve as an example.


##Requirements

###Python
[Python v2.7.X](https://www.python.org) with the following modules to be installed.

-[Numpy](http://www.numpy.org)  
-[Sympy](http://www.sympy.org/en/index.html)  
