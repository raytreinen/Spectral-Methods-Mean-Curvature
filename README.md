# Spectral-Methods-Mean-Curvature

This file corresponds to the paper "Spectral methods for prescribed mean curvature equations" by Jonas Haug, Rachel Jewell, and Ray Treinen.

General information on the programs:
We treat minimal surface, constant mean curvature surface, and capillary surface problems with Dirichlet boundary conditions and the nonlinear Neumann conditions known as capillary data on the domains of the rectangle, disk, and annulus.
The files are organized by a naming convention.  The following examples of names are sufficient to explain the codes in this directory.  MinDirAnn.m is a solver for a minimal surface with Dirichlet boundary data on the annulus.  CapillaryNeuRectangle.m is a solver for a capillary surface with nonlinear Neumann data on the rectangle.  CMCDirDisk.m is a solver for a CMC surface with Dirichlet data on the disk.

Reproducability:
Fig 1: To reproduce, modify MinDirAnn.m. Set the height at radius to be constant, ha = 1.28792 (line 16), and the height at radius b to be constant, hb = 0 (line 17). 

Fig 2: To reproduce, run CapillaryNeuRectangle.m as is. Physical parameters are the domain (lines 14-17), kappa (line 19), and boundary data (line 21). Computational parameters are the number of points (line 26).

Fig 3 top: To reproduce modify MinDirRectangle.m.  Use boundary data
g = @(x,y) 0.13*(sin(4*pi*x)).^2 + 0.1*(sin(4*pi*y));
on line 20.  On line 23 use N = 86 for the number of points to get the contour plot to smooth out.  The rest of the computational parameters are the same.

Fig 3 middle: To reproduce, run MinDirDisk.m as is.  Physical parameters are the boundary data on line 11 and the computational parameters are the number of points and the tolerances on lines 14-22.

Fig 3 bottom: To reproduce, run MinDirAnn.m as is.  Physical parameters are the boundary data on lines 17 and 18 involving trigonometric functions, and the previous lines should be commented out.   The computational parameters are the number of points and the tolerances on lines 21-27.

Fig 4 both left and right: To reproduce modify MinDirRectangle.m.  Use boundary data
g = @(x,y) 0.13*(sin(4*pi*x)).^2 + 0.1*(sin(4*pi*y));
on line 20.  On line 23 use N = 55 for the number of points.  Uncomment lines 175-184, 242-252, and 257-267 to get the error plots.  The absolute error plot may not be perfectly replicated depending on the machine it is run on.

Fig 5 top: To reproduce, modify CMCDirRectangle.m. Redefine the domain by ccc = -2 (line 16) and ddd = 2 (line 17). Physical parameters are the domain (lines 14-17), kappa (line 19), and boundary data (line 21). Computational parameters are the number of points (line 24).

Fig 5 middle: To reproduce, run CMCDirDisk.m as is. Physical parameters are lambda (line 12) and boundary data (line 14). Computational parameters are number of radial points (line 18) and number of angular points (line 20).

Fig 5 bottom: To reproduce, run CMCDirAnn.m as is. The physical parameters are lambda, the inner and outer radii, and the inner and outer boundary functions. Lambda is set to 0.5 (line 16). For the inner radius, a = 1 (line 13). For the outer radius, b = 2 (line 14). For the inner boundary function, ha(t) = 0.5 + 0.1*sin(2*t).^2 (line 18). For the outer boundary function, hb(t) = 0.5 - 0.1*cos(2*t).^2 (line 19). The computational parameters are the number of radial and angular grid points, the Newton and BVP tolerances, epsilon, and MM. For the radial grid points, N = 50 (line 22). For the angular grid points, M1 = 80 (line 23). For the Newtown tolerance, new_tol = 1e-13 (line 25). For the BVP tolerance, bvp_tol = 1e-10 (line 26). For epsilon, ep = 1e-8 (line 27). MM = 100 (line 28).

Fig 6 top: To reproduce, run capillaryDirRectangle.m as is. Physical parameters are the domain (lines 14-17), kappa (line 19), and boundary data (line 21). Computational parameters are the number of points (line 24).

Fig 6 middle: To reproduce, run CapillaryDirDisk.m as is. Physical parameters are kappa (line 12) and boundary data (line 13). Computational parameters are number of radial points (line 16) and number of angular points (line 18).

Fig 6 bottom: To reproduce, run CapillaryDirAnn.m as is. Physical parameters are radii (lines 13 and 14), kappa (line 16), and boundary data (lines 18 and 19). Computational parameters are number of radial points (line 22) and number of angular points (line 23).

Fig 7 top: To reproduce, run MinNeuDisk.m as is. Physical parameters are kappa (line 12) and boundary data (line 13). Computational parameters are number of radial points (line 16) and number of angular points (line 18).

Fig 7 bottom: To reproduce, run MinNeuAnn.m as is. Physical parameters are radii (lines 14 and 15) and boundary data (lines 18 and 19). Computational parameters are number of radial points (line 23) and number of angular points (line 24).

Fig 8 top: To reproduce, run CMCNeuDisk.m as is. Physical parameters are boundary data (line 13). Computational parameters are number of radial points (line 16) and number of angular points (line 18).

Fig 8 bottom: To reproduce, run CMCNeuAnn.m as is. Physical parameters are radii (lines 14 and 15) and boundary data (lines 18 and 19). Computational parameters are number of radial points (line 23) and number of angular points (line 24).

Fig 9 top: Modify CapillaryNeuRectangle.m and set ccc = -10 (line 16) and ddd = 10 (line 17).  This may take some time and require a faster machine.  On a 2023 Macbook Pro it takes approximately 12 minutes.  Use “axis square” in the plotting  commands on lines 336 and 344.

Fig 9 mid: To reproduce, run CapillaryNeuDisk.m as is. Physical parameters are kappa (line 12), and boundary data (line 14). Computational parameters are number of radial points (line 17) and number of angular points (line 19).

Fig 9 bottom: To reproduce, modify CapillaryNeuAnn.m. Set the boundary data at radius b to be bgamma = @(t) pi/3 + 0.2*cos(6*t) (line 19). Physical parameters are radii (lines 13 and 14), kappa (line 16), and boundary data (lines 18 and 19). Computational parameters are number of radial points (line 22) and number of angular points (line 23).

Fig 10 top: To reproduce, run MinMixAnn.m as is. Physical parameters are radii (lines 15 and 16), and boundary data (lines 18-24). Computational parameters are number of radial points (line 27) and number of angular points (line 28).

Fig 10 mid: To reproduce, run CMCMixAnn.m as is. Physical parameters are radii (lines 16 and 17), lambda (line 19), and boundary data (lines 21-27). Computational parameters are number of radial points (line 30) and number of angular points (line 31).

Fig 10 bottom: To reproduce, run CapillaryMixAnn.m as is. Physical parameters are radii (lines 16 and 17), kappa (line 19), and boundary data (lines 21-27). Computational parameters are number of radial points (line 30) and number of angular points (line 31).
![image](https://github.com/raytreinen/Spectral-Methods-Mean-Curvature/assets/105992653/da028727-c9be-46b7-a030-9eea94b9a4db)
