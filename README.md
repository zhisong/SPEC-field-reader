# SPEC-field-reader
This repository contains FORTRAN modules that reads the SPEC input and computes the magnetic field and metrics.

Note: the code is only tested in stellarator symmetric cases currently.

### SPEC geometry and field
Please see the following two links for references.

For magnetic field representation: [matrix.pdf](https://w3.pppl.gov/~shudson/Spec/matrix.pdf)

For geometry representation: [coords.pdf](https://w3.pppl.gov/~shudson/Spec/coords.pdf)

### Usage
1. Compile and include all the modules, declare a type(state) to carry the information from SPEC output.
```fortran
use spec_state
use spec_field
use spec_geometry
use spec_io
type(state) :: ss
```

2. Read the SPEC output, the .sp.h5 file and the .sp.A file. For example:
```fortran
call read_spec_h5('G3V02L1Fi.001.sp.h5',ss)
call read_spec_field('.G3V02L1Fi.001.sp.A',ss)
```

3. Compute the field or metrics at a desired location (s, theta, xi).
```fortran
integer :: lvol
real  :: s, theta, xi
real  :: a(3), gb(3), dgb(3,3)
real  :: jac, djac(3), x(3), gij(3,3), dgij(3,3,3)

lvol =  2  ! in volume 2
s =  -0.3
theta =  0.7
xi =  2.3
call get_spec_field(ss%A(lvol), ss%lrad(lvol), s, theta, xi, a, gb, dgb)
call get_spec_coord(ss%Ri, lvol, ss%mn, s, theta, xi, jac, djac, x, gij, dgij)
```
Please note that the output "gb" has three components, and it is actually $J B^i$, with "dgb" its derivatives. The output "a" is actually $A_i$. The output "jac" is the Jacobian $J$, "djac" is its derivatives, "x" is the coordinates, and the metric tensors "gij" are actually $g_{ij}$. Finally, "dgij(i,j,k)" is the partial derivatives of gij, i.e. $\partial_{x^k} g_{ij}$, where $(x^1, x^2, x^3) = (s, \theta, \xi)$.

The test example "G3V02L1Fi.001" is an elliptic stellarator case (fixed-boundary), one of the standard test cases of SPEC. The case has two inner volumes. The outside boundary of the plasma is shown in the figure below, while the Poincare cross section at $0, \frac{\pi}{10}, \frac{2\pi}{10}, \frac{3\pi}{10}$. 

![boundary](/images/example_boundary.png) 
![poincare0](/images/example_poincare0.png)
![poincare1](/images/example_poincare1.png)
![poincare2](/images/example_poincare2.png)
![poincare3](/images/example_poincare3.png)
