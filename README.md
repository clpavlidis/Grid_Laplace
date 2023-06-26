# Grid Generation

This code represents the generation of a body fitted, structured computational grid, in 2D space. The boundary nodes of the computational grid were created using elliptic-type partial differential equations (PDEs), while for the inner nodes, the method of 2D binomial interpolation was used. The fundamental 2D grid generation PDEs, were described by the Laplacian form. For the correctness of the code, a random internal node was chosen to be positioned outside the computational grid in order to restore it to its original location, via the iterative algorithm Gauss-Seidel (G-S). For the modeled case study, an appropriate code in FORTRAN90 programming language was developed.

## Prerequisites

- Appropriate FORTRAN90 compiler

- Excellent knowledge in Numerical Methods

- Applied knowledge in Partial Differential Equations (PDEs)

- Applied knowledge in Grid Generation Methods (Grid Metrices, etc.)

## Usage

- The numbering of the nodes of the computational grid, follows the string numbering starting from the left side and from bottom to top (i.e. the 1st node of the computational grid, is placed in the bottom left side)

- Due to the mixed derivative terms on the metrices of the computational grid, it results a 9-node stencil. Using clockwise rotation, the stencil's corner nodes are flipped to the main stencil cross. Thus, the node stencil from 9, is transformed to 5. For the node flipping, backward differentiation is selected between the central node of each side an its previous, considering zero RHS (Right Hand Side)

- All the 1st order PDEs, are modelled using central differentiation

- IM & JM are the total nodes on the x- & the y-direction correspondingly

- Code lines 69 & 70, describe the coordinates of the boundary nodes of the computational grid (the user is able to modify them as its will)

- Code lines 79 & 84, describe the coordinates of the inner nodes of the computational grid using 2D binomial interpolation (the user is able to modify them as its will)

- Code lines 91 & 92, describe the selected node and its coorsinates in order to be placed outside the computational grid (the user is able to modify them as its will)

- d1,...,d5, are the diagonals of the 5-node stencil (Penta-diagonal matrix formation)

## Compile-Run

Compile source code:

```bash
gfotran ./Grid_Laplace.f90 -o Grid_Laplace.exe
```

Run the generated executable:

```bash
chmod +x ./Grid_Laplace.exe
./Grid_Laplace.exe
```

## Output Files

- Initial computational grid (with the selected node, placed outside the computational grid)

- Boundary nodes coordinates

- Jacobian matrix values (area of each cell of the computational grid)

- Contravariant metrices

- Final computational grid (the selected node, is restored to its original position inside the computational grid)

- Min & max coordinates of the computational grid grid on each direction (computational grid length)

- Norm of the x-direction, after the solver is being converged

- Norm of the y-direction, after the solver is being converged

- Informational data, about the code provess (CPU needed time, etc.)

## Contact Data

- Author: Chariton L. Pavlidis (@clpavlidis), aka "Bogdan"
- Contact: clpavlidis@outlook.com
