# Meshing

Meshes in ExtendableFEMBase are represented using the `ExtendableGrid` type. For detailed information on grid structures and constructors, refer to the [ExtendableGrids.jl documentation](https://github.com/j-fu/ExtendableGrids.jl).

To generate simplex grids, you can use external tools such as [SimplexGridFactory.jl](https://github.com/j-fu/SimplexGridFactory.jl).

Each cell, face, and edge in a mesh is associated with an `AbstractElementGeometry` (as defined in [ExtendableGrids.jl](https://github.com/j-fu/ExtendableGrids.jl)). These geometries are used to dispatch key functionality, including local-to-global transformations, enumeration rules, basis function definitions, volume calculations, and mesh refinements. See below for a list of supported element geometries.


## Recognized Geometries and Reference Domains

The following list contains all subtypes of ExtendableGrids.AbstractElementGeometries and their reference domains for which the package offers finite elements on them.

----
##### Edge1D <: AbstractElementGeometry1D

    [1]-----[2]               [1] = [0]
                              [2] = [1]
                              
----
##### Triangle2D
    
    [3]                 
     | \   
     |   \                    [1] = [0,0]
     |     \                  [2] = [1,0]
     |       \                [3] = [0,1]
     |         \ 
    [1]--------[2]
            
----
##### Parallelogram2D <: Quadrilateral2D

    [4]--------[3]               
     |          |             [1] = [0,0]
     |          |             [2] = [1,0]
     |          |             [3] = [1,1]
     |          |             [4] = [0,1]
    [1]--------[2]

    Note: most finite elements only work as intended on Parallelogram2D
          since the local<>global map stays affine in this case


----
##### Tetrahedron3D

    [4]                 
     |\\   
     | \ \                    [1] = [0,0,0]
     |  \  \                  [2] = [1,0,0]
     |   \   \                [3] = [0,1,0]
     | _-[3]-_ \              [4] = [0,0,1]
    [1]--------[2]


----
##### Parallelepiped3D <: Hexahedron3D
                         
        [8]--------[7]        [1] = [0,0,0]
       / |        / |         [2] = [1,0,0]
    [5]--------[6]  |         [3] = [1,1,0]
     |   |      |   |         [4] = [0,1,0]
     |   |      |   |         [5] = [0,0,1]
     |  [4]-----|--[3]        [6] = [1,0,1]
     | /        | /           [7] = [1,1,1]
    [1]--------[2]            [8] = [0,1,1]

    Note: most finite elements only work as intended on Parallelepiped3D
          since the local<>global map stays affine in this case
