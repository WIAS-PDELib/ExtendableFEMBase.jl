# About the Examples

The provided examples are designed with the following goals in mind:
- They can be run directly from the Julia REPL.
- Each example is implemented as a Julia module, named similarly to the file's basename.
- Examples can serve as starting points for your own projects.
- Some examples include test cases that are integrated into the test suite.
- Assembly and solve times (especially for the first run) may be significantly higher due to Julia's just-in-time compilation.


## Running the examples

In order to run `ExampleXXX`, perform the following steps:

- Download the example file (e.g., using the source code link at the top of the page).
- Ensure all required packages are installed in your Julia environment.
- In the Julia REPL, run:

```julia
julia> include("ExampleXXX.jl")
julia> ExampleXXX.main()
```

- Some examples provide visual output via the optional argument `Plotter = PyPlot` or `Plotter = GLMakie` (these require the corresponding plotting package to be installed and loaded).
