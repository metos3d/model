Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D

-- Model --

Metos3D can be coupled to every (biogeochemical) model that conforms to the following interface:

subroutine metos3dbgc(n, nz, m, nb, nd, dt, q, t, y, u, b, d)
    integer :: n           ! tracer count
    integer :: nz          ! layer count
    integer :: m           ! parameter count
    integer :: nb          ! boundary condition count
    integer :: nd          ! domain condition count
    real*8  :: dt          ! ocean time step
    real*8  :: q(nz, n)    ! bgc model output
    real*8  :: t           ! point in time
    real*8  :: y(nz, n)    ! bgc model input
    real*8  :: u(m)        ! parameters
    real*8  :: b(nb)       ! boundary conditions
    real*8  :: d(nz, nd)   ! domain conditions
end subroutine

The interface decouples biogeochemical models and driver routines (ocean circulation, forcing, geometry) programmatically. It gives you the possibility to provide a free number of tracers, parameters, boundary and domain conditions. It suits well an optimization as well as an Automatic Differentiation (AD) context.

-- Documentation --

See the LaTeX document in the 'doc' directory for examples and more detailed information.
