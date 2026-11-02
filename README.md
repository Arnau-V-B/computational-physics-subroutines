# Computational physics subroutines

This is a small compilation of subroutines written in **Fortran 90** useful for numerically solving common physics problems. In particular, they were used for the lab sessions in the *Computational Physics* course imparted during the [Bachelor's degree in Physics at Universitat de Barcelona](https://web.ub.edu/en/web/estudis/w/bachelordegree-g1035).

## Content

- **[`function_roots_and_interpolation.f90`](./function_roots_and_interpolation.f90) (Ch. 2)**:

    Contains subroutines for zeroth-order and linear interpolation (`interpol_0` and `interpol_1`), as well as root-finding methods: bisection (`bisection`), Regula Falsi (`regula_falsi`) and Newton–Raphson (`newtonraphson`). The file also includes an example function and derivative used by the root-finding routines; adapt it to the function you want to study.

- **[`numerical_integration.f90`](./numerical_integration.f90) (Ch. 3)**:

    Contains subroutines to numerically integrate any given function. Specifically, it provides code for the chained trapezoidal rule (`trapezoidalrule`), the Simpson's 3/8 rule (`simpson38`) and a subroutine to derivate tabular data up to second-order accuracy (`tabderivative`).

- **[`random_numbers_and_montecarlo.f90`](./random_numbers_and_montecarlo.f90) (Ch. 4)**:

    Contains subroutines to generate random number distributions and to perform Monte Carlo integrations. It includes a histogram generator (`histogram`), sampling routines (acceptance-rejection `accept`, uniform `uniform`, exponential `exponential` and Gaussian via Box-Müller `gaussian`), generalized Monte Carlo integration (`montecarlo`) and hit-or-miss integration (`hitmiss`).

- **[`ordinary_differential_equations.f90`](./ordinary_differential_equations.f90) (Ch. 5)**:

    Contains subroutines to solve ordinary differential equations (ODEs) of the form $d^2y = f(x, y, dy)$. It includes a sample ODE structure (`odefunc`), one-point integrators (Euler 1st `euler_first` and 2nd order `euler_second`, and Runge-Kutta 4th order `rungekutta4order`), two-point integrators (Adams-Bashforth `adams_bashforth_multi`) and a shooting method (`shooting`) as a use case example to find eigenvalues $E_n$ to the stationary Schrödinger equation with the secant method.

- **[`partial_differential_equations.f90`](./partial_differential_equations.f90) (Ch. 6)**:

    Contains subroutines to solve partial differential equations (PDEs). It includes steady-state solvers for the 2D Poisson equation with Dirichlet (`poisson_2D_dir`) and Neumann (`poisson_2D_neu`) boundary conditions (using Jacobi, over-relaxation or Gauss-Seidel schemes), a 1D Poisson solver with losses (`poisson_1D_loss`), the Crank-Nicolson method for the 1D heat equation with losses (`crank_nicolson`) and a tridiagonal system solver (`tridiag`) used internally by Crank-Nicolson.

## Reminder: using Fortran 90

The `.f90` files in this repository contain subroutines and supporting procedures, not standalone programs with a `program` entry point. To run them, write a **main program** that calls the desired routines and compile it together with the relevant source file(s).

### Typical Fortran 90 script structure

```f90
program main
    implicit none

    ! Variable declaration
    integer :: i, n
    double precision :: a, b, x, y, res
    double precision :: my_func
    ...

    ! Common blocks (if needed by the subroutines)
    COMMON/PARAMETER/w
    external my_func
    ...

    ! Initialisation
    a = 0.d0
    b = 1.d0
    n = 100
    ...

    ! Code lines:
    ! call the subroutines
    call trapezoidalrule(a, b, n, myfunc, res)
    ...
    ! call functions
    y = my_func(x)
    ...

    ! Output
    write(*,*) "Integral =", res
    ...

end program main

! Definition of functions
double precision function myfunc(x)
    implicit none
    double precision, intent(in) :: x

    myfunc = x**2
    return
end function myfunc
...

! Definition of subroutines
subroutine trapezoidalrule(x1, x2, k, func, res)
    implicit none
    double precision :: x1, x2, func, res
    integer :: k

    ...

    return
end subroutine trapezoidalrule
```

### Compilation with `gfortran`

In order to execute a `.f90` script first it needs to be **compiled** into a `.exe` file by a dedicated compiler. For example, when using the [GNU Fortran compiler](https://fortran-lang.org/learn/os_setup/install_gfortran/) (`gfortran`), you will need to write the following lines in the terminal in order to **compile** and then **run** `main`:

```bash
gfortran main.f90 -o main
./main
```