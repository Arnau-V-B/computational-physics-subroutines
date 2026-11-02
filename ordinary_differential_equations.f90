! This script contains some useful subroutines to solve ordinary differential equations (ODE).


subroutine odefunc(neq, x, yin, dyout)
    implicit none
    integer :: neq
    double precision :: x, yin(neq), dyout(neq)
    double precision :: E, V

    ! DESCRIPTION:
    ! It receives an abscissa x and a vector of coordinates y and returns another vector with the same dimensions with the corresponding derivatives
    ! of each coordinate according to the ODE function of the problem d2y = f(x,y,dy). It allows generalizing any problem of this style since the
    ! definition of the functions is done inside it.

    ! VARIABLES:
    ! (input)
    ! neq   --> Number of dimensions of the physical problem (number of generalized coordinates in vector y)
    ! x     --> Value of the corresponding abscissa (where the derivative is to be found)
    ! yin   --> Vector with all the coordinates of the problem that are to be differentiated
    ! (output)
    ! dyout --> Vector with the differentiated coordinates y
    ! (internal)
    ! Everything else that specifically belongs to the description of the problem


    ! EXAMPLE: time-independent Schrödinger equation --> d2ψ = -2(E-V)ψ

    COMMON/ENERGY/E
    V = -2.4d0

    dyout(1) = yin(2)
    dyout(2) = -2*(E-V) * yin(1)

    return
end


subroutine euler_first(h, neqs, yin, odefunc, yout)
    implicit none
    integer :: neqs
    double precision :: h, yin(neqs), yout(neqs)
    double precision :: dy(neqs)
    double precision :: x

    COMMON/ABSCISSA/x

    ! DESCRIPTION:
    ! This subroutine executes the first-order Euler method to integrate a single step of the coordinate vector y with its derivatives dy from
    ! an ordinary differential equation of shape d2y = func(x,y,dy). It returns a vector dfouncout with the integrated values of y and dy.
    ! Scheme:
    !   K1 = f(x, y_in)
    !   y_out = y_in + h*K1

    ! VARIABLES:
    ! (input)
    ! h       --> Predefined step that the method will use
    ! neqs    --> Number of dimensions of the physical problem (number of generalized coordinates y)
    ! yin     --> Vector of neqs components with the coordinates to be integrated (y)
    ! odefunc --> Subroutine with the specific ODE of the problem that returns the derivatives of y and dy at each abscissa
    ! (output)
    ! yout    --> Vector of neqs components with the already integrated coordinates (y')
    ! (internal)
    ! x       --> Abscissa of the problem (parametric) for the coefficient K1
    ! dy      --> Vector with neqs dimensions that collects the value of the derivatives of y (i.e. K1)

    ! We call subroutine 'odefunc()' so that it computes the derivatives (dy) of the coordinate vector (yin) at the abscissa x
    call odefunc(neqs, x, yin, dy)

    ! We compute the vector of integrated coordinates (yout) with the Euler scheme
    yout = yin + h*dy

    return
end


subroutine euler_first_multi(odefunc, N, xlow, xhigh, Ndim, y0, yx)
    implicit none
    integer :: N, Ndim, i, j
    double precision :: xlow, xhigh, y0(Ndim), h
    double precision :: yx(Ndim,N+1)
    double precision :: x, y(Ndim), dy(Ndim)

    ! DESCRIPTION:
    ! This subroutine executes N steps of the first order Euler method to approximate the evolution of the coordinate vector y with its derivatives
    ! dy from an ordinary differential equation of shape d2y = func(x,y,dy). It returns a matrix yx with the evolution of y at each iteration.
    ! Scheme:
    !   K1 = f(x_n, y_n)
    !   y_n+1 = y_n + h*K1

    ! VARIABLES:
    ! (input)
    ! odefunc     --> Subroutine carrying the y --> dy transformations characteristic of the problem at each abscissa value
    ! N           --> Number of iterations of the method
    ! xlow, xhigh --> Interval of x in which the evolution of the coordinates y is to be determined
    ! Ndim        --> Dimensions of the coordinate vector y (nº of coordinates to integrate)
    ! y0          --> Vector with Ndim dimensions carrying the initial conditions of the coordinates y
    ! (output)
    ! yx          --> Matrix of dim(Ndim,N+1) that records the evolution of each coordinate y for each iteration of the method (includes y0)
    ! (internal)
    ! h           --> Step of the interval defined in x
    ! x           --> Independent variable of the problem (parametric)
    ! y           --> Vector with Ndim dimensions that collects the last computed value of the coordinates y at each iteration
    ! dy          --> Vector with Ndim dimensions that collects the last computed value of the derivatives dy at each iteration
    ! i, j        --> Counters

    ! We determine the equally spaced step with which the method will advance
    h = abs(xhigh - xlow) / dble(N)

    ! We introduce the initial conditions y0 into the final list matrix for each component y
    do j = 1, Ndim
        yx(j,1) = y0(j)
    enddo

    ! We initialize y with the initial conditions
    y = y0

    ! We execute the Euler method N times
    do i = 1, N
        ! We generate the x variable of the previous step
        x = xlow + (i-1)*h

        ! We call subroutine "odefun()" so that it computes the derivatives dy from the abscissa x and the coordinates y of the previous step
        call odefunc(Ndim, x, y, dy)

        ! In a loop that goes through each of the coordinates we compute its next value using the Euler scheme and record it in the yx matrix
        do j = 1, Ndim
            yx(j,i+1) = yx(j,i) + h*dy(j)
            ! We update the coordinates y for the next iteration/step
            y(j) = yx(j,i+1)
        enddo

    enddo

    return
end


subroutine euler_second(h, neqs, yin, odefunc, yout)
    implicit none
    integer :: neqs
    double precision :: h, yin(neqs), yout(neqs)
    double precision :: dy1(neqs), y2(neqs), dy2(neqs)
    double precision :: x, x2

    COMMON/ABSCISSA/x

    ! DESCRIPTION:
    ! This subroutine executes the second-order Euler method to integrate a single step of the coordinate vector y with its derivatives dy from an
    ! ordinary differential equation of shape d2y = func(x,y,dy). It returns a vector yout with the integrated values of y and dy.
    ! Scheme:
    !   K1 = f(x, y_in)
    !   K2 = f(x + h/2, y_in + h/2*K1)
    !   y_out = y_in + h*K2

    ! VARIABLES:
    ! (input)
    ! h        --> Predefined step that the method will use
    ! neqs     --> Number of dimensions of the physical problem (number of generalized coordinates y)
    ! yin      --> Vector of neqs components with the coordinates to be integrated (y)
    ! odefunc  --> Subroutine with the specific ODE of the problem that returns the derivatives of y and dy at each abscissa
    ! (output)
    ! yout     --> Vector of neqs components with the already integrated coordinates (y')
    ! (internal)
    ! x, x2    --> Abscissas of the problem (parametric) for the coefficients k1 and k2
    ! dy1, dy2 --> Vectors with neqs dimensions that collect the value of the derivatives of y1 (i.e. K1) and y2 (i.e. K2) respectively

    ! We call subroutine 'odefunc()' so that it computes the derivatives (dy1) of the coordinate vector (yin) at the abscissa x
    call odefunc(neqs, x, yin, dy1)

    ! We compute the modified abscissa and coordinates for the K2 parameter
    x2 = x + h/2.d0
    y2 = yin + h/2.d0 * dy1

    ! We call 'odefunc()' so that it returns K2 from x2 and y2
    call odefunc(neqs, x2, y2, dy2)

    ! We compute the vector of integrated coordinates (yout) with the Euler algorithm
    yout = yin + h*dy2

    return
end


subroutine euler_second_multi(odefunc, N, xlow, xhigh, Ndim, y0, yx)
    implicit none
    integer :: N, Ndim, i, j
    double precision :: xlow, xhigh, y0(Ndim), h
    double precision :: yx(Ndim,N+1)
    double precision :: x1, x2, y1(Ndim), y2(Ndim), dy1(Ndim), dy2(Ndim)

    ! DESCRIPTION:
    ! This subroutine executes N steps of the second-order Euler method to approximate the evolution of the coordinate vector y with its derivatives
    ! dy from an ordinary differential equation of shape d2y = func(x,y,dy). It returns a matrix yx with the evolution of y at each iteration.
    ! Scheme:
    !   K1 = f(x_n, y_n)
    !   K2 = f(x_n + h/2, y_n + h/2*K1)
    !   y_n+1 = y_n + h*K2

    ! VARIABLES:
    ! (input)
    ! odefunc     --> Subroutine carrying the y-->dy transformations characteristic of the problem at each abscissa value
    ! N           --> Number of iterations of the method
    ! xlow, xhigh --> Interval of x in which the evolution of the coordinates y is to be determined
    ! Ndim        --> Dimensions of the coordinate vector y (nº of coordinates to integrate)
    ! y0          --> Vector with Ndim dimensions carrying the initial conditions of the coordinates y
    ! (output)
    ! yx          --> Matrix of dim(Ndim,N+1) that records the evolution of each coordinate y for each iteration of the method (includes y0)
    ! (internal)
    ! h           --> Step of the interval defined in x
    ! x1, x2      --> Independent variables of the problem (parametric) for the coefficients K1 and K2 respectively
    ! y1, y2      --> Vectors with Ndim dimensions that collect the last computed value of the coordinates y1 (for K1) and y2 (for K2) at each iteration respectively
    ! dy1, dy2    --> Vectors with Ndim dimensions that collect the last computed value of the derivatives dy1 (i.e. K1) and dy2 (i.e. K2) at each iteration respectively
    ! i, j        --> Counters

    ! We determine the equally spaced step with which the method will advance
    h = abs(xhigh - xlow) / dble(N)

    ! We introduce the initial conditions y0 into the final list matrix for each component y
    do j = 1, Ndim
        yx(j,1) = y0(j)
    enddo

    ! We initialize y1 with the initial conditions
    y1 = y0

    ! We execute the second-order Euler method N times
    do i = 1, N
        ! We generate the abscissa x1 prior to the step we want to compute
        x1 = xlow + (i-1)*h
        ! We call subroutine "odefun()" so that it computes the derivatives dy1 (i.e. K1) from the abscissa x1 and the coordinates y1 of the previous step
        call odefunc(Ndim, x1, y1, dy1)

        ! We generate the abscissa x2 from x1
        x2 = x1 + h/2.d0
        ! We modify the coordinates y1 to generate y2
        y2 = y1 + h/2.d0 * dy1
        ! We call subroutine "odefun()" so that it computes the derivatives dy2 (i.e. K2) from the abscissa x2 and the coordinates y2
        call odefunc(Ndim, x2, y2, dy2)

        ! In a loop that goes through each of the coordinates we compute by Euler its next value and record it in the yx matrix
        do j = 1, Ndim
            yx(j,i+1) = yx(j,i) + h*dy2(j)
            ! We update the coordinates y1 for the next iteration/step
            y1(j) = yx(j,i+1)
        enddo

    enddo
    return
end


subroutine adams_bashforth_multi(odefunc, N, xlow, xhigh, Ndim, y0, yx)
    implicit none
    integer :: N, Ndim, i, j
    double precision :: xlow, xhigh, y0(Ndim), h
    double precision :: yx(Ndim,N+1)
    double precision :: x1, x2, y1(Ndim), y2(Ndim), dy1(Ndim), dy2(Ndim)

    ! DESCRIPTION:
    ! This subroutine executes N steps of the multi-step Adams-Bashforth method to approximate the evolution of the coordinates of vector y with its
    ! derivatives dy, with an initial Euler step, from an ordinary differential equation of shape d2y=func(x,y,dy). It returns a matrix yx with the
    ! evolution of y at each iteration.
    ! Scheme:
    !   K1 = f(x_n-1, y_n-1)
    !   K2 = f(x_n, y_n)
    !   y_n+1 = y_n + h/2*(-K1 + 3*K2)

    ! VARIABLES:
    ! (input)
    ! odefunc     --> Subroutine carrying the y-->dy transformations characteristic of the problem at each abscissa value
    ! N           --> Number of iterations of the method
    ! xlow, xhigh --> Interval of x in which the evolution of the coordinates y is to be determined
    ! Ndim        --> Dimensions of the coordinate vector y (nº of coordinates to integrate)
    ! y0          --> Vector with Ndim dimensions carrying the initial conditions of the coordinates y
    ! (output)
    ! yx          --> Matrix of dim(Ndim,N+1) that records the evolution of each coordinate y for each iteration of the method (includes y0)
    ! (internal)
    ! h           --> Step of the interval defined in x
    ! x1, x2      --> Independent variable of the problem (parametric) for points 1 and 2
    ! y1, y2      --> Vector with Ndim dimensions that collects the last computed value of the coordinates y at each iteration for points 1 and 2
    ! dy1, dy2    --> Vector with Ndim dimensions that collects the last computed value of the derivatives dy at each iteration for points 1 and 2
    ! i, j        --> Counters

    ! We determine the equally spaced step with which the method will advance
    h = abs(xhigh - xlow) / dble(N)

    ! We introduce the initial conditions y0 into the final list matrix for each component y
    do j = 1, Ndim
        yx(j,1) = y0(j)
    enddo

    ! We initialize x1 and y1 with the initial conditions
    x1 = xlow
    y1 = y0

    ! We generate dy1 from x1 and y1
    call odefunc(Ndim, x1, y1, dy1)

    ! We execute an initial Euler step and introduce the results into y2
    do j = 1,Ndim
        yx(j,2) = yx(j,1) + h*dy1(j)
        y2(j) = yx(j,2)
    enddo

    ! We execute the multi-step Adams-Bashforth method for the remaining N-1 iterations
    do i = 1, N-1
        ! We generate the previous abscissas x1 and x2
        x1 = xlow + (i-1)*h
        x2 = xlow + i*h

        ! We call subroutine "odefunc()" so that it computes the derivatives dy from the abscissa x and the coordinates y of the previous step
        call odefunc(Ndim, x1, y1, dy1)
        call odefunc(Ndim, x2, y2, dy2)

        ! In a loop that goes through each of the coordinates we compute by Adams-Bashforth its next value with dy1 and dy2 and record it in the yx matrix
        do j = 1, Ndim
            yx(j,i+2) = yx(j,i+1) + h/2.d0*(-dy1(j) + 3*dy2(j))
            ! We update the coordinates y1 and y2
            y1(j) = yx(j,i+1)
            y2(j) = yx(j,i+2)
        enddo

    enddo

    return
end


subroutine rungekutta4order(h, neqs, yin, odefunc, yout)
    implicit none
    integer :: neqs
    double precision :: h, yin(neqs), yout(neqs)
    double precision :: dy1(neqs), y2(neqs), dy2(neqs), y3(neqs), dy3(neqs), y4(neqs), dy4(neqs)
    double precision :: x, x2, x3, x4

    COMMON/ABSCISSA/x

    ! DESCRIPTION:
    ! This subroutine executes the 4th-order Runge-Kutta method to integrate a single step of the coordinate vector y with its derivatives dy from
    ! an ordinary differential equation of the type d2y = func(x,y,dy). It returns a vector yout with the integrated values of y and dy.
    ! Scheme:
    !   K1 = f(x, y_in)
    !   K2 = f(x + h/2, y_in + h/2*K1)
    !   K3 = f(x + h/2, y_in + h/2*K2)
    !   K4 = f(x + h, y_in + h*K3)
    !   y_out = y_in + h/6*(K1 + 2*K2 + 2*K3 + K4) 

    ! VARIABLES:
    ! (input)
    ! h                  --> Predefined step that the method will use
    ! neqs               --> Number of dimensions of the physical problem (number of generalized coordinates y)
    ! yin                --> Vector of neqs components with the coordinates to be integrated (y)
    ! odefunc            --> Subroutine with the specific ODE of the problem that returns the derivatives of y and dy at each abscissa
    ! (output)
    ! yout               --> Vector of neqs components with the already integrated coordinates (y')
    ! (internal)
    ! x, x2, x3, x4      --> Abscissas of the problem (parametric) for the coefficients K1, K2, K3 and K4 respectively
    ! y2, y3, y4         --> Vectors with neqs dimensions that carry the value of the coordinates y for K2, K3 and K4 respectively
    ! dy1, dy2, dy3, dy4 --> Vectors with neqs dimensions that collect the value of the derivatives of y1 (i.e. K1), y2 (i.e. K2), etc. respectively

    ! We call subroutine 'odefunc()' so that it computes the derivatives (dy1) of the coordinate vector (yin) at the abscissa x
    call odefunc(neqs, x, yin, dy1)

    ! We compute the modified abscissa and coordinates for the K2 parameter
    x2 = x + h/2.d0
    y2 = yin + h/2.d0 * dy1

    ! We call 'odefunc()' so that it returns K2 from x2 and y2
    call odefunc(neqs, x2, y2, dy2)

    ! We do the same for K3 and K4
    x3 = x + h/2.d0
    y3 = yin + h/2.d0 * dy2

    call odefunc(neqs, x3, y3, dy3)

    x4 = x + h
    y4 = yin + h*dy3

    call odefunc(neqs, x4, y4, dy4)

    ! We compute the vector of integrated coordinates (yout) with the 4th-order Runge-Kutta scheme
    yout = yin + h/6.d0*(dy1 + 2*dy2 + 2*dy3 + dy4)

    return
end


subroutine shooting(Nstep, xi, xf, fint, odefunc, Ndim, y0, E1, E2, En)
    implicit none
    integer :: i, Nstep, Ndim, count
    double precision :: xi, xf, h, x
    double precision :: E, E1, E2, En
    double precision :: y0(Ndim), y(Ndim), yx1(Ndim), yx2(Ndim), yxn(Ndim), res

    COMMON/ENERGY/E
    COMMON/ABSCISSA/x

    external odefunc

    ! DESCRIPTION:
    ! This subroutine executes the shooting method to find an eigenvalue 'En' of the stationary Schrödinger equation for an infinite-barrier
    ! potential from some initial conditions (y0) and boundary conditions (function(x=1) = 0). In order to refine the result, it uses the
    ! secant method up to a certain error margin.

    ! VARIABLES:
    ! (input)
    ! Nstep                     --> Number of steps in which the integrator subroutine is executed
    ! xi, xf                    --> Endpoints of the abscissa interval in which the ODE of the problem is integrated
    ! odefunc                   --> Subroutine that returns the derivatives of the coordinates y of the problem
    !                               (it is included within the integrator subroutine)
    ! fint(h,y,yx,Ndim,odefunc) --> Integrator subroutine that integrates one step h of the coordinate vector with Ndim dimensions y from
    !                               the functions of the physical problem provided by odefunc and returns the new coordinate vector yx
    ! Ndim                      --> Number of dimensions of the physical problem (i.e. number of coordinates and generalized momenta)
    ! y0                        --> Vector of Ndim dimensions with the initial conditions of the coordinates y of the problem
    ! E1, E2                    --> Energy values from which the secant method starts
    ! (output)
    ! En                        --> Eigenvalue found
    ! (internal)
    ! E, x                      --> Values of the energy and of the abscissa at each moment
    ! y                         --> Vector of Ndim dimensions that carries at each case the coordinates and generalized momenta
    ! yx1, yx2, yxn             --> Vectors of Ndim dimensions that carry the already integrated values of the coordinate vector y
    ! res                       --> Last computed value of the integrated coordinate yxn with the En value computed by the secant method.
    !                               The proximity of this value to the boundary condition of the problem is what determines whether the 
    !                               method continues executing to further refines the result
    ! count                     --> Counter of the number of times the secant method is executed until the desired eigenvalue is found
    ! i                         --> Internal counter

    ! We compute the step that we will use for the integration method of subroutine 'fint()'
    h=(xf-xi)/dble(Nstep)

    ! We initialize the execution counter of the method to 0
    count = 0
    ! We set an arbitrarily high value of the comparator parameter so that we can enter the while loop and start the shooting method
    res = 1.d5

    ! While the value of the function at the last abscissa (res = yxn(x=1)) is not below the indicated margin (1.d-5)
    ! the shooting method will be executed repeatedly
    do while (abs(res).ge.1.d-5)
        ! We assign the energy value E1 to the energy E and the initial coordinates y0 to the vector y
        E = E1
        y = y0

        ! We integrate the coordinates y Nstep times with the integrator subroutine 'fint()'
        do i = 1, Nstep
            ! We define the abscissa in each case
            x = xi + (i-1)*h

            ! We call the integrator subroutine with step h, the coordinates y, the dimensions of the problem Ndim and the equations of
            ! the problem 'odefunc()' so that it returns the integrated coordinates yx1
            call fint(h, Ndim, y, odefunc, yx1)

            ! We reassign to the vector y of the next step the integrated coordinates of the previous step
            y = yx1
        enddo
        ! With this, we keep only the last integrated value of yx1, i.e. that of the function at x = 1

        ! We do the same for the energy E2
        E = E2
        y = y0

        do i = 1, Nstep
            x = xi + (i-1)*h

            call fint(h, Ndim, y, odefunc, yx2)

            y = yx2
        enddo

        ! With the secant method, we compute the next energy value En
        En = (E1*yx2(1) - E2*yx1(1)) / (yx2(1) - yx1(1))

        ! We do the same as with E1 and E2
        E = En
        y = y0

        do i = 1, Nstep
            x = xi + (i-1)*h

            call fint(h, Ndim, y, odefunc, yxn)

            y = yxn
        enddo

        ! We assign the last computed value of yxn (function(x=1)) to the comparator parameter 'res'. If this is below the established error margin,
        ! the method stops and returns the resulting eigenvalue En from the secant method with which yxn has been computed
        res = yxn(1)

        ! We add 1 to the execution counter
        count = count+1

        ! We write a record in the data file of the evolution of the En value at each execution until it becomes the eigenvalue of the problem
        write(1,*) count, En

        ! If the method has not yet converged, we reassign the values E2-->E1 and En-->E2
        E1 = E2
        E2 = En
    enddo

    return
end