! This script contains some useful subroutines to perform interpolation or finding the roots (i.e. f(x_root) = 0) of functions.


! INTERPOLATION _____________________________________________________________________________________________________________________________________

subroutine interpol_0(ndim, x, y, xin, yout)
    implicit none
    integer :: ndim, i
    double precision :: x(ndim), y(ndim), xin, yout

    ! DESCRIPTION:
    ! This subroutine computes the zeroth-order interpolation of a sequence of ndim points (x,y) at the abscissa xin and returns the corresponding
    ! value yout. Understanding the zeroth-order interpolation of a function f(x) as that the value of the function in an interval x_i - x_(i+1)
    ! remains constant throughout the interval with the initial value f(x_i).

    ! VARIABLES:
    ! (input)
    ! ndim --> Number of known points
    ! x, y --> Vectors of dimension ndim containing the abscissas and ordinates of the known points respectively
    ! xin  --> Value of the abscissa at which the interpolated value is to be found
    ! (output)
    ! yout --> Interpolated ordinate of the introduced abscissa
    ! (internal)
    ! i    --> Counter


    ! We loop over all the grid points. We subtract 1 from the dimension to avoid an overshoot
    do i=1,ndim-1
        ! We look for the relative position of the introduced abscissa within the known grid of points and perform the interpolation
        if ((xin.ge.x(i)).and.(xin.le.x(i+1))) then
            yout=y(i)
        endif

    enddo

    return
end


subroutine interpol_1(ndim, x, y, xin, yout)
    implicit none
    integer :: ndim, i
    double precision :: x(ndim), y(ndim), xin, yout

    ! DESCRIPTION:
    ! This subroutine computes the linear interpolation (via Lagrange polynomials) of a sequence of ndim points (x,y) at the abscissa xin and returns
    ! the corresponding value yout.

    ! VARIABLES:
    ! (input)
    ! ndim --> Number of known points
    ! x, y --> Vectors of dimension ndim containing the abscissas and ordinates of the known points respectively
    ! xin  --> Value of the abscissa at which the interpolated value is to be found
    ! (output)
    ! yout --> Interpolated ordinate of the introduced abscissa
    ! (internal)
    ! i    --> Counter

    
    ! We loop over all the grid points. We subtract 1 from the dimension to avoid an overshoot
    do i=1,ndim-1

        ! We look for the relative position of the introduced abscissa within the known grid of points and perform the interpolation
        if ((xin.ge.x(i)).and.(xin.le.x(i+1))) then
            yout=y(i)+((y(i+1)-y(i))/(x(i+1)-x(i)))*(xin-x(i))
        endif

    enddo

    return
end


! ROOT FINDING ______________________________________________________________________________________________________________________________________

subroutine func(x, valf, valdf)
    implicit none
    double precision :: x, valf, valdf
    double precision :: pi

    ! DESCRIPTION: Common function for the root-finding subroutines

    pi = 2*asin(1.d0)

    valf = sinh(x) * ((-57.d0*pi)/160.d0 + (57.d0/80.d0 + (17.d0*pi)/20.d0)*x - (17.d0/10.d0 + pi/2.d0)*x**2 + x**3)

    valdf = cosh(x) * ((-57.d0*pi)/160.d0 + (57.d0/80.d0 + (17.d0*pi)/20.d0)*x - (17.d0/10.d0 + pi/2.d0)*x**2 + x**3) + &
            sinh(x) * (57.d0/80.d0 + (17.d0*pi)/20.d0 - 2*(17.d0/10.d0 + pi/2.d0)*x + 3*x**2)

    return
end


subroutine bisection(func, A, B, prec, niter, valroot)
    implicit none
    double precision :: A, B, C
    double precision :: prec, valroot
    double precision :: fA, fB, fC, df
    integer :: niter

    ! DESCRIPTION:
    ! This subroutine executes a bisection method to find the solution of a function 'func' between two endpoints A and B. The iterations stop either
    ! when the solution is found directly or when the imposed precision is reached.

    ! VARIABLES:
    ! (input)
    ! func(x,f,df) --> external subroutine that returns the value of a function and its derivative at the introduced abscissa
    ! A, B         --> Values of the endpoints within which the solution (root) is to be found
    ! prec         --> Desired precision
    ! (output)
    ! niter        --> Number of iterations performed until the solution is found
    ! valroot      --> Abscissa where the method considers the solution/root to be located
    ! (internal)
    ! C            --> Midpoint of the interval that is recalculated at each iteration
    ! fA, fB, fC   --> Values of the function returned by subroutine func() when the abscissas A, B and C are introduced respectively
    ! df           --> Value of the derivative of the function at the introduced abscissa (irrelevant for this method)


    ! We call the function values at the endpoints of the interval
    call func(A, fA, df)
    call func(B, fB, df)

    ! We start the iteration counter at 0
    niter = 0

    ! We first check whether either of the introduced values is in fact already one of the roots
    if (fA.eq.0.d0) then
        valroot = A
        return
    else if (fB.eq.0.d0) then
        valroot = B
        return
    endif

    ! Then we check whether in the initially introduced interval there is no sign change and, therefore, no root
    if (fA*fB.gt.0.d0) then
        write(*,*) "There is no root of f(x) in this interval"
        stop
    endif

    ! Finally, while the interval is outside the established precision, the approximation will continue to be iterated
    do while (abs(A-B).gt.prec)
        ! We split the interval in two
        C = (A + B) / 2.d0
        ! We call the function value at the midpoint
        call func(C, fC, df)

        ! We add one iteration to the counter
        niter = niter + 1

        ! We exit the loop if the root is found directly by chance
        if (fC.eq.0.d0) then
            valroot = C
            exit
        endif

        ! We look for which of the two subintervals contains the sign change and, therefore, the root
        if (fA*fC.lt.0.d0) then
            B = C
            fB = fC
        else
            A = C
            fA = fC
        endif

        ! We assign the value of the last computed root at each iteration
        valroot = C
    enddo

    return
end


subroutine regula_falsi(func, A, B, prec, niter, valroot)
    implicit none
    double precision :: A, B, C
    double precision :: prec, err, valroot
    double precision :: fA, fB, fC, df
    integer :: niter

    ! DESCRIPTION:
    ! This subroutine executes a Regula-Falsi method to find the solution of a function 'func' between two endpoints A and B. The iterations stop
    ! either when the solution is found directly or when the imposed precision is reached.

    ! VARIABLES:
    ! (input)
    ! func(x,f,df) --> external subroutine that returns the value of a function and its derivative at the introduced abscissa
    ! A, B         --> Values of the endpoints within which the solution (root) is to be found
    ! prec         --> Desired precision
    ! (output)
    ! niter        --> Number of iterations performed until the solution is found
    ! valroot      --> Abscissa where the method considers the solution/root to be located
    ! (internal)
    ! C            --> Midpoint of the interval that is recalculated at each iteration
    ! err          --> Variable carrying the precision (error in the abscissa) at each iteration
    ! fA, fB, fC   --> Values of the function returned by subroutine func() when the abscissas A, B and C are introduced respectively
    ! df           --> Value of the derivative of the function at the introduced abscissa (irrelevant for this method)


    ! We call the function values at the endpoints of the interval
    call func(A, fA, df)
    call func(B, fB, df)

    ! We start the iteration counter at 0
    niter = 0

    ! We first check whether either of the introduced values is in fact already one of the roots
    if (fA.eq.0.d0) then
        valroot = A
        return
    else if (fB.eq.0.d0) then
        valroot = B
        return
    endif

    ! Then we check whether in the initially introduced interval there is no sign change and, therefore, no root
    if (fA*fB.gt.0.d0) then
        write(*,*) "There is no root of f(x) in this interval"
        stop
    endif

    ! We initialize an arbitrarily large error so we can enter the loop for the first time
    err = 1.d6

    ! Finally, we check whether there really is a sign change.
    ! While this is outside the established precision, the approximation will continue to be iterated
    do while (err.gt.prec)

        ! We compute the midpoint C
        C = (A*fB - B*fA) / (fB - fA)

        ! We compute the error in the abscissa
        err = min(C-A, B-C)

        ! We call the function value at C
        call func(C, fC, df)

        ! We add one iteration to the counter
        niter = niter + 1

        ! We exit the loop if the root is found directly by chance
        if (fC.eq.0.d0) then
            valroot = C
            exit
        endif

        ! We look for which of the two subintervals contains the sign change and, therefore, the root
        if (fA*fC.lt.0.d0) then
            B = C
            fB = fC
        else
            A = C
            fA = fC
        endif

        ! We assign the value of the last computed root at each iteration
        valroot = C
    enddo

    return
end


subroutine newtonraphson(func, xini, prec, niter, valroot)
    implicit none
    double precision :: A, B
    double precision :: xini, prec, valroot
    double precision :: fA, dfA
    integer :: niter

    ! DESCRIPTION:
    ! This subroutine executes a Newton-Raphson method to find the solution of a function 'func' between two endpoints A and B. The iterations stop
    ! either when the solution is found directly or when the imposed precision is reached.

    ! VARIABLES:
    ! (input)
    ! func(x,f,df) --> external subroutine that returns the value of a function and its derivative at the introduced abscissa
    ! xini         --> Abscissa of the initial seed from which the method starts
    ! prec         --> Desired precision
    ! (output)
    ! niter        --> Number of iterations performed until the solution is found
    ! valroot      --> Abscissa where the method considers the solution/root to be located
    ! (internal)
    ! A, B         --> Seeds recalculated at each iteration
    ! fA, dfA      --> Value of the function and its derivative returned by subroutine func() when the abscissa A is introduced

    
    ! We assign the initial seed point and compute the value of the function there
    A = xini

    call func(A, fA, dfA)

    ! We start the iteration counter at 0
    niter = 0

    ! We check whether the introduced seed value is directly the root
    if (fA.eq.0.d0) then
        valroot = A
        write(*,*) "You have directly introduced a zero! Root at x = ", valroot
        return
    endif

    ! If nothing happens, the iteration count begins
    niter = 1

    ! We look for the value of the next seed point
    B = A - (fA/dfA)
    ! While at it, we write the iteration and the provisional root value to a data file
    write(1,*) niter, B

    ! While the imposed precision is not reached, the process will keep iterating indefinitely (up to a limit of 1000 iterations)
    do while ((abs(A-B).gt.prec).and.(niter.le.1000))
        ! We reassign the value of the seed
        A = B
        ! We repeat the previous process again
        call func(A, fA, dfA)

        B = A - (fA/dfA)

        niter=niter+1

        write(1,*) niter, B
    enddo

    ! Finally, once all possible procedures are finished, we assign the last value of the root
    valroot = B
    return
end