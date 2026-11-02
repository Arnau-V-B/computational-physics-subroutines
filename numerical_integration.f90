! This script contains some useful subroutines to numerically integrate functions.


subroutine trapezoidalrule(x1, x2, k, func, res)
    implicit none
    double precision :: x1, x2, func, res
    double precision :: h
    integer :: k, i

    ! DESCRIPTION:
    ! This subroutine computes the numerical integral of a function func() through the chained trapezoidal rule (2-point) between x1 and x2 in
    ! 3^k equally spaced intervals.

    ! VARIABLES:
    ! (input)
    ! x1, x2 --> Endpoints between which the function is to be integrated
    ! k      --> Number of intervals/steps (3^k) to be used
    ! func() --> External function to be integrated
    ! (output)
    ! res    --> Value of the integral of the function
    ! (internal)
    ! h      --> Equally spaced step of the method
    ! i      --> Counter


    ! We compute the step value
    h = abs(x2 - x1) / dble(3**k)

    ! We add the contribution of the endpoints of the interval (half)
    res = (func(x1) + func(x2)) / 2.d0

    ! For the interior points of the interval we sum each contribution
    do i = 1, 3**k - 1
        res = res + func(x1 + i*h)
    enddo

    ! After summing the contributions of the endpoints and the interior points, we multiply the result by the step
    res = h*res

    return
end


subroutine simpson38(x1, x2, k, func, res)
    implicit none
    double precision :: x1, x2, func, res
    double precision :: h
    integer :: k, i

    ! DESCRIPTION:
    ! This subroutine computes the numerical integral of a function func() through the chained Simpson-3/8 rule (4-point) between x1 and x2 in
    ! 3^k equally spaced intervals.

    ! VARIABLES:
    ! (input)
    ! x1, x2 --> Endpoints between which the function is to be integrated
    ! k      --> Number of intervals/steps (3^k) to be used
    ! func() --> External function to be integrated
    ! (output)
    ! res    --> Value of the integral of the function
    ! (internal)
    ! h      --> Equally spaced step of the method
    ! i      --> Counter


    ! We compute the step
    h = abs(x2 - x1) / dble(3**k)

    ! We initialize the res variable to 0
    res = 0.d0

    ! We add to itself the Simpson 3/8 sequence as many times as the interval has steps in groups of three: for i=0,3,6,...,3^k-3
    ! For each "i" the components "i", "i+1", "i+2" and "i+3" are involved; in the next step "i+3" becomes "i" (hence we go in steps of three)
    do i=0, 3**k - 3, 3
        res = res + func(x1 + i*h) + 3*func(x1 + (i+1)*h) + 3*func(x1 + (i+2)*h) + func(x1 + (i+3)*h)
    enddo

    ! We apply the 3/8 factor of the step to the final result to complete the Simpson-3/8 formula
    res = ((3*h) / 8.d0) * res
    return
end


subroutine tabderivative(ndata, valsx, f, df)
    implicit none
    integer :: ndata, i
    double precision :: valsx(ndata), f(ndata), df(ndata)
    double precision :: h

    ! DESCRIPTION:
    ! This subroutine computes the two-point numerical derivative from two lists of ndata values: one with the abscissas and the other with
    ! the values of the function at each abscissa. It returns a list of equal dimension with the values of the derivative at each point.

    ! VARIABLES:
    ! (input)
    ! ndata --> Number of dimensions of the input point lists
    ! valsx --> Vector of dimension ndata carrying the abscissas
    ! f     --> Vector of dimension ndata carrying the values of the function at each abscissa
    ! (output)
    ! df    --> Vector of dimension ndata carrying the values of the derivative of the function at each abscissa
    ! (internal)
    ! h     --> Equidistant step between abscissas of the valsx list
    ! i     --> Counter


    ! We compute the step value (assuming it is homogeneous)
    h = abs(valsx(2) - valsx(1))

    ! We compute the derivative at the initial endpoint (forward)
    df(1) = (f(2) - f(1)) / h

    ! We compute the derivatives at the middle points (centered)
    do i = 2, ndata - 1
        df(i) = (f(i+1) - f(i-1)) / (2*h)
    enddo

    ! We compute the derivative at the final endpoint (backward)
    df(ndata) = (f(ndata) - f(ndata-1)) / h
    return
end
