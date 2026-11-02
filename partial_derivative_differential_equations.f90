! This script contains some useful subroutines to solve differential equations with partial derivatives.


subroutine poisson_2D_dir(icontrol, Nx, Ny, h, T0, Tint, S, tol, Txy)
    implicit none
    integer :: icontrol, Nx, Ny, Px, Py, count
    double precision :: h, T0(4), Tint, S, tol
    double precision :: Txy_old(Nx,Ny), Txy(Nx,Ny), Terr(Nx-2,Ny-2)
    double precision :: x, y, err, w
    integer :: i, j

    ! DESCRIPTION:
    ! This subroutine characterizes the steady-state temperature profile of a two-dimensional rectangular plate from an initial temperature distribution,
    ! with thermal sources and Dirichlet boundary conditions, performing all the calculations inside a mesh/grid of points. To do so, it can use the
    ! Jacobi (0), over-relaxation (1) or Gauss-Seidel (2) methods depending on the value of the parameter 'icontrol' introduced.

    ! VARIABLES:
    ! (input)
    ! icontrol --> 0, 1 or 2 depending on whether the Jacobi, over-relaxation or Gauss-Seidel method is to be executed respectively
    ! Nx, Ny   --> Total number of points in X and Y of the grid with which we want to work
    ! h        --> Step between each point (separation). It controls, therefore, the spatial resolution of the calculation
    ! T0       --> Vector of four dimensions with which the conditions of the four boundaries are introduced. T0 = (T(0,y),T(Lx,y),T(x,0),T(x,Ly))
    ! Tint     --> Initial uniform temperature value at the interior points with which the calculation starts
    ! S(x,y)   --> 2D function that characterizes all the sources on the plate for any point in space
    ! tol      --> Error tolerance value sufficient for the program to stop
    ! (output)
    ! Txy      --> Matrix of dimensions (Nx,Ny) carrying the temperature values at each point of the plate once the steady state is reached
    ! (internal)
    ! Px, Py   --> Allow tracking the temperature evolution at a specific point of the plate as a function of the iterations
    ! count    --> Iteration counter
    ! Txy_old  --> Matrix of dimensions (Nx,Ny) carrying the temperature values at each point of the plate computed in the previous iteration (Jacobi method)
    ! Terr     --> Matrix of dimensions (Nx-2,Ny-2) (interior points only) that measures the error in the temperature calculation at each point for each iteration
    ! x, y     --> Values of the abscissa and ordinate at each moment
    ! err      --> Variable carrying the maximum error value within the Terr matrix. When this is less than the input tolerance, the method stops
    ! w        --> Parameter defining the convergence of the over-relaxation method (is defined externally to the subroutine)
    ! i, j     --> Counters in the X and Y dimensions respectively

    COMMON/PARAMETER/w

    ! We determine the coordinates of the point whose temperature evolution we want to track
    Px = int(7.5d0/h) + 1
    Py = int(23.5d0/h) + 1

    ! We generate the Txy_old and Txy matrices with the defined boundary and initial conditions
    do i = 1, Nx
        do j = 1, Ny
            if (i.eq.1) then        ! Boundary (0,y)
                Txy_old(i,j) = T0(1)
                Txy(i,j) = T0(1)
            else if (i.eq.Nx) then  ! Boundary (Lx,y)
                Txy_old(i,j) = T0(2)
                Txy(i,j) = T0(2)
            else if (j.eq.1) then   ! Boundary (x,0)
                Txy_old(i,j) = T0(3)
                Txy(i,j) = T0(3)
            else if (j.eq.Ny) then  ! Boundary (x,Ly)
                Txy_old(i,j) = T0(4)
                Txy(i,j) = T0(4)
            else                    ! Interior points
                Txy_old(i,j) = Tint
                Txy(i,j) = Tint
            endif 
        enddo
    enddo

    ! We start with an arbitrarily large error value so as to enter the while loop
    err = 1.d6
    ! We start the iteration counter of the method
    count = 0

    ! We print to the text file the initial temperature value of the tracking point
    write(1,*) count, Txy_old(Px,Py)

    ! We start a while loop that won't stop the calculations until the desired error tolerance is reached or until the limit of iterations is exceeded
    do while ((err.gt.tol).and.(count.eq.100000))

        ! We go through all the interior points of the grid
        do i = 2, Nx-1
            do j = 2, Ny-1

                ! We define the coordinates of each point
                x = (i-1)*h
                y = (j-1)*h

                ! We execute the Jacobi method if icontrol=0
                if (icontrol.eq.0) then
                    ! Calculation of the T value
                    Txy(i,j) = (Txy_old(i+1,j) + Txy_old(i-1,j) + Txy_old(i,j+1) + Txy_old(i,j-1) + h**2 * S(x,y)) / 4.d0
                    ! Error in the calculation (temperature difference with respect to the previous step)
                    Terr(i-1,j-1) = abs(Txy(i,j) - Txy_old(i,j))

                ! We execute the over-relaxation method if icontrol=1
                else if (icontrol.eq.1) then
                    ! Error in the calculation
                    ! (it is at the beginning since it depends on the current T value and with over-relaxation the next one is computed)
                    Terr(i-1,j-1) = abs(w*(Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y) - 4*Txy(i,j))) / 4.d0
                    ! Calculation of the next T value
                    Txy(i,j) = Txy(i,j) + w*(Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y) - 4*Txy(i,j)) / 4.d0

                ! We execute the Gauss-Seidel method if icontrol=2
                else if (icontrol.eq.2) then
                    ! Error in the calculation (same as over-relaxation but without the w parameter)
                    Terr(i-1,j-1) = abs(Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y) - 4*Txy(i,j)) / 4.d0
                    ! Calculation of the next T value
                    Txy(i,j) = (Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y)) / 4.d0

                endif

            enddo
        enddo

        ! At each iteration, we write the temperature value at the tracking point to the text file
        write(1,*) count, Txy(Px,Py)

        ! We take the maximum error over the whole grid
        err = maxval(Terr)
        ! We add 1 to the iteration counter
        count = count+1
        ! We update the Txy_old matrix for the next step
        Txy_old = Txy

    enddo

    ! We write to screen the total number of iterations the program needed to reach the steady state
    write(1,*) 
    write(1,*) 
    write(*,*) "Number of iterations:", count

    return
end


subroutine poisson_2D_neu(icontrol, Nx, Ny, h, T0, Tint, S, tol, Txy)
    implicit none
    integer :: icontrol, Nx, Ny, Px, Py, count
    double precision :: h, T0(4), Tint, S, tol
    double precision :: Txy_old(Nx+2,Ny), Txy(Nx+2,Ny), Terr(Nx,Ny-2)   ! Add +1 to the number of points according to the number of Neumann BCs
    double precision :: x, y, err, w
    integer :: i, j

    COMMON/PARAMETER/w

    ! DESCRIPTION:
    ! This subroutine characterizes the steady-state temperature profile of a two-dimensional rectangular plate from an initial temperature distribution,
    ! with thermal sources and either Dirichlet or Neumann boundary conditions, performing all the calculations inside a mesh/grid of points. To do so, 
    ! it can use the Jacobi (0), over-relaxation (1) or Gauss-Seidel (2) methods depending on the value of the parameter 'icontrol' introduced.

    ! VARIABLES:
    ! (input)
    ! icontrol --> 0, 1 or 2 depending on whether the Jacobi, over-relaxation or Gauss-Seidel method is to be executed respectively
    ! Nx, Ny   --> Total number of points in X and Y of the grid with which we want to work
    ! h        --> Step between each point (separation). It controls, therefore, the spatial resolution of the calculation
    ! T0       --> Vector of four dimensions with which the conditions of the four boundaries are introduced. T0 = (T(0,y),T(Lx,y),T(x,0),T(x,Ly))
    ! Tint     --> Initial uniform temperature value at the interior points with which the calculation starts
    ! S(x,y)   --> 2D function that characterizes all the sources on the plate for any point in space
    ! tol      --> Error tolerance value sufficient for the program to stop
    ! (output)
    ! Txy      --> Matrix of dimensions (Nx,Ny) carrying the temperature values at each point of the plate once the steady state is reached
    ! (internal)
    ! Px, Py   --> Allow tracking the temperature evolution at a specific point of the plate as a function of the iterations
    ! count    --> Iteration counter
    ! Txy_old  --> Matrix of dimensions (Nx,Ny) carrying the temperature values at each point of the plate computed in the previous iteration (Jacobi method)
    ! Terr     --> Matrix of dimensions (Nx-2,Ny-2) (interior points only) that measures the error in the temperature calculation at each point for each iteration
    ! x, y     --> Values of the abscissa and ordinate at each moment
    ! err      --> Variable carrying the maximum error value within the Terr matrix. When this is less than the input tolerance, the method stops
    ! w        --> Parameter defining the convergence of the over-relaxation method (is defined externally to the subroutine)
    ! i, j     --> Counters in the X and Y dimensions respectively

    ! Add +1 if the Neumann condition is at x=0 or y=0 (axis shift)
    Px = int(8.d0/h) + 2
    Py = int(13.5d0/h) + 1

    ! Add +1 to the loop as many times as there are Neumann BCs (maximum case: +2 on each axis with respect to the original loop)

    ! Interior points
    do i = 2, Nx+1
        do j = 2, Ny-1
            Txy_old(i,j) = Tint
            Txy(i,j) = Tint
        enddo
    enddo

    ! Boundaries (0,y) and (Lx,y)
    do j=1,Ny
        ! Dirichlet BCs
        ! Txy_old(1,j) = T0(1)
        ! Txy(1,j) = T0(1)
        ! Txy_old(Nx,j) = T0(2)
        ! Txy(Nx,j) = T0(2)
        ! Neumann BCs
        Txy_old(1,j) = Txy_old(3,j) + 2*h * T0(1)
        Txy(1,j) = Txy(3,j) + 2*h * T0(1)
        Txy_old(Nx+2,j) = Txy_old(Nx,j) + 2*h * T0(2)
        Txy(Nx+2,j) = Txy(Nx,j) + 2*h * T0(2)
    enddo

    ! Boundaries (x,0) and (x,Ly)
    do i=1,Nx+2
        ! Dirichlet BCs
        Txy_old(i,1) = T0(3)
        Txy(i,1) = T0(3)
        Txy_old(i,Ny) = T0(4)
        Txy(i,Ny) = T0(4)
        ! Neumann BCs
        ! Txy_old(i,1) = Txy_old(i,3) + 2*h * T0(3)
        ! Txy(i,1) = Txy(i,3) + 2*h * T0(3)
        ! Txy_old(i,Ny+2) = Txy_old(i,Ny) + 2*h * T0(4)
        ! Txy(i,Ny+2) = Txy(i,Ny) + 2*h * T0(4)
    enddo    

    ! We start with an arbitrarily large error value so as to enter the while loop
    err = 1.d6
    ! We start the iteration counter of the method
    count = 0

    ! We print to the text file the initial temperature value of the tracking point
    write(1,*) count, Txy_old(Px,Py)

    ! We start a while loop that won't stop the calculations until the desired error tolerance is reached or until the limit of iterations is exceeded
    do while ((err.gt.tol).or.(count.eq.100000))
        ! Add +1 to Nx-1/Ny-1 depending on the number of Neumann BCs
        do i = 2, Nx+1
            do j = 2, Ny-1

                ! Subtract -1 if the Neumann BC is at x=0 or y=0 (axis shift)
                x = (i-2)*h
                y = (j-1)*h

                ! We execute the Jacobi method if icontrol=0
                if (icontrol.eq.0) then
                    ! Calculation of the T value
                    Txy(i,j) = (Txy_old(i+1,j) + Txy_old(i-1,j) + Txy_old(i,j+1) + Txy_old(i,j-1) + h**2 * S(x,y)) / 4.d0
                    ! Error in the calculation (temperature difference with respect to the previous step)
                    Terr(i-1,j-1) = abs(Txy(i,j) - Txy_old(i,j))

                ! We execute the over-relaxation method if icontrol=1
                else if (icontrol.eq.1) then
                    ! Error in the calculation
                    ! (it is at the beginning since it depends on the current T value and with over-relaxation the next one is computed)
                    Terr(i-1,j-1) = abs(w*(Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y) - 4*Txy(i,j))) / 4.d0
                    ! Calculation of the next T value
                    Txy(i,j) = Txy(i,j) + w*(Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y) - 4*Txy(i,j)) / 4.d0

                ! We execute the Gauss-Seidel method if icontrol=2
                else if (icontrol.eq.2) then
                    ! Error in the calculation (same as over-relaxation but without the w parameter)
                    Terr(i-1,j-1) = abs(Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y) - 4*Txy(i,j)) / 4.d0
                    ! Calculation of the next T value
                    Txy(i,j) = (Txy(i+1,j) + Txy(i-1,j) + Txy(i,j+1) + Txy(i,j-1) + h**2 * S(x,y)) / 4.d0

                endif

                ! We apply the Neumann BCs to the y endpoints
                Txy(1,j) = Txy(3,j) + 2*h * T0(1)
                Txy(Nx+2,j) = Txy(Nx,j) + 2*h * T0(2)    
            enddo
            ! We apply the Neumann BCs to the x endpoints
            !Txy(i,1) = Txy(i,3) + 2*h * T0(3)
            !Txy(i,Ny+2) = Txy(i,Ny) + 2*h * T0(4)
        enddo

        ! At each iteration, we write the temperature value at the tracking point to the text file
        write(1,*) count, Txy(Px,Py)

        ! We take the maximum error over the whole grid
        err = maxval(Terr)
        ! We add 1 to the iteration counter
        count = count + 1
        ! We update the Txy_old matrix for the next step
        Txy_old = Txy

    enddo

    ! We write to screen the total number of iterations the program needed to reach the steady state
    write(1,*) 
    write(1,*) 
    write(*,*) "Number of iterations:", count

    return
end


subroutine poisson_1D_loss(Nx, h, Tc, Tint, tol, Tx_new)
    implicit none
    integer :: Nx
    double precision :: h, Tc(2), Tint, tol
    double precision :: Tx_new(0:Nx), Tx_old(0:Nx), Terr(1:Nx-1)
    double precision :: D, err
    integer :: i

    COMMON/PARAMETER/D

    ! DESCRIPTION:
    ! This subroutine characterizes the steady-state temperature profile of a 1D metal bar with losses from an initial temperature distribution
    ! assuming Dirichlet BCs. To do so, it solves the 1D Poisson equation in partial derivatives using a method similar to Jacobi's.

    ! VARIABLES:
    ! (input)
    ! Nx     --> Total number of grid points in X with which we want to work
    ! h      --> Step between each point (separation). It controls, therefore, the spatial resolution of the calculation
    ! Tc     --> Vector of two dimensions with which the conditions of the two boundaries of the bar are introduced Tc = (T(0),T(L))
    ! Tint   --> Initial uniform temperature value at the interior points of the bar with which the calculation starts
    ! tol    --> Error tolerance value sufficient for the program to stop
    ! (output)
    ! Tx_new --> Vector of dimension (0:Nx) carrying the temperature values at each point of the bar once the steady state is reached
    ! (internal)
    ! Tx_old --> Vector of dimension (0:Nx) carrying the temperature values at each point of the bar computed in the previous iteration
    ! Terr   --> Vector of dimension (1,Nx-1) (interior points only) that measures the error in the temperature calculation at each point for each iteration
    ! err    --> Variable carrying the maximum error value within the Terr vector. When this is less than the input tolerance, the method stops
    ! D      --> Parameter that weights the importance of the thermal losses of the bar
    ! i      --> Counter in the X dimension


    ! Initial and boundary conditions
    Tx_old = Tint
    Tx_old(0) = Tc(1)
    Tx_new(0) = Tc(1)
    Tx_old(Nx) = Tc(2)
    Tx_new(Nx) = Tc(2)

    ! We start the error with an arbitrarily large value so that it can enter the while loop
    err = 1.d6

    ! While the desired error tolerance is not reached, the method will execute indefinitely
    do while (err.gt.tol)

        ! We compute the temperature profile at each iteration, together with the error associated with the calculation (only for the interior points)
        do i = 1, Nx-1
            Tx_new(i) = (Tx_old(i+1) + Tx_old(i-1)) / (D*h**2 + 2.d0)
            Terr(i) = abs(Tx_new(i) - Tx_old(i))
        enddo

        ! We look for the maximum error value
        err = maxval(Terr)

        ! We update the T_old vector with the T_new found to start the next iteration
        Tx_old = Tx_new
    enddo

    return
end


subroutine crank_nicolson(icontrol, Nx, dx, Nt, dt, alpha, T_new)
    implicit none
    integer :: icontrol, Nx, Nt
    double precision :: x, dx, t, dt, alpha, beta, r
    double precision :: T_old(0:Nx), T_new(0:Nx), T_right(1:Nx-1), T_tri(1:Nx-1)
    double precision :: a(1:Nx-1), b(1:Nx-1), c(1:Nx-1), B_right(1:Nx-1,1:Nx-1)
    double precision :: x1, x2, x3, sum, Tmean, Tamb, L
    integer :: i, j, k

    COMMON/DATA/beta,Tamb
    COMMON/DIMENSIONS/L


    ! DESCRIPTION:
    ! This subroutine executes the Crank-Nicolson method for a 1D metal bar with losses from an initial temperature distribution, assuming Dirichlet BCS.
    ! It is an implicit method of fast resolution since it converts the equation dT/dt = alpha(x) d2T/dx2 - beta(T - Tamb) into a system of tridiagonal
    ! matrices on the left and right (A·T_n+1 = B·T_n). Thus it needs an external subroutine tridiag() that efficiently solves this system.

    ! VARIABLES:
    ! (input)
    ! icontrol      --> Depending on the value introduced (0, 1, 2 or 3) the subroutine prints a series of information to screen or to a text file
    ! Nx, Nt        --> Total number of points in x (to characterize the bar) and in t (to characterize time)
    ! dx, dt        --> Steps between each point of x and t respectively
    ! alpha         --> Function dependent on x that determines the value of the thermal diffusivity of the bar material at each point
    ! T_new         --> Vector of dimension (0:Nx) with the initial temperature distribution of the bar
    ! (output)
    ! T_new         --> Vector of dimension (0:Nx) carrying the temperature values at each point of the bar at each time step. Once the subroutine
    !                   finishes executing, this vector that initially carried the initial conditions of the bar will now carry the final conditions.
    ! (internal)
    ! x, t          --> Variables carrying the spatial and temporal dimensions respectively
    ! T_old         --> Vector of dimension (0:Nx) carrying the temperature values at each point of the bar in the initial conditions
    !                   (this vector remains unchanged)
    ! T_right       --> Vector of dimension (1:Nx-1) resulting from the scalar product (B_right·T_new) on the RHS for the interior points of the bar
    ! T_tri         --> Vector of dimension (1:Nx-1) carrying the temperature values at each point of the bar after the tridiagonal calculation
    ! a, b, c       --> Vectors of dimension (1,Nx-1) that symbolize the three diagonals (left, central and right) of the tridiagonal matrix on the LHS (i.e. A)
    ! B_right       --> Tridiagonal matrix of dimensions (1:Nx-1,1:Nx-1) on the RHS (i.e. B)
    ! beta, r       --> Physical parameters of the calculations (loss coefficient and Courant number, respectively)
    ! x1, x2, x3    --> Position values of the specific points of the bar to be studied
    ! Tamb          --> Ambient temperature value of the medium to correct the computed temperature values
    ! sum, Tmean, L --> Parameters for the statistical study of the average temperature of the bar
    ! i, j, k       --> Counters


    ! We build the operator matrices A and B
    B_right = 0.d0

    ! We iterate through the interior points
    do i = 1, Nx-1
        x = i*dx
        r = alpha(x) * dt / dx**2

        a(i) = -1.d0*r
        b(i) = 2.d0*(1.d0 + r) + beta*dt
        c(i) = -1.d0*r

        B_right(i,i) = 2.d0*(1.d0 - r) - beta*dt
        if (i.lt.Nx-1) then
            B_right(i,i+1) = r
            B_right(i+1,i) = r
        endif
    enddo

    ! We eliminate the contributions at the positions next to the boundaries
    a(1) = 0.d0
    c(Nx-1) = 0.d0

    ! We save the initial conditions to not lose them in the calculation
    T_old = T_new

    ! We start the time loop
    do k = 1, nt
        t = k*dt

        ! We compute the T_right vector on the RHS as the scalar product between matrix B_right and vector T_new (i.e B_right·T_new)
        ! This product will have two fewer dimensions than the original vectors T_old and T_new, since the values at the boundaries are fixed
        do i = 1, Nx-1
            T_right(i) = 0.d0
            do j = 1, Nx-1
                T_right(i) = T_right(i) + B_right(i,j) * T_new(j)
            enddo
        enddo

        ! We add the boundary conditions at the points contiguous to the fixed boundaries
        T_right(1) = T_right(1) + 2.d0*(alpha(0.d0) * dt/dx**2) * T_old(0)
        T_right(Nx-1) = T_right(Nx-1) + 2.d0*(alpha((Nx-1) * dx) * dt/dx**2) * T_old(Nx)

        ! To work with the tridiagonal matrix we only need the interior points of the T_new vector
        ! We work, therefore, with an auxiliary unknown vector T_tri with two fewer dimensions
        do i = 1, Nx-1
            T_tri(i) = T_new(i)
        enddo

        ! We call subroutine tridiag() so that it computes the components of the new temperature vector T_tri
        call tridiag(a, b, c, T_right, T_tri, Nx-1)

        ! We undo the change T_new --> T_tri that we made before tridiagonalizing
        do i=1,Nx-1
            T_new(i) = T_tri(i)
        enddo
        ! This obtained T_new vector already contains the updated temperature distribution at time t=k*dt


        ! Depending on the value introduced for the icontrol variable, the subroutine will write certain information

        ! Write nothing if icontrol=0
        if (icontrol.eq.0) then
            ! We write a message to screen informing that the subroutine has been executed
            if (k.eq.nt) write(*,*) "Crank-Nicolson completed!"

        ! Keep track of the evolution of temperature at specific points of the bar if icontrol=1
        else if (icontrol.eq.1) then
            ! Point positions
            x1 = 0.06d0
            x2 = 0.42d0
            x3 = 1.26d0

            ! We include the initial state
            if (k.eq.1) then
                write(1,*) 0.d0, T_old(int(x1/dx)) + Tamb, T_old(int(x2/dx)) + Tamb, T_old(int(x3/dx)) + Tamb
            endif

            ! For each t we write the temperature value at the different points
            write(1,*) t, T_new(int(x1/dx)) + Tamb, T_new(int(x2/dx)) + Tamb, T_new(int(x3/dx)) + Tamb

        ! Keep track of the evolution of the average temperature of the bar if icontrol=2
        else if (icontrol.eq.2) then
            ! In this case, we use the trapezoidal rule to numerically compute the integral T(x,t)dx

            ! We include the initial state
            if (k.eq.1) then
                sum = (T_old(0) + T_old(Nx)) / 2.d0
                do i=1,Nx-1
                    sum = sum + T_old(i)
                enddo
                sum = dx * sum

                Tmean = (1/L) * sum

                write(1,*) 0.d0, Tmean + Tamb
            endif

            ! For each t, we compute the average temperature Tmean
            sum = (T_new(0) + T_new(Nx)) / 2.d0
            do i = 1, Nx-1
                sum = sum + T_new(i)
            enddo
            sum = dx * sum

            Tmean = (1/L) * sum

            write(1,*) t, Tmean + Tamb

        ! Keep track of the evolution of the temperature profile of the bar if icontrol=3
        else if (icontrol.eq.3) then
            ! We include the initial state
            if (k.eq.1) then
                do i = 0, Nx
                    x = i*dx
                    write(1,*) t, x, T_old(i) + Tamb
                enddo
                write(1,*) 
                write(1,*) 
            endif

            ! Every 50 time steps, we write the temperature profile of the bar
            if (mod(k,50).eq.0) then
                do i = 0, Nx
                    x = i*dx
                    write(1,*) t, x, T_new(i) + Tamb
                enddo
                write(1,*) 
                write(1,*) 
            endif
        endif

    enddo

    write(1,*) 
    write(1,*) 

    return
end


subroutine tridiag(A, B, C, R, PSI, IMAX)
    implicit double precision :: (A-H,K,O-Z)
    implicit integer :: (I-J , L-N)
    double precision :: BET
    double precision ::  GAM(4001)
    double precision :: A(IMAX),B(IMAX),C(IMAX),R(IMAX),PSI(IMAX)

    ! This is the subroutine that solves the problem tridiagonal system T·PSI = R;
    ! where T is a tridiagonal matrix and A (lower), B (central), C (upper) its subdiagonals
    !  A   0 a1, a2, a3, ...., aIMAX
    !  B   b1 b2, b3, b4, ...., bIMAX
    !  C   c1, c2, c3, c4, ...., 0 

    if(B(1).EQ.0.) pause
    BET = B(1)
    PSI(1) = R(1) / BET

    do J = 2, IMAX
        GAM(J) = C(J-1) / BET
        BET = B(J) - A(J) * GAM(J)
        if(BET.EQ.0) pause
        PSI(J) = (R(J) - A(J) * PSI(J-1)) / BET
    enddo

    do J = IMAX-1, 1, -1
        PSI(J) = PSI(J) - GAM(J+1) * PSI(J+1)
    enddo

   return
end
