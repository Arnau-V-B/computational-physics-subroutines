! This script contains some useful subroutines to generate random number distributions and perform Monte Carlo integrations.

! Reminder: Always start the program with a seed to keep reproducibility:
!   ISEED = .....
!   CALL SRAND(ISEED)


! RANDOM NUMBERS ____________________________________________________________________________________________________________________________________

subroutine histogram(ndata, xdata, xa, xb, nbox, xhis, vhis, errhis, boxsize)
    implicit none
    integer :: ndata, nbox
    double precision :: xdata(ndata), xhis(nbox), vhis(nbox), errhis(nbox)
    double precision :: xa, xb, boxsize
    integer :: i, j, count, sum

    ! DESCRIPTION:
    ! This subroutine generates a normalized histogram of a set of values introduced in a vector within a given interval. It returns the central
    ! positions of each box, the corresponding heights (normalized with the total number of values and the box width), the errors associated
    ! with these heights and the box width.

    ! VARIABLES:
    ! (input)
    ! ndata   --> Number of values to analyze
    ! xdata   --> Vector of dimension ndata carrying the values
    ! xa,xb   --> Interval of abscissas in which the histogram is to be generated
    ! nbox    --> Number of boxes desired for the histogram
    ! (output)
    ! xhis    --> Vector of dimension nbox with the central positions of the boxes
    ! vhis    --> Vector of dimension nbox with the normalized heights of the boxes
    ! errhis  --> Vector of dimension nbox with the error in the calculation of each of the heights
    ! boxsize --> Uniform width of the boxes
    ! (internal)
    ! count   --> Keeps the count of values in each box (reset when changing box)
    ! sum     --> Keeps the total count of values accepted by the histogram (not reset throughout the process)
    ! i, j    --> Counters


    ! We define the width of the histogram boxes
    boxsize = abs(xb - xa) / dble(nbox)

    ! We initialize the sum variable that will keep a count of the total number of data accepted into the histogram
    sum = 0

    ! We start a loop that counts for each box of the histogram
    do i = 1, nbox

        ! We initialize the counter of variables accepted in box i
        count = 0

        ! We compute the central position of the box
        xhis(i) = ((xa + (i-1)*boxsize) + (xa + i*boxsize)) / 2.d0

        ! We start a second loop that goes through the entire introduced list of variables
        do j = 1, ndata
            ! If variable j is within the interval corresponding to box i where we are counting, we add 1 to the counter
            if ((xdata(j).ge.(xhis(i)-boxsize/2.d0)).and.(xdata(j).le.(xhis(i)+boxsize/2.d0))) then
                count = count + 1
            endif
        enddo

        ! As we can see, the sum variable keeps the total count of variables accepted in all boxes analyzed so far
        sum = sum + count

        ! We compute the height of box i as the total number of variables accepted in it (count) and normalize it by dividing by the
        ! total number of data introduced initially and the width of the box intervals
        vhis(i) = count / (ndata*boxsize)

        ! We compute the error in the calculation of this box height
        errhis(i) = sqrt(abs((count / dble(ndata)) * (1.d0 - (count / dble(ndata))))) / (boxsize * sqrt(dble(ndata)))
    enddo

    ! We compare the total number of variables introduced into the histogram with those accepted by it
    ! The two values need not be equal since in certain cases the introduced interval will not coincide with that of the variable list itself
    write(*,*) "Number of data =", ndata, "Accepted =", sum

    return
end


subroutine accept(ndata, func, xlow, xhigh, sup, num)
    implicit none
    integer :: ndata
    double precision :: num(ndata), xlow, xhigh, sup, func
    double precision :: x1, x2, x, y
    integer :: count, i
    double precision :: ev, ev2, var, std

    ! DESCRIPTION:
    ! This subroutine executes an acceptance-rejection method in order to create a set of randomly generated values distributed according to a
    ! given probability density in the form of a mathematical function. Thus, it continuously generates random values in intervals bounded
    ! in xlow and xhigh and only accepts those that fall under the curve (i.e. y < func(x)).

    ! VARIABLES:
    ! (input)
    ! ndata             --> Amount of numbers distributed according to the probability density to be generated
    ! func              --> Mathematical function characterizing the probability density that will govern the accepted numbers
    ! xlow, xhigh       --> Interval of abscissas in which the numbers are to be generated (where the probability function is defined)
    ! sup               --> Maximum value that the random numbers corresponding to the variable y will reach (i.e. maximum of the probability function)
    ! (output)
    ! num               --> Vector of dimension ndata carrying the generated numbers
    ! (internal)
    ! x1, x2            --> Base random numbers (between 0 and 1)
    ! x, y              --> x and y coordinates respectively of the points (appropriately scaled from x1 and x2) that will fit the probability function
    ! count             --> Keeps the total count of numbers accepted by the method
    ! i                 --> Counter
    ! ev, ev2, var, std --> Variables for the statistical study of the resulting number distribution: expected value and its square, variance and
    !                       standard deviation respectively


    ! We initialize the counter that will limit the number of random variables generated by the subroutine
    count = 0

    ! We build a loop that executes the acceptance process repeatedly until ndata random variables have been generated
    do while (count.le.ndata)
        ! We generate two random values uniformly distributed between 0 and 1
        x1 = rand()
        x2 = rand()

        ! We scale the x coordinate from x1 so that it is uniformly distributed within the desired interval
        x = abs(xhigh - xlow)*x1 + xlow
        ! We scale the y coordinate from x2 so that it is uniformly distributed between 0 and the introduced maximum bound
        y = sup * x2

        ! When the randomly generated point (x,y) lies under the curve of the introduced function, we accept it
        if (y.le.func(x)) then
            ! We add 1 to the counter of accepted (generated) variables
            count = count + 1
            ! We add the accepted x value to the output list of generated variables
            num(count) = x
        endif
    enddo

    ! Statistical study of the generated variables
    ev = 0.d0
    ev2 = 0.d0

    ! We compute the simple and quadratic expected values (as an average of all generated values)
    do i = 1, ndata
        ev = ev + num(i) / dble(ndata)
        ev2 = ev2 + (num(i))**2 / dble(ndata)
    enddo

    ! We compute the variance
    var = ev2 - ev**2
    ! We compute the standard deviation from the variance
    std = sqrt(var)

    ! We write the computed statistical values to the main program's text file
    write(1,*) "# Statistical values of the probability density p(x) computed with the acceptance-rejection method:"
    write(1,*) 
    write(1,*) "Expected value =", ev, "Variance =", var, "Standard deviation =", std
    write(1,*) 
    write(1,*) 

    return
end


subroutine uniforme(ndata, xlow, xhigh, num)
    implicit none
    integer :: ndata, i
    double precision :: xlow, xhigh, num(ndata)

    ! DESCRIPTION:
    ! This subroutine generates a uniform distribution of random numbers in a given interval

    ! VARIABLES:
    ! (input)
    ! ndata       --> Number of uniformly distributed numbers to be generated
    ! xlow, xhigh --> Interval in which the numbers are to be generated
    ! (output)
    ! num         --> Vector of dimension ndata carrying the generated numbers
    ! (internal)
    ! i           --> Counter


    do i = 1, ndata
        ! Variable change method for a uniform distribution
        num(i) = abs(xhigh - xlow) * rand() + xlow
    enddo

    return
end


subroutine exponential(ndata,xlam,num)
    implicit none
    integer :: ndata, i
    double precision :: xlam, num(ndata)

    ! DESCRIPTION:
    ! This subroutine generates, through a direct variable change of the type x --> lnt, a set of random variables distributed
    ! according to an exponential probability (with a given lambda) from values randomly distributed between [0,1).

    ! VARIABLES:
    ! (input)
    ! ndata --> Number of exponentially distributed numbers to be generated
    ! xlam  --> Value of the lambda parameter of the distribution
    ! (output)
    ! num   --> Vector of dimension ndata carrying the generated numbers
    ! (internal)
    ! i     --> Counter


    ! We generate ndata exponentially distributed values with the variable change from seed values uniformly distributed between 0 and 1
    do i = 1, ndata
        num(i) = (-1/xlam) * log(1.d0 - rand())
    enddo

    return
end


subroutine gaussian(ndata, val_mu, val_sigma, num)
    implicit none
    integer :: ndata, i
    double precision :: val_mu, val_sigma, pi, num(ndata)
    double precision :: x1, x2, r, theta, val_norm

    COMMON/DADES/pi

    ! DESCRIPTION:
    ! This subroutine generates a Gaussian distribution of random values using the Box-Müller method, adapting mu and sigma to our interest.

    ! VARIABLES:
    ! (input)
    ! ndata     --> Number of Gaussian-distributed numbers to be generated
    ! val_mu    --> Value of the mu parameter of the distribution (center)
    ! val_sigma --> Value of the sigma parameter of the distribution (deviation/width)
    ! (output)
    ! num       --> Vector of dimension ndata carrying the generated numbers
    ! (internal)
    ! pi        --> Value of the constant pi (imported from the main program)
    ! x1, x2    --> Variables carrying the base random numbers (distributed between [0,1))
    ! r, theta  --> Secondary variables that transform x1 and x2 into Gaussians
    ! val_norm  --> Value of each generated point before applying mu and sigma (i.e. normal distribution)
    ! i         --> Counter


    ! We start a loop of ndata
    do i = 1, ndata
        ! We generate the pair of random variables [0,1)
        x1 = rand()
        x2 = rand()

        ! We redefine x1 and x2 into r and theta so that they are Gaussian
        r = sqrt(-2 * log(x1))
        theta = 2*pi * x2

        ! We generate a Gaussian value with normal distribution
        val_norm = r * cos(theta)

        ! We transform the normal Box-Müller distribution into ours, with the given mu and sigma
        num(i) = val_mu + val_norm * val_sigma
    enddo

    return
end


! MONTE CARLO _______________________________________________________________________________________________________________________________________

subroutine montecarlo(N, valsx, func, dens, I, errI)
    implicit none
    integer :: N, j
    double precision :: valsx(N), func, dens
    double precision :: I, I2, errI

    ! DESCRIPTION:
    ! This subroutine computes the integral of a function from a set of values distributed according to a given density following the
    ! generalized Monte Carlo method. Thus, it can be applied to any case of Monte Carlo integration.

    ! VARIABLES:
    ! (input)
    ! N     --> Amount of numbers distributed according to the introduced probability density
    ! valsx --> List with the numbers distributed according to the probability density
    ! func  --> Function we want to integrate (1D function)
    ! dens  --> Normalized probability density in which the N introduced valsx values are distributed (1D function)
    ! (output)
    ! I     --> Final value of the integral
    ! errI  --> Error associated with the numerical calculation of the integral
    ! (internal)
    ! I2    --> Value of the integral squared, needed for the error calculation

    ! We initialize the values of the integral and its square to zero
    I = 0.d0
    I2 = 0.d0

    ! We execute a summation to compute the value of the integral and its square (useful for computing the associated error)
    do j = 1, N
        I = I + (1/dble(N)) * (func(valsx(j)) / dens(valsx(j)))
        I2 = I2 + (1/dble(N)) * (func(valsx(j)) / dens(valsx(j)))**2
    enddo

    ! Once the integral is computed, we determine the estimated error based on the number of values used
    errI = (1/sqrt(dble(N))) * sqrt(I2 - I**2)

    return
end


subroutine hitmiss(ndata, xlow, xhigh, sup, func, I, errI)
    implicit none
    integer :: ndata, nhit, i
    double precision :: xlow, xhigh, sup, func, I, errI
    double precision :: x1, x2, x, y

    ! DESCRIPTION:
    ! This subroutine executes a hit-or-miss method to numerically integrate a given function. It repeatedly generates random values
    ! appropriately scaled in a fixed region and approximates the value of the integral, based on the proportion between those that 
    ! fall under the curve of the function and the total amount.

    ! VARIABLES:
    ! (input)
    ! ndata       --> Total amount of random numbers generated over the fixed region
    ! xlow, xhigh --> Interval of abscissas that horizontally delimit the fixed region
    ! sup         --> Vertical limit of the fixed region (the closer to the maximum of the function to be integrated, the better)
    ! func        --> 1D function to be integrated
    ! (output)
    ! I           --> Final value of the integral
    ! errI        --> Error associated with the numerical calculation of the integral
    ! (internal)
    ! nhit        --> Number of values that fall under the curve (number of hits)
    ! x1, x2      --> Primary random variables (distributed between [0,1))
    ! x, y        --> Appropriately scaled coordinates of the generated values
    ! i           --> Counter

    ! We initialize the hit counter (values that land under the curve)
    nhit = 0

    ! We start a loop that generates ndata values (x,y) uniformly distributed within the range marked by the function curve
    do i = 1, ndata
        ! We generate two random values uniformly distributed between 0 and 1
        x1 = rand()
        x2 = rand()

        ! We scale the x coordinate from x1 so that it is uniformly distributed within the desired interval
        x = abs(xhigh - xlow) * x1 + xlow
        ! We scale the y coordinate from x2 so that it is uniformly distributed between 0 and the introduced maximum bound
        y = sup * x2

        ! When the randomly generated value (x,y) lies under the curve of the introduced function, we accept it
        if (y.le.f(x)) then
            ! We add 1 to the hit counter
            nhit = nhit + 1
        endif

    enddo

    ! We print to screen the record of generated and hit variables
    write(*,*) "Generated values =", ndata, "Hit values =", nhit

    ! We compute the value of the integral of the function from the proportion of hit values with respect to the generated ones
    I = sup * abs(xhigh - xlow) * nhit / dble(ndata)

    ! We compute the error associated with the numerical calculation of the integral
    errI = sup * abs(xhigh - xlow) / sqrt(dble(ndata)) * sqrt(nhit/dble(ndata) * (1.d0 - nhit/dble(ndata)))
    
    return
end