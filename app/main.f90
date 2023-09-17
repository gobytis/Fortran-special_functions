program main
    use special_functions
    implicit none

    print *, "beta"
    print *, beta(3, 4)
    print *, beta(3.0, 4.0)
    print *, beta(3.0d0, 4.0d0)


    print *, "polygamma k = -1", polygamma(-1, 1.0d0)
    print *, "k = 0"
    print *, "x=1.0d-10", polygamma(0, 1.0d-10)
    print *, "x=1.0d-5", polygamma(0, 1.0d-5)
    print *, "x=1.0d-2", polygamma(0, 1.0d-2)
    print *, "x=1.0d2", polygamma(0, 1.0d2)
    print *, "x=4.0d3", polygamma(0, 4.0d3)
    print *, "x=1.0d5", polygamma(0, 1.0d5)
    print *, "x=1.0d10", polygamma(0, 1.0d10)
    print *, "x=0", polygamma(0, 0.0d0)
    print *, "x<0", polygamma(0, -1.0d0)

    print *, "trigamma (k=1)"
    print *, "x=1.0d-10", polygamma(1, 1.0d-10)
    print *, "x=1.0d-5", polygamma(1, 1.0d-5)
    print *, "x=1.0d-2", polygamma(1, 1.0d-2)
    print *, "x=1.0d2", polygamma(1, 1.0d2)
    print *, "x=4.0d3", polygamma(1, 4.0d3)
    print *, "x=1.0d5", polygamma(1, 1.0d5)
    print *, "x=1.0d10", polygamma(1, 1.0d10)
    print *, "x=0", polygamma(1, 0.0d0)
    print *, "x<0", polygamma(1, -1.0d0)
end program main
