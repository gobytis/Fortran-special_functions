!> ポリガンマ関数の実装です。
!> This is a modern fortran implementation of a following algorithm.
!> [Reference]
!> 石岡恒憲 (1993) ポリ・ガンマ関数の C 言語，および Fortran 77 言語による算譜.
!>     応用統計学, 22(1), 23--37.
!>     Ishioka T. (1993) Plygamma functions subroutine programmed in C and Fortran 77
!> 
submodule (special_functions) polygamma_function
    use, intrinsic :: ieee_arithmetic
    implicit none

    integer, parameter :: i_bernoulli_trm = 10
    real(dp), parameter :: bernoulli_number(i_bernoulli_trm) = [ &
    1.6666666666666666d-1, 3.3333333333333333d-2, &
    2.3809523809523809d-2, 3.3333333333333333d-2, &
    7.5757575757575757d-2, 2.5311355311355311d-1, &
    1.6666666666666667d0,  7.0921568627450980d0,  &
    5.4971177944862155d1,  5.2912424242424242d2   &
    ]


    contains


    !> ディガンマ関数用のラッパーです。
    impure module function dp_digamma(x) result(retval)
        !> 入力
        real(dp), intent(in) :: x
        !> 戻り値
        real(dp) :: retval
        integer, parameter :: k = 0

        retval = dp_polygamma(k, x)
    end function


    !> トリガンマ関数用のラッパーです。
    impure module function dp_trigamma(x) result(retval)
        !> 入力
        real(dp), intent(in) :: x
        !> 戻り値
        real(dp) :: retval
        integer, parameter :: k = 1

        retval = dp_polygamma(k, x)
    end function


    !> ディガンマ関数の値 (psi^k(x)) を返します。
    impure recursive module function dp_polygamma(k, x) result(retval)
        !> 関数パラメータ
        !> k = 0, 1, 2, corresponds to digamma, trigamma, tetragamma.... respectively.
        integer, intent(in) :: k
        !> 入力
        real(dp), intent(in) :: x
        !> 戻り値
        real(dp) :: retval


        if(k < 0 .or. x < 0.0_dp) then
            retval = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if(x == 0.0_dp) then
            retval = ieee_value(0.0_dp, ieee_positive_inf)
            return
        end if

        block
            real(dp) :: large_number
            real(dp) :: k_fractorial

            large_number = get_large_number(k)
            k_fractorial = gamma(dble(k + 1))

            if(x >= large_number) then
                retval = get_polygamma_for_large_x(k, x)
            else
                retval = get_polygamma_for_small_x(k, x, large_number)
            end if
        end block
    end function

    
    impure module function get_large_number(k) result(large_number)
        integer, intent(in) :: k
        real(dp) :: large_number
        
        real(dp), parameter :: default_value = 13.06_dp
        real(dp), parameter :: power = 1.0_dp / 18.0_dp
        real(dp), parameter :: f_coeff = 174611.0_dp / 55.0_dp
        real(dp), parameter :: large_number_coeff = 6.81921_dp


        if(k < 0) then
            large_number = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if(k <= 3) then
            large_number = default_value
            return
        end if

        block
            real(dp) :: f
            integer :: i
            
            do i = 21, (k + 9)
                f = f * dble(i)
            end do
            do i = 3, (k + 1)
                f = f / dble(i)
            end do

            f = f * f_coeff
            large_number = large_number_coeff * f ** power
            if(large_number < default_value) then
                large_number = default_value
            end if
        end block
    end function 


    impure module function get_polygamma_for_large_x(k, x) result(retval)
        integer, intent(in) :: k
        real(dp), intent(in) :: x
        real(dp) :: retval
        integer :: i

        if(k == 0) then
            retval = dp_digamma_for_large_x(x)
        else
            retval = dp_high_order_polygamma_for_large_x(k, x)
        end if
    end function


    impure module function get_polygamma_for_small_x(k, x, large_value) result(polygamma)
        integer, intent(in) :: k
        real(dp), intent(in) :: x
        real(dp), intent(in) :: large_value
        real(dp) :: polygamma

        real(dp) :: y
        real(dp) :: y_powered
        integer :: isgn
        integer :: n
        integer :: i

        n = int(large_value - x)
        y = dble(n) + x + 1.0_dp

        if(x < 0.0 .and. y == dble(int(y))) then
            polygamma = ieee_value(0.0_dp, ieee_positive_inf)
            return
        end if

        polygamma = dp_polygamma(k, y)
    
        if(mod(k, 2) == 0) then
            isgn = -1
        else
            isgn = 1
        end if

        block
            real(dp) :: k_fractorial

            k_fractorial = gamma(dble(k + 1))

            do i = 1, n+1
                y = y - 1.0_dp
                if(abs(y) < 1.0d-3) then
                    if(x > 0.0d0) then
                        y = x - dble(int(x + 0.5))
                    else
                        y = x - dble(int(x - 0.5))
                    end if
                end if

                y_powered = y**k
                if(y_powered * y == 0.0d0) then
                    polygamma = ieee_value(0.0_dp, ieee_negative_inf)
                    return
                end if

                polygamma = polygamma + dble(isgn) * k_fractorial / y_powered / y
            end do
        end block
    end function


    pure module function dp_digamma_for_large_x(x) result(digamma)
        real(dp), intent(in) :: x
        real(dp) :: digamma

        real(dp) :: x_squared
        integer :: isgn
        integer :: ii
        integer :: i2
        integer :: i

        x_squared = x * x
        digamma = 0.0_dp
        isgn = 1
        do i = 1, i_bernoulli_trm
            ii = i_bernoulli_trm - i + 1
            i2 = ii + ii
            digamma = digamma + bernoulli_number(ii) / dble(i2 * isgn)
            digamma = digamma / x_squared
            isgn = - isgn
        end do
        digamma = digamma + log(x) - 0.5_dp / x
    end function


    pure module function dp_high_order_polygamma_for_large_x(k, x) result(polygamma)
        integer, intent(in) :: k
        real(dp), intent(in) :: x
        real(dp) :: polygamma

        real(dp) :: c
        real(dp) :: x_squared
        integer :: isgn
        integer :: i, j

        x_squared = x * x
        polygamma = 0.0_dp
        c = 1.0_dp

        if(mod(k, 2) == 0) then
            isgn = 1
        else
            isgn = -1
        end if

        block
            integer :: ii
            integer :: i2
            do i = 1, i_bernoulli_trm
                ii = i_bernoulli_trm - i + 1
                i2 = ii + ii
                do j = i2 + 1, i2 + k - 1
                    c = c * dble(j)
                end do

                polygamma = polygamma + bernoulli_number(ii) * c * dble(isgn)
                polygamma = polygamma / x_squared
                isgn = - isgn
            end do
        end block

        do i = 1, k
            polygamma = polygamma / x
        end do
        
        block
            real(dp) :: x_powered
            real(dp) :: k_fractorial

            x_powered = x ** k
            k_fractorial = gamma(dble(k + 1))
            polygamma = polygamma - 0.5_dp * k_fractorial / x_powered / x * dble(isgn)
            c = gamma(dble(k))
            polygamma = polygamma - c / x_powered * dble(isgn)
        end block
    end function

end submodule
