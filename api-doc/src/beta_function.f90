!> ベータ関数の実装です。
submodule (special_functions) beta_function
    implicit none

    contains

    !> 整数に対する beta 関数の値を返します。
    !> `beta(a, b)` = $int_0^1 t^{a-1} (1 - t)^{a - 1} dt$
    pure module function int_beta(a, b) result(retval)
        !> argument 1
        integer, intent(in) :: a
        !> argument 2
        integer, intent(in) :: b
        !> 戻り値
        real(dp) :: retval
        real(dp) :: dp_a
        real(dp) :: dp_b

        dp_a = dble(a)
        dp_b = dble(b)

        retval = exp(log_gamma(dp_a) + log_gamma(dp_b) - log_gamma(dp_a + dp_b))
    end function


    !> 単精度実数に対する beta 関数の値を返します。
    !> `beta(a, b)` = $int_0^1 t^{a-1} (1 - t)^{a - 1} dt$
    pure module function sp_beta(a, b) result(retval)
        !> argument 1
        real(sp), intent(in) :: a
        !> argument 2
        real(sp), intent(in) :: b
        !> 戻り値
        real(sp) :: retval

        retval = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function


    !> 倍精度実数に対する beta 関数の値を返します。
    !> `beta(a, b)` = $int_0^1 t^{a-1} (1 - t)^{a - 1} dt$
    pure module function dp_beta(a, b) result(retval)
        !> argument 1
        real(dp), intent(in) :: a
        !> argument 2
        real(dp), intent(in) :: b
        !> 戻り値
        real(dp) :: retval

        retval = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
    end function

end submodule
