!===============================================================
! Random number generation utilities for MC/MD simulations
!
! Utilidades de generación de números aleatorios para simulaciones MC/MD
!===============================================================
module rng_mod
    use iso_fortran_env, only: real64
    implicit none

    private
    public :: set_seed
    public :: random_uniform
    public :: random_uniform_range
    public :: random_gaussian
    public :: random_gaussian_truncated
    public :: random_vector
    public :: random_vector_gaussian
    public :: random_vector_gaussian_truncated
    public :: random_matrix_gaussian

    !---------------------------------------------------------------
    ! Private variables for Box-Muller
    ! Variables privadas para Box-Muller
    !---------------------------------------------------------------
    logical :: has_spare = .false.
    real(real64) :: spare_gauss = 0.0_real64

contains

    !---------------------------------------------------------------
    ! Initializes the random number generator with a given seed
    ! If no seed is provided, generates one from the system clock
    !
    ! Inicializa el generador de números aleatorios con una semilla dada
    ! Si no se proporciona semilla, se genera a partir del reloj del sistema
    !---------------------------------------------------------------
    subroutine set_seed(seed)
        integer, intent(in), optional :: seed
        integer :: i, n
        integer, allocatable :: seed_array(:)
        integer :: s
        integer :: t, count

        call random_seed(size = n)
        allocate(seed_array(n))

        if (present(seed)) then
            s = mod(seed, 2147483647)
            if (s == 0) s = 1
        else
            call system_clock(count, t)
            s = mod(abs(int(count * 1234567 + int(1000*real(t)))), 2147483647)
            if (s == 0) s = 1
        end if

        do i = 1, n
            seed_array(i) = mod(s + i - 1, 2147483647)
            if (seed_array(i) == 0) seed_array(i) = 1
        end do

        call random_seed(put = seed_array)
        deallocate(seed_array)
        print *, "RNG initialized with seed =", s
    end subroutine set_seed


    !---------------------------------------------------------------
    ! Returns a uniformly distributed random number in [0,1)
    !
    ! Devuelve un número aleatorio uniformemente distribuido en [0,1)
    !---------------------------------------------------------------
    function random_uniform() result(r)
        real(real64) :: r
        call random_number(r)
    end function random_uniform

    !---------------------------------------------------------------
    ! Returns a uniformly distributed random number in [a,b)
    !
    ! Devuelve un número aleatorio uniformemente distribuido en [a,b)
    !---------------------------------------------------------------
    function random_uniform_range(a,b) result(r)
        real(real64), intent(in) :: a, b
        real(real64) :: r
        if (b <= a) error stop "random_uniform_range: b debe ser mayor que a"
        r = a + (b-a)*random_uniform()
    end function random_uniform_range


    !---------------------------------------------------------------
    ! Generates a Gaussian (normal) random number using Box-Muller
    !
    ! Genera un número aleatorio gaussiano (normal) usando Box-Muller
    !---------------------------------------------------------------
    function random_gaussian(mean, sigma) result(r)
        real(real64), intent(in) :: mean, sigma
        real(real64) :: r, u1, u2
        if (has_spare) then
            r = spare_gauss
            has_spare = .false.
        else
            call random_number(u1)
            call random_number(u2)
            r = sqrt(-2.0_real64*log(u1)) * cos(2.0_real64*3.141592653589793*u2)
            spare_gauss = sqrt(-2.0_real64*log(u1)) * sin(2.0_real64*3.141592653589793*u2)
            has_spare = .true.
        end if
        r = r*sigma + mean
    end function random_gaussian


    !---------------------------------------------------------------
    ! Generates a Gaussian random number truncated to [min_val,max_val]
    !
    ! Genera un número aleatorio gaussiano truncado en [min_val,max_val]
    !---------------------------------------------------------------
    function random_gaussian_truncated(mean, sigma, min_val, max_val) result(r)
        real(real64), intent(in) :: mean, sigma, min_val, max_val
        real(real64) :: r
        if (max_val <= min_val) error stop "random_gaussian_truncated: max_val debe ser mayor que min_val"
        do
            r = random_gaussian(mean,sigma)
            if (r >= min_val .and. r <= max_val) exit
        end do
    end function random_gaussian_truncated


    !---------------------------------------------------------------
    ! Generates a vector of n uniform random numbers in [0,1)
    !
    ! Genera un vector de n números aleatorios uniformes en [0,1)
    !---------------------------------------------------------------
    function random_vector(n) result(vec)
        integer, intent(in) :: n
        real(real64), allocatable :: vec(:)
        allocate(vec(n))
        call random_number(vec)
    end function random_vector


    !---------------------------------------------------------------
    ! Generates a vector of n Gaussian random numbers
    ! Vectorized using Box-Muller for efficiency
    !
    ! Genera un vector de n números aleatorios gaussianos
    ! Vectorizado usando Box-Muller para eficiencia
    !---------------------------------------------------------------
    function random_vector_gaussian(n, mean, sigma) result(vec)
        integer, intent(in) :: n
        real(real64), intent(in) :: mean, sigma
        real(real64), allocatable :: vec(:)
        integer :: i
        allocate(vec(n))
        do i = 1, n
            vec(i) = random_gaussian(mean, sigma)
        end do
    end function random_vector_gaussian


    !---------------------------------------------------------------
    ! Generates a vector of n truncated Gaussian random numbers
    !
    ! Genera un vector de n números aleatorios gaussianos truncados
    !---------------------------------------------------------------
    function random_vector_gaussian_truncated(n, mean, sigma, min_val, max_val) result(vec)
        integer, intent(in) :: n
        real(real64), intent(in) :: mean, sigma, min_val, max_val
        real(real64), allocatable :: vec(:)
        integer :: i
        if (max_val <= min_val) error stop "random_vector_gaussian_truncated: max_val debe ser mayor que min_val"
        allocate(vec(n))
        do i = 1, n
            vec(i) = random_gaussian_truncated(mean,sigma,min_val,max_val)
        end do
    end function random_vector_gaussian_truncated


    !---------------------------------------------------------------
    ! Generates a dim1 x dim2 matrix of Gaussian random numbers
    ! Vectorized using random_vector_gaussian
    !
    ! Genera una matriz dim1 x dim2 de números aleatorios gaussianos
    ! Vectorizada usando random_vector_gaussian
    !---------------------------------------------------------------
    function random_matrix_gaussian(dim1, dim2, mean, sigma) result(mat)
        integer, intent(in) :: dim1, dim2
        real(real64), intent(in) :: mean, sigma
        real(real64), allocatable :: mat(:,:)
        integer :: i
        allocate(mat(dim1,dim2))
        do i = 1, dim1
            mat(i,:) = random_vector_gaussian(dim2, mean, sigma)
        end do
    end function random_matrix_gaussian

end module rng_mod
