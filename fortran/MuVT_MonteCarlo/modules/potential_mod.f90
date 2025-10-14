!---------------------------------------------------------------
! Defines functions to compute interaction energies
! Scalable Lennard-Jones potential for N particles
!
! Define funciones para calcular energías de interacción
! Potencial de Lennard-Jones escalable para N partículas
!---------------------------------------------------------------
module potential_mod
    use iso_fortran_env, only: real64
    implicit none

    private
    public :: lj_energy
    public :: lj_energy_array

contains

    !---------------------------------------------------------------
    ! Computes the Lennard-Jones potential energy between two particles
    !
    ! Calcula la energía del potencial de Lennard-Jones entre dos partículas
    !---------------------------------------------------------------
    function lj_energy(r, epsilon, sigma) result(E)
        real(real64), intent(in) :: r
        real(real64), intent(in), optional :: epsilon, sigma
        real(real64) :: E
        real(real64) :: eps_val, sigma_val
        real(real64) :: r2, r6

        ! Default values / Valores por defecto
        eps_val   = 1.0_real64
        sigma_val = 1.0_real64
        if (present(epsilon)) eps_val = epsilon
        if (present(sigma))   sigma_val = sigma

        ! Avoid division by zero / Protección contra r=0
        if (r <= 1.0e-12_real64) then
            E = huge(1.0_real64)
        else
            ! Optimización de potencias
            r2 = (sigma_val/r)**2
            r6 = r2**3
            E = 4.0_real64*eps_val*(r6**2 - r6)
        end if
    end function lj_energy

!---------------------------------------------------------------
! Computes the total Lennard-Jones potential energy for N particles
! pos is 3 x N, consistent with system_mod
!
! Calcula la energía total del potencial de Lennard-Jones para N partículas
! pos es 3 x N, coherente con system_mod
!---------------------------------------------------------------
    function lj_energy_array(pos, epsilon, sigma) result(E_tot)
        real(real64), intent(in) :: pos(:, :)  ! 3 x N
        real(real64), intent(in), optional :: epsilon, sigma
        real(real64) :: E_tot
        integer :: i, j, N
        real(real64) :: dx, dy, dz, r
        real(real64) :: eps_val, sigma_val

        ! Number of particles / Número de partículas
        N = size(pos,2)

        ! Default values / Valores por defecto
        eps_val = merge(epsilon,1.0_real64,present(epsilon))
        sigma_val = merge(sigma,1.0_real64,present(sigma))

        E_tot = 0.0_real64
        do i = 1, N-1
            do j = i+1, N
                dx = pos(1,i) - pos(1,j)
                dy = pos(2,i) - pos(2,j)
                dz = pos(3,i) - pos(3,j)
                r = sqrt(dx*dx + dy*dy + dz*dz)
                E_tot = E_tot + lj_energy(r, eps_val, sigma_val)
            end do
        end do
    end function lj_energy_array

end module potential_mod
