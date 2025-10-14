!===============================================================
! Handles writing simulation data in binary (unformatted, stream) format.
!
! Features:
!   - Saves system properties: number of particles, energy, density, max displacement.
!   - Optionally writes Umbrella Sampling data.
!   - Optionally writes Metadynamics data.
!   - Writes particle positions in a compact binary format.
!
! Encargado de la escritura de datos de la simulación en formato
! binario (unformatted, stream).
!
! Características:
!   - Guarda información del sistema: número de partículas,
!     energía, densidad, desplazamiento máximo.
!   - Opcionalmente guarda datos de Umbrella Sampling.
!   - Opcionalmente guarda datos de Metadinámica.
!   - Escribe posiciones de partículas en formato binario compacto.
!===============================================================
module io_mod
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: write_binary

contains

    !---------------------------------------------------------------
    ! Writes all relevant μVT Monte Carlo simulation data to a binary file.
    !
    ! Arguments:
    !   file        - Output filename
    !   positions   - Particle positions (3 x N)
    !   n_particles - Number of particles
    !   energy      - Total system energy
    !   density     - Instantaneous density
    !   max_disp    - Maximum allowed displacement
    !
    ! Optional arguments:
    !   umbrella_on - Logical, enables Umbrella Sampling
    !   xi_target   - Control variable (if umbrella_on = .true.)
    !   k_bias      - Bias force constant
    !   meta_on     - Logical, enables Metadynamics
    !   sigma_meta  - Gaussian width for Metadynamics
    !   height_meta - Gaussian height for Metadynamics
    !   meta_bias   - Accumulated bias potential
    !
    ! Notes:
    !   - The file is opened as "unformatted stream" for system-independent binary I/O.
    !   - If Umbrella/Metadynamics are not active, .false. is written.
    !   - This binary format is compact but not portable across architectures with different endianness.
    !
    ! Escribe en un archivo binario toda la información relevante
    ! de la simulación Monte Carlo μVT.
    !
    ! Argumentos:
    !   file        - Nombre del archivo de salida
    !   positions   - Posiciones de las partículas (3 x N)
    !   n_particles - Número de partículas
    !   energy      - Energía total del sistema
    !   density     - Densidad instantánea
    !   max_disp    - Desplazamiento máximo permitido
    !
    ! Opcionales:
    !   umbrella_on - Lógico, activa Umbrella Sampling
    !   xi_target   - Variable de control (si umbrella_on = .true.)
    !   k_bias      - Constante de fuerza del sesgo
    !   meta_on     - Lógico, activa Metadinámica
    !   sigma_meta  - Anchura de gaussiana para metadinámica
    !   height_meta - Altura de gaussiana para metadinámica
    !   meta_bias   - Potencial de sesgo acumulado
    !
    ! Notas:
    !   - El archivo se abre en modo "unformatted stream", lo que
    !     asegura compatibilidad entre sistemas y lectura binaria directa.
    !   - Si Umbrella/Metadinámica no están activos, se escribe .false.
    !   - Este formato binario es compacto pero no portátil entre
    !     arquitecturas con distinto endianness.
    !---------------------------------------------------------------
    subroutine write_binary(file, positions, n_particles, energy, density, max_disp, &
                            umbrella_on, xi_target, k_bias, meta_on, sigma_meta, height_meta, meta_bias)
        character(len=*), intent(in) :: file
        real(real64), intent(in) :: positions(:,:)   ! Matriz de posiciones (3 x N)
        integer, intent(in) :: n_particles           ! Número de partículas
        real(real64), intent(in) :: energy           ! Energía total
        real(real64), intent(in) :: density          ! Densidad instantánea
        real(real64), intent(in) :: max_disp         ! Desplazamiento máximo
        logical, intent(in), optional :: umbrella_on, meta_on
        real(real64), intent(in), optional :: xi_target, k_bias
        real(real64), intent(in), optional :: sigma_meta, height_meta
        real(real64), intent(in), optional :: meta_bias(:)
        integer :: iunit

        ! Abrir archivo en modo binario, sobrescribiendo si ya existe
        open(newunit=iunit, file=file, status='replace', &
             form='unformatted', access='stream')

        !-----------------------------------------------------------
        ! Header with general system properties
        ! Cabecera con propiedades generales del sistema
        !-----------------------------------------------------------
        write(iunit) n_particles
        write(iunit) energy
        write(iunit) density
        write(iunit) max_disp

        !-----------------------------------------------------------
        ! Umbrella Sampling block
        ! Bloque de Umbrella Sampling
        !-----------------------------------------------------------
        if (present(umbrella_on)) then
            write(iunit) umbrella_on
            if (umbrella_on .and. present(xi_target) .and. present(k_bias)) then
                write(iunit) xi_target
                write(iunit) k_bias
            end if
        else
            write(iunit) .false.
        end if

        !-----------------------------------------------------------
        ! Metadynamics block
        ! Bloque de Metadinámica
        !-----------------------------------------------------------
        if (present(meta_on)) then
            write(iunit) meta_on
            if (meta_on .and. present(sigma_meta) .and. present(height_meta)) then
                write(iunit) sigma_meta
                write(iunit) height_meta
                if (present(meta_bias)) write(iunit) meta_bias
            end if
        else
            write(iunit) .false.
        end if

        !-----------------------------------------------------------
        ! Particle positions
        ! Posiciones de partículas
        !-----------------------------------------------------------
        write(iunit) positions(:,1:n_particles)

        ! Cerrar archivo
        close(iunit)
    end subroutine write_binary

end module io_mod
