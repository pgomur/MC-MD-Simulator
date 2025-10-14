!===============================================================
! Implements core Monte Carlo moves, acceptance tracking, Metadynamics,
! Replica Exchange, and related utilities for μVT simulations.
!
! Implementa movimientos Monte Carlo, seguimiento de aceptación,
! Metadinámica, Replica Exchange y utilidades relacionadas para simulaciones μVT.
!===============================================================
module mc_mod
    use iso_fortran_env, only: real64
    use system_mod
    use rng_mod
    use potential_mod
    implicit none

    private
    public :: monte_carlo_step, move_particle, insert_particle, delete_particle
    public :: adjust_max_disp, reset_counters, print_acceptance
    public :: set_metadynamics, update_metadynamics
    public :: init_replicas, exchange_replicas
    public :: meta_active, meta_bias, sigma_meta, height_meta

    !-----------------------------------------------------------
    ! Counters for the current window / Contadores para la ventana actual (reseteables)
    !-----------------------------------------------------------
    integer :: n_moves_attempted = 0, n_moves_accepted = 0
    integer :: n_insert_attempted = 0, n_insert_accepted = 0
    integer :: n_delete_attempted = 0, n_delete_accepted = 0

    !-----------------------------------------------------------
    ! Global counters for reporting / Contadores globales para reportes (no reseteables)
    !-----------------------------------------------------------
    integer :: tot_moves_attempted = 0, tot_moves_accepted = 0
    integer :: tot_insert_attempted = 0, tot_insert_accepted = 0
    integer :: tot_delete_attempted = 0, tot_delete_accepted = 0

    !-----------------------------------------------------------
    ! Metadynamics / Metadinámica
    !-----------------------------------------------------------
    logical :: meta_active = .false.
    real(real64), allocatable :: meta_bias(:)
    real(real64) :: sigma_meta = 0.0_real64
    real(real64) :: height_meta = 0.0_real64
    integer :: meta_steps = 0

    !-----------------------------------------------------------
    ! Replica Exchange Monte Carlo (REMC) / Intercambio de réplicas
    !-----------------------------------------------------------
    integer :: n_replicas = 1
    real(real64), allocatable :: replica_temps(:)
    real(real64), allocatable :: replica_energies(:)

contains

    !---------------------------------------------------------------
    ! Performs a single Monte Carlo step, selecting randomly among
    ! particle move, insertion, or deletion. Updates system state
    ! and applies Metadynamics bias if active.
    !
    ! Realiza un paso Monte Carlo, seleccionando aleatoriamente entre
    ! mover, insertar o eliminar una partícula. Actualiza el estado
    ! del sistema y aplica el sesgo de metadinámica si está activo.
    !
    ! Optional Arguments / Argumentos opcionales:
    !   replica_id - ID of the replica for REMC / ID de la réplica para REMC
    !---------------------------------------------------------------
    subroutine monte_carlo_step(replica_id)
        integer, intent(in), optional :: replica_id
        real(real64) :: r
        r = random_uniform()
        if (r < 0.7_real64) then
            call move_particle()
        else if (r < 0.85_real64) then
            call insert_particle()
        else
            call delete_particle()
        end if
        call update_system()

        if (meta_active) call update_metadynamics()
    end subroutine monte_carlo_step

    !---------------------------------------------------------------
    ! Performs a standard Monte Carlo particle displacement.
    ! Selects a random particle, proposes a trial displacement, and
    ! applies the Metropolis criterion to accept or reject the move.
    !
    ! Realiza un desplazamiento estándar de partícula Monte Carlo.
    ! Selecciona una partícula aleatoria, propone un desplazamiento de prueba
    ! y aplica el criterio de Metropolis para aceptar o rechazar el movimiento.
    !---------------------------------------------------------------
    subroutine move_particle()
        integer :: i
        real(real64) :: dr(3), E_old, E_new, deltaE, accept_prob

        if (n_particles == 0) return

        !-----------------------------------------------------------
        ! Update counters / Actualizar contadores
        !-----------------------------------------------------------
        n_moves_attempted = n_moves_attempted + 1
        tot_moves_attempted = tot_moves_attempted + 1

        !-----------------------------------------------------------
        ! Select particle and propose displacement / Seleccionar partícula y proponer desplazamiento
        !-----------------------------------------------------------
        i = int(random_uniform()*n_particles) + 1
        dr = (random_uniform() - 0.5_real64) * 2.0_real64 * max_disp

        E_old = local_energy(i)

        !-----------------------------------------------------------
        ! Apply trial move / Aplicar movimiento de prueba
        !-----------------------------------------------------------
        positions(:,i) = positions(:,i) + dr

        E_new = local_energy(i)
        deltaE = E_new - E_old

        !-----------------------------------------------------------
        ! Metropolis acceptance criterion / Criterio de aceptación de Metropolis
        !-----------------------------------------------------------
        accept_prob = exp(-deltaE / temperature)
        if (accept_prob > 1.0_real64) accept_prob = 1.0_real64

        if (random_uniform() < accept_prob) then
            ! accept
            n_moves_accepted = n_moves_accepted + 1
            tot_moves_accepted = tot_moves_accepted + 1
            ! update total energy (system_mod: energy_total variable expected)
            energy_total = energy_total + deltaE
        else
            ! reject: revert
            positions(:,i) = positions(:,i) - dr
        end if
    end subroutine move_particle

    !---------------------------------------------------------------
    ! Attempts a particle insertion in the grand-canonical ensemble.
    ! Samples a random position, computes interaction energy with
    ! existing particles, applies the Metropolis criterion, and
    ! updates positions, meta_bias, and system energy if accepted.
    !
    ! Intenta la inserción de una partícula en el ensamble gran-canónico.
    ! Genera una posición aleatoria, calcula la energía de interacción
    ! con las partículas existentes, aplica el criterio de Metropolis
    ! y actualiza posiciones, meta_bias y energía del sistema si se acepta.
    !---------------------------------------------------------------
    subroutine insert_particle()
        integer :: j, n_old
        real(real64) :: new_particle(3)
        real(real64), allocatable :: new_positions(:,:)
        real(real64), allocatable :: new_meta(:)
        real(real64) :: E_new, accept_prob
        real(real64) :: factor

        n_insert_attempted = n_insert_attempted + 1
        tot_insert_attempted = tot_insert_attempted + 1

        !-----------------------------------------------------------
        ! Sample particle uniformly in the simulation box / Generar posición aleatoria en la caja
        !-----------------------------------------------------------
        call random_number(new_particle)
        new_particle = new_particle * volume

        !-----------------------------------------------------------
        ! Compute interaction energy with existing particles / Calcular energía de interacción con partículas existentes
        !-----------------------------------------------------------
        E_new = 0.0_real64
        if (n_particles > 0) then
            ! Serial loop for reproducibility
            do j = 1, n_particles
                E_new = E_new + lj_energy(norm2(positions(:,j) - new_particle))
            end do
        end if

        !-----------------------------------------------------------
        ! Grand-canonical acceptance factor / Factor de aceptación gran-canónico
        ! Simplified, ignoring thermal wavelength / Simplificado, ignorando longitud de onda térmica
        !-----------------------------------------------------------
        factor = volume / real(n_particles + 1, real64)

        accept_prob = factor * exp(-(E_new - chemical_potential)/temperature)
        if (accept_prob > 1.0_real64) accept_prob = 1.0_real64

        !-----------------------------------------------------------
        ! Metropolis criterion: accept or reject / Criterio de Metropolis: aceptar o rechazar
        !-----------------------------------------------------------
        if (random_uniform() < accept_prob) then
            n_old = n_particles
            allocate(new_positions(3, n_old+1))
            if (n_old > 0) new_positions(:,1:n_old) = positions
            new_positions(:, n_old+1) = new_particle

            !-------------------------------------------------------
            ! Update meta_bias if Metadynamics is active / Actualizar meta_bias si metadinámica activa
            !-------------------------------------------------------
            if (meta_active) then
                if (allocated(meta_bias)) then
                    allocate(new_meta(n_old+1))
                    if (n_old > 0) new_meta(1:n_old) = meta_bias
                    new_meta(n_old+1) = 0.0_real64
                    call move_alloc(new_meta, meta_bias)
                else
                    allocate(meta_bias(n_old+1))
                    meta_bias = 0.0_real64
                end if
            end if

            call move_alloc(new_positions, positions)
            n_particles = n_old + 1
            n_insert_accepted = n_insert_accepted + 1
            tot_insert_accepted = tot_insert_accepted + 1

            !-------------------------------------------------------
            ! Update total system energy / Actualizar energía total del sistema
            !-------------------------------------------------------
            energy_total = energy_total + E_new
        end if
    end subroutine insert_particle

    !---------------------------------------------------------------
    ! Attempts a particle deletion in the grand-canonical ensemble.
    ! Selects a random particle, computes its interaction energy with
    ! the rest, applies the Metropolis criterion, and updates positions,
    ! meta_bias, and system energy if deletion is accepted.
    !
    ! Intenta la eliminación de una partícula en el ensamble gran-canónico.
    ! Selecciona una partícula aleatoria, calcula su energía de interacción
    ! con el resto, aplica el criterio de Metropolis y actualiza posiciones,
    ! meta_bias y energía del sistema si la eliminación es aceptada.
    !---------------------------------------------------------------
    subroutine delete_particle()
        integer :: i, j, n_old
        real(real64) :: E_old, accept_prob
        real(real64), allocatable :: tmp_positions(:,:)
        real(real64), allocatable :: tmp_bias(:)

        if (n_particles == 0) return

        n_delete_attempted = n_delete_attempted + 1
        tot_delete_attempted = tot_delete_attempted + 1

        i = int(random_uniform()*n_particles) + 1

        E_old = 0.0_real64
        if (n_particles > 1) then
            !$omp parallel do reduction(+:E_old) schedule(static)
            do j = 1, n_particles
                if (j /= i) E_old = E_old + lj_energy(norm2(positions(:,i) - positions(:,j)))
            end do
        end if

        !-----------------------------------------------------------
        ! Grand-canonical deletion acceptance factor / Factor de aceptación gran-canónico para eliminación
        !-----------------------------------------------------------
        accept_prob = real(n_particles, real64) / volume * exp((E_old - chemical_potential)/temperature)
        if (accept_prob > 1.0_real64) accept_prob = 1.0_real64

        if (random_uniform() < accept_prob) then
            n_old = n_particles
            if (n_old > 1) then
                allocate(tmp_positions(3, n_old-1))
                if (i > 1) tmp_positions(:,1:i-1) = positions(:,1:i-1)
                if (i < n_old) tmp_positions(:,i:n_old-1) = positions(:,i+1:n_old)
            else
                allocate(tmp_positions(3,0))
            end if

            !-------------------------------------------------------
            ! Adjust meta_bias if Metadynamics is active / Ajustar meta_bias si metadinámica activa
            !-------------------------------------------------------
            if (meta_active .and. allocated(meta_bias)) then
                if (n_old > 1) then
                    allocate(tmp_bias(n_old-1))
                    if (i > 1) tmp_bias(1:i-1) = meta_bias(1:i-1)
                    if (i < n_old) tmp_bias(i:n_old-1) = meta_bias(i+1:n_old)
                    call move_alloc(tmp_bias, meta_bias)
                else
                    if (allocated(meta_bias)) then
                        deallocate(meta_bias)
                    end if
                end if
            end if

            call move_alloc(tmp_positions, positions)
            n_particles = n_old - 1
            n_delete_accepted = n_delete_accepted + 1
            tot_delete_accepted = tot_delete_accepted + 1
            !-------------------------------------------------------
            ! Update total system energy / Actualizar energía total del sistema
            !-------------------------------------------------------
            energy_total = energy_total - E_old
        end if
    end subroutine delete_particle

    !---------------------------------------------------------------
    ! Computes the interaction energy of particle i with all other
    ! particles in the system using the Lennard-Jones potential.
    !
    ! Calcula la energía de interacción de la partícula i con todas
    ! las demás partículas del sistema usando el potencial de
    ! Lennard-Jones.
    !---------------------------------------------------------------
    function local_energy(i) result(E)
        integer, intent(in) :: i
        integer :: j
        real(real64) :: E
        E = 0.0_real64
        if (n_particles > 1) then
            ! Serial loop for reproducibility
            do j = 1, n_particles
                if (j /= i) E = E + lj_energy(norm2(positions(:,i) - positions(:,j)))
            end do
        end if
    end function local_energy

    !---------------------------------------------------------------
    ! Function: norm2
    ! Computes the Euclidean norm (magnitude) of a 3D vector r.
    !
    ! Función: norm2
    ! Calcula la norma euclidiana (magnitud) de un vector 3D r.
    !---------------------------------------------------------------
    function norm2(r) result(d)
        real(real64), intent(in) :: r(3)
        real(real64) :: d
        d = sqrt(sum(r**2))
    end function norm2

    !---------------------------------------------------------------
    ! Adjusts the maximum displacement (max_disp) for particle moves
    ! based on the acceptance ratio to maintain an optimal acceptance
    ! rate (~50%) in Monte Carlo simulations.
    !
    ! Ajusta el desplazamiento máximo (max_disp) para los movimientos
    ! de partículas basado en la tasa de aceptación, con el fin de
    ! mantener una aceptación óptima (~50%) en simulaciones Monte Carlo.
    !---------------------------------------------------------------
    subroutine adjust_max_disp()
        real(real64) :: ratio
        ratio = real(n_moves_accepted, real64) / max(1, n_moves_attempted)
        if (ratio > 0.55_real64) then
            max_disp = max_disp * 1.05_real64
        else if (ratio < 0.45_real64) then
            max_disp = max_disp * 0.95_real64
        end if
        call reset_counters()
    end subroutine adjust_max_disp

    !---------------------------------------------------------------
    ! Resets all windowed Monte Carlo attempt and acceptance counters
    ! to zero. Used after adjusting max_disp or at the start of a
    ! new measurement interval.
    !
    ! Reinicia todos los contadores de intentos y aceptaciones de
    ! Monte Carlo (ventana) a cero. Se utiliza después de ajustar
    ! max_disp o al inicio de un nuevo intervalo de medición.
    !---------------------------------------------------------------
    subroutine reset_counters()
        n_moves_attempted = 0
        n_moves_accepted  = 0
        n_insert_attempted = 0
        n_insert_accepted  = 0
        n_delete_attempted = 0
        n_delete_accepted  = 0
    end subroutine reset_counters

    !---------------------------------------------------------------
    ! Prints the Monte Carlo acceptance rates for moves, insertions,
    ! and deletions. Displays both the rates for the current window
    ! (since last reset) and the total accumulated rates.
    !
    ! Imprime las tasas de aceptación de Monte Carlo para movimientos,
    ! inserciones y eliminaciones. Muestra tanto las tasas del
    ! intervalo actual (desde el último reinicio) como las tasas
    ! acumuladas totales.
    !---------------------------------------------------------------
    subroutine print_acceptance()
        real(real64) :: rm_w, ri_w, rd_w
        real(real64) :: rm_t, ri_t, rd_t

        rm_w = real(n_moves_accepted, real64) / max(1, n_moves_attempted)
        ri_w = real(n_insert_accepted, real64) / max(1, n_insert_attempted)
        rd_w = real(n_delete_accepted, real64) / max(1, n_delete_attempted)

        rm_t = real(tot_moves_accepted, real64) / max(1, tot_moves_attempted)
        ri_t = real(tot_insert_accepted, real64) / max(1, tot_insert_attempted)
        rd_t = real(tot_delete_accepted, real64) / max(1, tot_delete_attempted)

        print *, "Acceptance rates (window):"
        print *, " Move:   ", rm_w
        print *, " Insert: ", ri_w
        print *, " Delete: ", rd_w

        print *, "Acceptance rates (total):"
        print *, " Move:   ", rm_t
        print *, " Insert: ", ri_t
        print *, " Delete: ", rd_t
    end subroutine print_acceptance

    !---------------------------------------------------------------
    ! Initializes the metadynamics bias potential for the Monte Carlo
    ! simulation. Allocates the meta_bias array and sets the Gaussian
    ! width (sigma) and height for the bias potential.
    !
    ! Inicializa el potencial de sesgo de metadinámica para la simulación
    ! Monte Carlo. Asigna el array meta_bias y establece la anchura
    ! (sigma) y altura de las gaussianas para el potencial de sesgo.
    !---------------------------------------------------------------
    subroutine set_metadynamics(sigma, height)
        real(real64), intent(in) :: sigma, height
        integer :: n_local

        sigma_meta = sigma
        height_meta = height
        meta_active = .true.

        n_local = n_particles
        if (n_local > 0) then
            if (.not. allocated(meta_bias)) then
                allocate(meta_bias(n_local))
                meta_bias = 0.0_real64
            else
                if (size(meta_bias) /= n_local) then
                    deallocate(meta_bias)
                    allocate(meta_bias(n_local))
                    meta_bias = 0.0_real64
                end if
            end if
        end if
        meta_steps = 0
    end subroutine set_metadynamics

    !---------------------------------------------------------------
    ! Updates the metadynamics bias potential by adding a Gaussian
    ! hill based on the current particle positions. The bias is
    ! centered on the average x-coordinate of all particles.
    !
    ! Actualiza el potencial de sesgo de metadinámica agregando
    ! una "colina" gaussiana basada en las posiciones actuales
    ! de las partículas. El sesgo se centra en la coordenada x
    ! promedio de todas las partículas.
    !---------------------------------------------------------------
    subroutine update_metadynamics()
        integer :: i, n_local
        real(real64), allocatable :: pos_snapshot(:,:)
        real(real64), allocatable :: delta(:)
        real(real64) :: center_x

        if (.not. meta_active) return
        n_local = n_particles
        if (n_local <= 0) return

        allocate(pos_snapshot(3, n_local))
        pos_snapshot(:,1:n_local) = positions(:,1:n_local)

        allocate(delta(n_local))
        delta = 0.0_real64

        center_x = 0.0_real64
        do i = 1, n_local
            center_x = center_x + pos_snapshot(1, i)
        end do
        center_x = center_x / real(n_local, real64)

        ! Compute Gaussian contribution for each particle / Calcular contribución gaussiana
        !$omp parallel do schedule(static) default(shared) private(i)
        do i = 1, n_local
            delta(i) = height_meta * exp( -0.5_real64 * ((pos_snapshot(1,i) - center_x)**2) / (sigma_meta**2) )
        end do

        if (.not. allocated(meta_bias)) then
            allocate(meta_bias(n_local))
            meta_bias = 0.0_real64
        else
            if (size(meta_bias) /= n_local) then
                deallocate(meta_bias)
                allocate(meta_bias(n_local))
                meta_bias = 0.0_real64
            end if
        end if

        do i = 1, n_local
            meta_bias(i) = meta_bias(i) + delta(i)
        end do

        deallocate(pos_snapshot)
        deallocate(delta)

        meta_steps = meta_steps + 1
    end subroutine update_metadynamics

    !---------------------------------------------------------------
    ! Initializes Replica Exchange Monte Carlo (REMC) replicas.
    ! Allocates temperature and energy arrays for each replica
    ! and sets initial energies to zero.
    !
    ! Inicializa réplicas para Monte Carlo con intercambio (REMC).
    ! Asigna arrays de temperaturas y energías para cada réplica
    ! y establece las energías iniciales en cero.
    !---------------------------------------------------------------
    subroutine init_replicas(nr, temps)
        integer, intent(in) :: nr
        real(real64), intent(in) :: temps(nr)
        if (allocated(replica_temps)) deallocate(replica_temps)
        if (allocated(replica_energies)) deallocate(replica_energies)
        n_replicas = nr
        allocate(replica_temps(nr))
        allocate(replica_energies(nr))
        replica_temps = temps
        replica_energies = 0.0_real64
    end subroutine init_replicas

    !---------------------------------------------------------------
    ! Attempts to exchange two REMC replicas using the Metropolis
    ! acceptance criterion based on temperatures and energies.
    !
    ! Intenta intercambiar dos réplicas REMC usando el criterio de
    ! aceptación de Metropolis basado en temperaturas y energías.
    !---------------------------------------------------------------
    subroutine exchange_replicas(i, j)
        integer, intent(in) :: i, j
        real(real64) :: d, p_accept
        d = (1.0_real64/replica_temps(i) - 1.0_real64/replica_temps(j)) * (replica_energies(j) - replica_energies(i))
        p_accept = min(1.0_real64, exp(d))
        if (random_uniform() < p_accept) then
            call swap_systems(i, j)
        end if
    end subroutine exchange_replicas

end module mc_mod
