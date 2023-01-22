subroutine pdmc(a, dt, nmax, energy, accep, tau, E_ref)
  implicit none
  double precision, intent(in)  :: a, dt, tau
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy, accep
  double precision, intent(in)  :: E_ref

  integer*8        :: istep
  integer*8        :: n_accep
  double precision :: sq_dt, chi(3), d2_old, prod, u
  double precision :: psi_old, psi_new, d2_new, argexpo, q
  double precision :: r_old(3), r_new(3)
  double precision :: d_old(3), d_new(3)
  double precision :: e, w, normalization, tau_current

  double precision, external :: e_loc, psi

  sq_dt = dsqrt(dt)

  ! Initialization
  energy  = 0.d0
  n_accep = 0_8
  normalization = 0.d0

  w           = 1.d0
  tau_current = 0.d0

  call random_gauss(r_old,3)

  call drift(a,r_old,d_old)
  d2_old  = d_old(1)*d_old(1) + &
            d_old(2)*d_old(2) + &
            d_old(3)*d_old(3)

  psi_old = psi(a,r_old)

  do istep = 1,nmax
     e = e_loc(a,r_old)
     w = w * dexp(-dt*(e - E_ref))

     normalization = normalization + w
     energy = energy + w*e
     
     tau_current = tau_current + dt

     ! Reset when tau is reached
     if (tau_current > tau) then
        w           = 1.d0
        tau_current = 0.d0
     endif

     call random_gauss(chi,3)
     r_new(:) = r_old(:) + dt*d_old(:) + chi(:)*sq_dt

     call drift(a,r_new,d_new)
     d2_new = d_new(1)*d_new(1) + &
              d_new(2)*d_new(2) + &
              d_new(3)*d_new(3)

     psi_new = psi(a,r_new)

     ! Metropolis
     prod = (d_new(1) + d_old(1))*(r_new(1) - r_old(1)) + &
            (d_new(2) + d_old(2))*(r_new(2) - r_old(2)) + &
            (d_new(3) + d_old(3))*(r_new(3) - r_old(3))

     argexpo = 0.5d0 * (d2_new - d2_old)*dt + prod

     q = psi_new / psi_old
     q = dexp(-argexpo) * q*q

     call random_number(u)

     if (u <= q) then

        n_accep = n_accep + 1_8

        r_old(:) = r_new(:)
        d_old(:) = d_new(:)
        d2_old   = d2_new
        psi_old  = psi_new

     end if

  end do

  energy = energy / normalization
  accep  = dble(n_accep) / dble(nmax)

end subroutine pdmc

program qmc
  implicit none
  double precision, parameter :: a     = 1.2d0
  double precision, parameter :: dt    = 0.05d0
  double precision, parameter :: E_ref = -0.5d0
  double precision, parameter :: tau   = 100.d0
  integer*8       , parameter :: nmax  = 100000
  integer         , parameter :: nruns = 30

  integer          :: irun
  double precision :: X(nruns), accep(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call pdmc(a, dt, nmax, X(irun), accep(irun), tau, E_ref)
  enddo

  call ave_error(X,nruns,ave,err)
  print *, 'E = ', ave, '+/-', err

  call ave_error(accep,nruns,ave,err)
  print *, 'A = ', ave, '+/-', err

end program qmc
