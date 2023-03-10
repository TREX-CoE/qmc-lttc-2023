subroutine uniform_montecarlo(a,nmax,energy)
  implicit none
  double precision, intent(in)  :: a
  integer*8       , intent(in)  :: nmax 
  double precision, intent(out) :: energy

  integer*8        :: istep
  double precision :: normalization, r(3), w, f

  double precision, external :: e_loc, psi

  energy = 0.d0
  normalization = 0.d0

  do istep = 1,nmax

     call random_number(r)
     r(:) = -5.d0 + 10.d0*r(:)

     f = psi(a,r)
     w = f*f

     energy = energy + w * e_loc(a,r)
     normalization = normalization + w

  end do

  energy = energy / normalization

end subroutine uniform_montecarlo

program qmc
  implicit none
  double precision, parameter :: a     = 1.2d0
  integer*8       , parameter :: nmax  = 100000
  integer         , parameter :: nruns = 30

  integer          :: irun
  double precision :: X(nruns)
  double precision :: ave, err

  do irun=1,nruns
     call uniform_montecarlo(a, nmax, X(irun))
  enddo

  call ave_error(X, nruns, ave, err)

  print *, 'E = ', ave, '+/-', err
end program qmc
