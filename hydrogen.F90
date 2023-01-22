double precision function potential(r)
  implicit none
  double precision, intent(in) :: r(3)

  double precision             :: distance

  distance = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) )

  if (distance > 0.d0) then
     potential = -1.d0 / distance
  else
     print *, 'Warning: potential at r=0.d0 diverges'
     potential = -huge(1.d0)
  end if

end function potential

subroutine test_potential
    implicit none
    double precision :: r(3)
    double precision :: expected_output
    double precision, external :: potential

    expected_output = -1.d0/15.d0

    r(:) = (/ 2.d0, 5.d0, 14.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 5.d0, 14.d0, 2.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ -2.d0, 5.d0, -14.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 5.d0, -14.d0, -2.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 0.d0, 9.d0, 12.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = (/ 9.d0, -12.d0, 0.d0 /)
    if (potential(r) /= expected_output) stop 'Failed'

    r(:) = 0.d0
    expected_output = -huge(1.d0)
    if (potential(r) /= expected_output) stop 'Failed r=0'
    print *, 'potential ok'
    
end subroutine test_potential

#ifdef TEST_H
program test_h
  call test_potential
end program test_h
#endif

double precision function psi(a, r)
  implicit none
  double precision, intent(in) :: a, r(3)

  psi = dexp(-a * dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) ))
end function psi

double precision function kinetic(a,r)
  implicit none
  double precision, intent(in) :: a, r(3)

  double precision             :: distance

  distance = dsqrt( r(1)*r(1) + r(2)*r(2) + r(3)*r(3) ) 

  if (distance > 0.d0) then
     kinetic =  a * (1.d0 / distance - 0.5d0 * a)
  else
     print *, 'Warning: kinetic energy diverges at r=0'
     kinetic =  a * (huge(1.d0) - 0.5d0 * a)
  end if

end function kinetic

double precision function e_loc(a,r)
  implicit none
  double precision, intent(in) :: a, r(3)

  double precision, external :: kinetic
  double precision, external :: potential

  e_loc = kinetic(a,r) + potential(r)

end function e_loc

subroutine drift(a,r,b)
  implicit none
  double precision, intent(in)  :: a, r(3)
  double precision, intent(out) :: b(3)

  double precision :: ar_inv

  ar_inv = -a / dsqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
  b(:)   = r(:) * ar_inv

end subroutine drift
