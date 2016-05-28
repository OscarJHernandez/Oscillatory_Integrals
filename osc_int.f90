program oscillatory_integrals
implicit none
real(8),parameter::a=0.d0 ! the starting point of the integration
integer,parameter:: integr = 2 ! for Sin(w*x) function
real(8),parameter::epsabs = 1.0d-14
real(8)::abserr
integer::ier
integer::neval
real(8)::result
real(8)::r,dh
real(8),external::A11
integer::i
real(8),parameter::hc =197.327

dh = 0.1d0

do i=0,1000
r = dh*dfloat(i)
call qawf(A11,a,r/(hc),integr,epsabs, result, abserr, neval, ier )
print *, r,result
end do

end program


real(8) function A11(p)
implicit none
real(8)::p
real(8),external::regulator
real(8),parameter::mpi =140.d0

A11 = (p*regulator(p))/(p**2+ mpi**2)**2

end function

real(8) function regulator(p)
implicit none
real(8)::p
real(8),parameter::Lambda=10.d0

regulator = dexp(-(p/Lambda)**2)
end function
