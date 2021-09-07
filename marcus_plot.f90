    Program main
!=============================================
! Integration of a function using Simpson rule 
!=============================================
implicit none
double precision f,g, a, b, integral
real*8 rf,rb
integer n, i, kk,atomic_time
double precision, parameter:: pi = 3.1415926
external f
external g

atomic_time=50

a = -0.1
b = 0.1

n = 100




do i=1,16
   call simpson(f,a,b,integral,n)
!   write (*,101) n, integral
   n = n*2
end do
rf=integral


n=100
do i=1,16
   call simpson1(g,a,b,integral,n)
!   write (*,101) n, integral
   n = n*2
end do

rb=integral

!write(*,*) rb,rf

call Marcus_plot(rf,rb)



end

  Function f(x)
!----------------------------------------
! Function for integration
!----------------------------------------
implicit none
double precision f, x,pi,Er,gama

pi=3.1416
Er=2000*((2.0E-04)**2)*(5.6097**2)/2
gama=1.0E-04
f=gama*(exp(x/9.5E-04)/(1+exp(x/9.5E-04)))*exp(-(Er-3.8E-03+x)**2/(4*Er*9.5E-04))/sqrt(4*pi*Er*9.5E-04)
!f=gama*(1/(1+exp(x/9.5E-04)))*exp(-(Er-3.8E-03-x)**2/(4*Er*9.5E-04))/sqrt(4*pi*Er*9.5E-04)

return
end

  Function g(x)
!----------------------------------------
! Function for integration
!----------------------------------------
implicit none
double precision g, x,pi,Er,gama

pi=3.1416
Er=2000*((2.0E-04)**2)*(5.6097**2)/2
gama=1.0E-04
!f=gama*(exp(x/9.5E-04)/(1+exp(x/9.5E-04)))*exp(-(Er+3.8E-03+x)**2/(4*Er*9.5E-04))/sqrt(4*pi*Er*9.5E-04)
g=gama*(1/(1+exp(x/9.5E-04)))*exp(-(Er+3.8E-03-x)**2/(4*Er*9.5E-04))/sqrt(4*pi*Er*9.5E-04)

return
end














 Subroutine simpson(f,a,b,integral,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================
implicit none
double precision f, a, b, integral,s
double precision h, x
integer nint
integer n, i

! if n is odd we add +1 to make it even
if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
s = 0.0
h = (b-a)/dfloat(n)
do i=2, n-2, 2
   x   = a+dfloat(i)*h
   s = s + 2.0*f(x) + 4.0*f(x+h)
end do
integral = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0
return
end subroutine simpson

 Subroutine simpson1(g,a,b,integral,n)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
! written by: Alex Godunov (October 2009)
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of integration
!==========================================================
implicit none
double precision g, a, b, integral,s
double precision h, x
integer nint
integer n, i

! if n is odd we add +1 to make it even
if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
s = 0.0
h = (b-a)/dfloat(n)
do i=2, n-2, 2
   x   = a+dfloat(i)*h
   s = s + 2.0*g(x) + 4.0*g(x+h)
end do
integral = (s + g(a) + g(b) + 4.0*g(a+h))*h/3.0
return
end subroutine simpson1


subroutine Marcus_plot(kf,kb)
implicit none        
double precision time,population
double precision Kt
integer :: tim,omega_t
real*8, intent(in) :: kf,kb

!write(*,*) kf,kb
population=0
kt=kf+kb
write(*,*) kf+kb
omega_t=100

do tim=0,5000*omega_t,100
   time=real(tim)
   population=exp(-(kt*time-log(kf)))/kt+(kb/kt)
   write(45,*) time,population
enddo   
end subroutine












