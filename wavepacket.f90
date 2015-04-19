program proj6
IMPLICIT NONE
INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
INTEGER :: j,nt, f
INTEGER, PARAMETER :: jmax=1001, ntmax=1000
COMPLEX(KIND=Prec14), DIMENSION(0:jmax-1) :: a,b,c,chi
COMPLEX(KIND=Prec14), DIMENSION(0:jmax-1) :: psi
COMPLEX(KIND=Prec14) :: ar,ai
REAL(KIND=Prec14), PARAMETER :: L=200.0d0,sigma=4.0d0
REAL(KIND=Prec14), PARAMETER :: x0=L/2.0d0,k0=2.0d0
REAL(KIND=Prec14), PARAMETER :: dx=L/jmax
REAL(KIND=Prec14), PARAMETER :: dt=0.1d0
REAL(KIND=Prec14), PARAMETER :: alpha=dt/(2.0d0*dx**2)
REAL(KIND=Prec14) :: x,sspi,pi, p, p_sum

open(10, file='wave0.dat')
open(11, file='wave1.dat')
open(12, file='wave2.dat')
open(13, file='wave3.dat')
open(20, file='prob.dat')

pi = 4.0d0*DATAN(1.0d0)
sspi = pi**(1/4)
ai=(0.0d0,1.0d0) ! i=sqrt(-1) (imaginary
ar=(1.0d0,0.0d0) ! 1 (real)

	a=-ai*alpha/4.0d0
	b=0.50d0*(ar+ai*alpha)
	c=-ai*alpha/4.0d0

do j=0,jmax-1
	x=dx*j
	psi(j)=(ar/(dsqrt(sigma)*sspi))*dexp(-(x-x0)**2/(2.0d0*sigma**2))*(ar*dcos(k0*x)+ai*dsin(k0*x))
	write(10,*) j, dreal(psi(j))
enddo

do nt=1, ntmax
	p = 0.0d0
	p_sum = 0.0d0
	call tridag(a,b,c,psi,chi,jmax)
	psi = chi - psi
	do j = 0, jmax - 1
		!integration for probability (should be equal to unity since integrated over all x)
		p = psi(j+1)*conjg(psi(j+1)) + psi(j)*conjg(psi(j))
		p_sum = p_sum + (p/2.0d0)*dt
		if((nt == 500) .or. (nt == ntmax))then
			f = 10 + int(2*(real(nt)/real(ntmax)))
			write(f,*) j, dreal(psi(j))
		endif	
	end do
		
	write(20,*) nt, p_sum
enddo

end program proj6



SUBROUTINE TRIDAG(A,B,C,R,U,N)
INTEGER :: N,J
INTEGER, PARAMETER :: NMAX=10000
INTEGER, PARAMETER :: DP= SELECTED_REAL_KIND(14)
COMPLEX(KIND=DP) :: BET
COMPLEX(KIND=DP), DIMENSION(NMAX) :: GAM
COMPLEX(KIND=DP), DIMENSION(N) :: A,B,C,R,U
IF(B(1).EQ.0.)PAUSE
BET=B(1)
U(1)=R(1)/BET
DO J=2,N
GAM(J)=C(J-1)/BET
BET=B(J)-A(J)*GAM(J)
IF(BET.EQ.0.)PAUSE
U(J)=(R(J)-A(J)*U(J-1))/BET
ENDDO
DO J=N-1,1,-1
U(J)=U(J)-GAM(J+1)*U(J+1)
ENDDO
end subroutine TRIDAG
