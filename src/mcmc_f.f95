module dnIGEmcmc
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public :: runMCMC_f
    ! REAL, PRIVATE      :: zero = 0.0, half = 0.5, one = 1.0, two = 2.0,   &
    !                       vsmall = TINY(1.0), vlarge = HUGE(1.0)

contains


  subroutine runMCMC_f(cte, xvector, xvector_shifted) bind(C, name = "runMCMC_f_")
      !real(kind = c_double), intent(in)               :: mean, sd
      integer(kind = c_int), intent(in)               :: cte

      real(kind = c_double), intent(in), dimension(100) :: xvector
      real(kind = c_double), intent(out), dimension(100) :: xvector_shifted

      !llc = 0.0_c_double
      print *, cte/100.0d0
      xvector_shifted = xvector + 666.0_c_double


  end subroutine runMCMC_f


  ! subroutine runMCMC_f(N, result) bind(C, name = "runMCMC_f_")
  !     integer(kind = c_int), intent(in)               :: N
  !     real(kind = c_double), intent(out), dimension(N) :: result
  !
  !     call dnIGE(N, result)
  !
  ! end subroutine runMCMC_f

  ! subroutine dnIGE(N, result)
  !   integer, intent(in):: N
  !   double precision, intent(out), dimension(N) :: result
  !   integer :: j
  !
  !   do j=1, N
  !     result(j) = random_normal(mean = 0.d0, sd = 1.d0)
  !   end do
  ! end subroutine dnIGE





!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                           functions from dnIGEmodule
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!
!
! function betaSum(N,sires,sire,ID,s,group,tau,Ag,Af,eg,ef,mu_g,mu_f,inf,T,Vmu_g)
!   !Arguments
!   integer, intent(in):: N,sires
!   integer, intent(in), dimension(N) :: group,ID,s,sire,inf
!   double precision, intent(in), dimension(N) :: tau
!   double precision, intent(in), dimension(sires) :: Ag,Af
!   double precision, intent(in), dimension(N) :: eg,ef
!   !double precision, dimension(1):: phi_Ag
!   double precision :: betaSum,mu_g,mu_f,T, Vmu_g
!   ! locals
!   integer :: i,j
!   double precision, dimension(N,1):: sum1,sum2
!   logical, dimension(N):: cond_k
!   sum1=0.d0
!   animal: do j=1,N
!   	nidx: if (s(j)==1) then      !likelihood evaluated for non-index cases only
!      cond_k=(group==group(j) .AND. tau .LE. tau(j)) !previous infecteds from same group
!       if (inf(j)==1) then     !infecteds
!        	sum1(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((tau(j)-tau(pack(ID,cond_k)))*&
!                   exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!       else                   !non infecteds
!           sum1(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((T-tau(pack(ID,cond_k)))&
!                     *exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!       end if
!   	end if nidx
!   end do animal
!   betaSum=sum(sum1)
!   !print *, betaSum
! end function
!
!
!
! !--------------------------------------------------------------------------------!
! !                   Conditional for BVs (based on likelihood)                    !
! !                    Assumes vector ordering, that is j=id(j), j=1,N             !
! !--------------------------------------------------------------------------------!
!       !lAg_cand(j)=cond_Ag(N,sires,sire,ID,s,group,tau_true,Ag_cand,Af_cur,eg_cur,ef_cur,&
!                       ! mu_g,mu_f,phi_Ag,Ag_cand(j),id_sire(j),inf,T)  !id_sire(j)
! !susceptibility
! function cond_Ag(N,sires,sire,ID,s,group,tau,Ag,Af,eg,ef,mu_g,mu_f,phi_Ag,Ag_j,sire_i,inf,T,beta)
!   !Arguments
!   integer, intent(in):: N,sires,sire_i
!   integer, intent(in), dimension(N) :: group,ID,s,sire,inf
!   double precision, intent(in), dimension(N) :: tau
!   double precision, intent(in), dimension(sires) :: Ag,Af
!   double precision, intent(in), dimension(N) :: eg,ef
!   double precision, dimension(1):: phi_Ag
!   double precision :: cond_Ag,mu_g,mu_f,Ag_j,T,beta
!   ! locals
!   integer :: i,j
!   double precision, dimension(N,1):: sum1,sum2
!   logical, dimension(N):: cond_k
!   sum1=0.d0
!   sum2=0.d0
!   animal: do j=1,N
!    l_g: if (sire(j)==sire_i .AND. s(j)==1) then      !non-index case offspring of sire i
!       cond_k=(group==group(j) .AND. tau .LT. tau(j)) !previous infecteds from same group
!       if (inf(j)==1) then     !infecteds
!        	sum1(j,1)=Ag(sire(j))
!       	sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((tau(j)-tau(pack(ID,cond_k)))*&
!                   exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!       else                   !non infecteds
!           sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((T-tau(pack(ID,cond_k)))&
!                     *exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!      end if
!   end if l_g
!   end do animal
!   cond_Ag=sum(sum1)-beta*sum(sum2)-(phi_Ag(1)*0.5d0)*(Ag_j**2)
!  ! print *, cond_Ag
! end function
!
!
! !Infectivity
! !lAf_cand(j)=cond_Af(N,sires,sire,ID,s,group,tau_true,Ag_cur,Af_cand,eg_cur,ef_cur,&
! !                 mu_g,mu_f,phi_Af,ngr,ng_off(j),Af_cand(j),id_sire(j),inf,T)  !id_sire(j)
! function cond_Af(N,sires,sire,ID,s,group,tau,Ag,Af,eg,ef,mu_g,mu_f, &
!                  phi_Af,ngr,ng_off_i,Af_j,sire_i,inf,T,beta)
!   !Arguments
!   integer, intent(in):: N,sires,ng_off_i,sire_i
!   integer, intent(in), dimension(N) :: group,ID,s,sire,inf
!   integer, intent(in), dimension(ng_off_i) :: ngr
!   double precision, intent(in), dimension(N) :: tau
!   double precision, intent(in), dimension(sires) :: Ag,Af
!   double precision, intent(in), dimension(N) :: eg,ef
!   double precision, dimension(1):: phi_Af
!   double precision :: cond_Af,mu_g,mu_f,Af_j,T,beta
!   ! locals
!   integer :: i,j
!   double precision, dimension(N,1):: sum1,sum2
!   logical, dimension(N):: cond_k
!   sum1=0.d0
!   sum2=0.d0
!   animal: do j=1,N
!    if (any(group(j)==ngr .AND. s(j)==1)) then  !if animal j is in a group that has some offspring from sire_i
!      if (tau(j) .GT. minval(tau(pack(ID,sire==sire_i .AND. group==group(j))))) then  !if animal j was infected after any offspring of sire_i in its group
!        cond_k=(group==group(j) .AND. tau .LT. tau(j))
!        if (inf(j)==1) then     ! infecteds
!        sum1(j,1)=log(sum(exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k)))))
!         sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((tau(j)-tau(pack(ID,cond_k)))*&
!                   exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!         else                   ! non-infecteds
!             sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((T-tau(pack(ID,cond_k)))&
!                       *exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!        end if
!       end if
!     end if
!   end do animal
!  cond_Af=sum(sum1)-beta*sum(sum2)-(phi_Af(1)*0.5d0)*(Af_j**2)
!   !print *, cond_Af
! end function
!
! !--------------------------------------------------------------------------------!
! !      Conditional for environmental effects (based on likelihood)               !
! !                    Assumes vector ordering, that is j=id(j), j=1,N             !
! !--------------------------------------------------------------------------------!
! !    leg_cur(j)=cond_Ag(N,sires,sire,ID,j,s,group,tau_true,Ag_cur,Af_cur,eg_cur(j),ef_cur,&
! !                        mu_g,mu_f,phi_eg,inf,T)
!
! !susceptibility
! function cond_eg(N,sires,sire,ID,id_j,s,group,tau,Ag,Af,eg_j,ef,mu_g,mu_f,phi_eg,inf,T,beta)
!   !Arguments
!   integer, intent(in):: N,id_j,sires
!   integer, intent(in), dimension(N) :: group,ID,s,sire,inf
!   double precision, intent(in), dimension(N) :: tau
!   double precision, intent(in), dimension(sires) :: Ag,Af
!   double precision, intent(in), dimension(N) ::ef
!   double precision, dimension(1):: phi_eg
!   double precision :: cond_eg,mu_g,mu_f,T, eg_j,beta
!   ! locals
!   !integer :: i,j
!   double precision :: sum1,sum2
!   logical, dimension(N):: cond_k
!   sum1=0.d0
!   sum2=0.d0
!       cond_k=(group==group(id_j) .AND. tau .LT. tau(id_j)) !previous infecteds from same group
!       if (inf(id_j)==1) then     ! infected animal
!         sum1=eg_j
!         sum2=exp(mu_g+Ag(sire(id_j))+eg_j)*sum((tau(id_j)-tau(pack(ID,cond_k)))*&
!                   exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!       else                       ! non-infected animal
!           sum2=exp(mu_g+Ag(sire(id_j))+eg_j)*sum((T-tau(pack(ID,cond_k)))&
!                     *exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!      end if
!   ! end do animal
!   cond_eg=sum1-beta*sum2-(phi_eg(1)*0.5d0)*(eg_j**2)
! !  print *, cond_eg
! end function
!
! !infectivity
! function cond_ef(N,sires,sire,ID,id_j,s,group,tau,Ag,Af,eg,ef,mu_g,mu_f, &
!                  phi_ef,inf,T,beta)
!   !Arguments
!   integer, intent(in):: N,sires,id_j
!   integer, intent(in), dimension(N) :: group,ID,s,sire,inf
!   double precision, intent(in), dimension(N) :: tau
!   double precision, intent(in), dimension(sires) :: Ag,Af
!   double precision, intent(in), dimension(N) :: eg,ef
!   double precision, dimension(1):: phi_ef
!   double precision :: cond_ef,mu_g,mu_f,T,beta
!   ! locals
!   integer :: i,j
!   double precision, dimension(N,1):: sum1,sum2
!   logical, dimension(N):: cond_k
!   sum1=0.d0
!   sum2=0.d0
!   animal: do j=1,N
!      if (tau(j) .GT. tau(id_j) .AND. group(j)==group(id_j)) then  !if animal was infected after id_j in its group
!         cond_k=(group==group(j) .AND. tau .LT. tau(j)) !infected prior to animal j
!        if (inf(j)==1) then     ! infected animal
!         sum1(j,1)=log(sum(exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k)))))
!         sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((tau(j)-tau(pack(ID,cond_k)))*&
!                   exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!         else                   ! non-infecteds
!             sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((T-tau(pack(ID,cond_k)))&
!                       *exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!        end if
!       end if
!   end do animal
!  cond_ef=sum(sum1)-beta*sum(sum2)-(phi_ef(1)*0.5d0)*(ef(id_j)**2)
!  ! print *, cond_ef
! end function
!
!
! !tau
! function condTau(N,sires,sire,ID,id_j,s,group,tau,Ag,Af,eg,ef,mu_g,mu_f, &
!                  phi_ef,inf,T,beta)
!   !Arguments
!   integer, intent(in):: N,sires,id_j
!   integer, intent(in), dimension(N) :: group,ID,s,sire,inf
!   double precision, intent(in), dimension(N) :: tau
!   double precision, intent(in), dimension(sires) :: Ag,Af
!   double precision, intent(in), dimension(N) :: eg,ef
!   double precision, dimension(1):: phi_ef
!   double precision :: condTau,mu_g,mu_f,T,beta
!   ! locals
!   integer :: i,j
!   double precision, dimension(N,1):: sum1,sum2
!   logical, dimension(N):: cond_k
!   sum1=0.d0
!   sum2=0.d0
!   animal: do j=1,N
!      if (tau(j) .GT. tau(id_j) .AND. group(j)==group(id_j)) then  !if animal was infected after id_j in its group
!         cond_k=(group==group(j) .AND. tau .LT. tau(j)) !infected prior to animal j
!        if (inf(j)==1) then     ! infected animal
!         sum1(j,1)=log(sum(exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k)))))
!         sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((tau(j)-tau(pack(ID,cond_k)))*&
!                   exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!         else                   ! non-infecteds
!             sum2(j,1)=exp(mu_g+Ag(sire(j))+eg(j))*sum((T-tau(pack(ID,cond_k)))&
!                       *exp(mu_f+Af(sire(pack(ID,cond_k)))+ef(pack(ID,cond_k))))
!        end if
!       end if
!   end do animal
!  condTau=sum(sum1)-beta*sum(sum2)
!  ! print *, condTau
! end function
!
!
!
!
!
! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
! !                           RV generation functions
! !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
! function random_normal(mean, sd)
!   !ADAPTED FROM ALAN MILLER's random_normal() (random.f90)
!     double precision ,intent(in)           :: mean, sd
!     double precision             :: random_normal
!
!     !     Local variables
!     double precision      :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
!                 r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
!
!     !    Generate P = (u,v) uniform in rectangle enclosing acceptance region
!
!     DO
!       CALL RANDOM_NUMBER(u)
!       CALL RANDOM_NUMBER(v)
!       v = 1.7156 * (v - 0.5)
!
!     !     Evaluate the quadratic form
!       x = u - s
!       y = ABS(v) - t
!       q = x**2 + y*(a*y - b*x)
!
!     !     Accept P if inside inner ellipse
!       IF (q < r1) EXIT
!     !     Reject P if outside outer ellipse
!       IF (q > r2) CYCLE
!     !     Reject P if outside acceptance region
!       IF (v**2 < -4.0*LOG(u)*u**2) EXIT
!     END DO
!
!     !     Return ratio of P's coordinates as the normal deviate
!     ! random_normal = (v/u)
!     random_normal = (v/u) * sd + mean
! end function
!
!
! FUNCTION random_gamma(s, scale, first) RESULT(fn_val)
!
! ! Adapted from Fortran 77 code from the book:
! !     Dagpunar, J. 'Principles of random variate generation'
! !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! !     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
! !     CALLS EITHER random_gamma1 (S > 1.0)
! !     OR random_exponential (S = 1.0)
! !     OR random_gamma2 (S < 1.0).
!
! !     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).
!
! double precision, INTENT(IN)    :: s, scale
! LOGICAL, INTENT(IN) :: first
! REAL                :: fn_val
! !     Local variables
!
!
! IF (s <= zero) THEN
!   WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
!   STOP
! END IF
!
! IF (s > one) THEN
!   fn_val = random_gamma1(s, first)
! ELSE IF (s < one) THEN
!   fn_val = random_gamma2(s, first)
! ELSE
!   fn_val = random_exponential()
! END IF
!
! fn_val = fn_val*scale
! RETURN
! END FUNCTION random_gamma
!
!
!
! FUNCTION random_gamma1(s, first) RESULT(fn_val)
!
!
! ! Uses the algorithm in
! ! Marsaglia, G. and Tsang, W.W. (2000) `A simple method for generating
! ! gamma variables', Trans. om Math. Software (TOMS), vol.26(3), pp.363-372.
!
! ! Generates a random gamma deviate for shape parameter s >= 1.
!
! double precision, INTENT(IN)    :: s
! LOGICAL, INTENT(IN) :: first
! REAL                :: fn_val
!
! ! Local variables
! REAL, SAVE  :: c, d
! REAL        :: u, v, x
!
! IF (first) THEN
!   d = s - one/3.
!   c = one/SQRT(9.0*d)
! END IF
!
! ! Start of main loop
! DO
!
! ! Generate v = (1+cx)^3 where x is random normal; repeat if v <= 0.
!
!   DO
!     x = random_normal(mean = 0.d0, sd = 1.d0)
!     v = (one + c*x)**3
!     IF (v > zero) EXIT
!   END DO
!
! ! Generate uniform variable U
!
!   CALL RANDOM_NUMBER(u)
!   IF (u < one - 0.0331*x**4) THEN
!     fn_val = d*v
!     EXIT
!   ELSE IF (LOG(u) < half*x**2 + d*(one - v + LOG(v))) THEN
!     fn_val = d*v
!     EXIT
!   END IF
! END DO
!
! RETURN
! END FUNCTION random_gamma1
!
!
!
! FUNCTION random_gamma2(s, first) RESULT(fn_val)
!
! ! Adapted from Fortran 77 code from the book:
! !     Dagpunar, J. 'Principles of random variate generation'
! !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! ! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! ! GAMMA2**(S-1) * EXP(-GAMMA2),
! ! USING A SWITCHING METHOD.
!
! !    S = SHAPE PARAMETER OF DISTRIBUTION
! !          (REAL < 1.0)
!
! double precision, INTENT(IN)    :: s
! LOGICAL, INTENT(IN) :: first
! REAL                :: fn_val
!
! !     Local variables
! REAL       :: r, x, w
! REAL, SAVE :: a, p, c, uf, vr, d
!
! IF (s <= zero .OR. s >= one) THEN
!   WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
!   STOP
! END IF
!
! IF (first) THEN                        ! Initialization, if necessary
!   a = one - s
!   p = a/(a + s*EXP(-a))
!   IF (s < vsmall) THEN
!     WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
!     STOP
!   END IF
!   c = one/s
!   uf = p*(vsmall/a)**s
!   vr = one - vsmall
!   d = a*LOG(a)
! END IF
!
! DO
!   CALL RANDOM_NUMBER(r)
!   IF (r >= vr) THEN
!     CYCLE
!   ELSE IF (r > p) THEN
!     x = a - LOG((one - r)/(one - p))
!     w = a*LOG(x)-d
!   ELSE IF (r > uf) THEN
!     x = a*(r/p)**c
!     w = x
!   ELSE
!     fn_val = zero
!     RETURN
!   END IF
!
!   CALL RANDOM_NUMBER(r)
!   IF (one-r <= w .AND. r > zero) THEN
!     IF (r*(w + one) >= one) CYCLE
!     IF (-LOG(r) <= w) CYCLE
!   END IF
!   EXIT
! END DO
!
! fn_val = x
! RETURN
!
! END FUNCTION random_gamma2
!
! FUNCTION random_exponential() RESULT(fn_val)
!
! ! Adapted from Fortran 77 code from the book:
! !     Dagpunar, J. 'Principles of random variate generation'
! !     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
!
! ! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! ! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! ! TO EXP(-random_exponential), USING INVERSION.
!
! REAL  :: fn_val
!
! !     Local variable
! REAL  :: r
!
! DO
!   CALL RANDOM_NUMBER(r)
!   IF (r > zero) EXIT
! END DO
!
! fn_val = -LOG(r)
! RETURN
!
! END FUNCTION random_exponential
!
!


end module dnIGEmcmc
