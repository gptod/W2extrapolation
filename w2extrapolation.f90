module extrapolation

    use para
    use ot_sinkhorn

    implicit none
    
    contains

subroutine extra_fortran(n,m,d,mu,nu,X,Y,f0,Z0,t,eps,deb,itmax,tol,tau,verb,f,g,Z,P,niter,errP,errG,obj)
    
    implicit none

    ! inputs
    integer :: n, m, d, itmax, verb
    real(wp) :: mu(n), nu(m), X(n,d), Y(m,d), f0(n), Z0(m,d)
    real(wp) :: t, eps, tol, tau
    logical :: deb
    !f2py intent(in) :: n, m, d, itmax, deb
    !f2py intent(in) :: mu, nu, X, Y, f0, Z0
    !f2py intent(in) :: t, eps, tol, tau, verb
    !f2py depend(n) :: mu, f0
    !f2py depend(m) :: nu
    !f2py depend(n,d) :: X
    !f2py depend(m,d) :: Y, Z0

    ! outputs
    integer :: niter
    real(wp) :: f(n), g(m), Z(m,d), P(n,m)
    real(wp) :: errP(itmax), errG(itmax), obj(itmax)
    !f2py intent(out) :: f, g, Z, P, niter, errP, errG, obj

    ! locals
    integer :: i
    real(wp) :: C(n,m), grad(m,d), norm_nu
    real(wp) :: grad_self(m,d), fin_self(m), fout_self(m), obj_self
    

    f=f0
    Z=Z0
    fin_self=0.0
    errP=0.0
    errG=0.0
    obj=0.0
    norm_nu=sqrt(sum(nu**2))
    niter=itmax

    do i=1,itmax
        
        call cost(n,m,d,X,Z,t,C)
        call mina_stab(n,m,C,f,mu,eps,g)
        call minb_stab(n,m,C,g,nu,eps,f)
        
        call get_P(n,m,f,g,C,mu,nu,eps,P)
        errP(i)=sqrt(sum((sum(P,dim=1)-nu)**2))/norm_nu

        call gradient(n,m,d,P,nu,X,Y,Z,t,grad)

        if (deb) then
            call gradient_deb_W2(m,d,nu,Z,fin_self,t,eps,itmax,tol,grad_self,fout_self,obj_self)
            fin_self=fout_self
        else
            grad_self=0.0
            obj_self=0.0
        end if
        grad=grad+grad_self
        Z=Z-tau*grad
        errG(i) = sqrt(sum(grad**2))

        ! approximate objective (for the true one I should recompute the plan with the new Z)
        obj(i)=-sum(f*mu)-sum(g*nu)+eps*sum(P)+sum(nu*sum((Z-Y)**2,2))/(2*(t-1))+obj_self
        
        if (verb==1) then
            write(*,'(A6,I6,2X,A13,E12.6,2X,A13,E12.6,2X,A11,E11.4)') &
            'Iter: ', i, 'Error marg.: ', errP(i), 'Error grad.: ', errG(i), 'Objective: ', obj(i)
        end if
        
        if (max(errP(i),errG(i))<tol) then
            niter=i
            exit
        end if

    end do

end subroutine extra_fortran

subroutine extra_nto1_fortran(Nm,nb,n,m,d,mu,nu,X,Y,f0,Z0,t,lm,eps,deb,itmax,tol,tau,verb,f,g,Z,P,niter,errP,errG,obj)
    
    implicit none

    ! inputs
    integer :: Nm, nb
    integer :: n(Nm), m, d, itmax, verb
    real(wp) :: mu(nb,Nm), nu(m), X(nb,d,Nm), Y(m,d), f0(nb,Nm), Z0(m,d)
    real(wp) :: t, lm(Nm), eps, tol, tau 
    logical :: deb
    !f2py intent(in) :: Nm, nb
    !f2py intent(in) :: n, m, d, itmax, deb
    !f2py intent(in) :: mu, nu, X, Y, f0, Z0
    !f2py intent(in) :: t, lm, eps, tol, tau, verb
    !f2py depend(Nm) :: n, lm
    !f2py depend(nb,Nm) :: mu, f0
    !f2py depend(m) :: nu
    !f2py depend(nb,d,Nm) :: X
    !f2py depend(m,d) :: Y, Z0

    ! outputs
    integer :: niter
    real(wp) :: f(nb,Nm), g(m,Nm), Z(m,d), P(nb,m,Nm)
    real(wp) :: errP(itmax), errG(itmax), obj(itmax)
    !f2py intent(out) :: f, g, Z, P, niter, errP, errG, obj

    ! locals
    integer :: i, j, k
    real(wp) :: C(nb,m,Nm), grad(m,d), norm_nu
    real(wp) :: errP_k(Nm), obj_k(Nm), grad_k(m,d,Nm), obj_kk
    real(wp) :: grad_self(m,d), fin_self(m), fout_self(m), obj_self
    
    C(:,:,:)=0.0
    f=f0
    Z=Z0
    fin_self=0.0
    errP=0.0
    errG=0.0
    obj=0.0
    norm_nu=sqrt(sum(nu**2))
    niter=itmax

    do i=1,itmax
        
        do k=1,Nm

            ! compute C_k
            call cost(n(k),m,d,X(1:n(k),:,k),Z,t/lm(k),C(1:n(k),:,k))

            ! compute lse_k
            call mina_stab(n(k),m,C(1:n(k),:,k),f(1:n(k),k),mu(1:n(k),k),eps,g(:,k))
            call minb_stab(n(k),m,C(1:n(k),:,k),g(:,k),nu,eps,f(1:n(k),k))

            ! compute P_k and err_k
            call get_P(n(k),m,f(1:n(k),k),g(:,k),C(1:n(k),:,k),mu(1:n(k),k),nu,eps,P(1:n(k),:,k))
            errP_k(k)=sqrt(sum((sum(P(1:n(k),:,k),dim=1)-nu)**2))/norm_nu

            ! compute grad_k
            call gradient_k(n(k),m,d,P(1:n(k),:,k),X(1:n(k),:,k),Z,t/lm(k),grad_k(:,:,k))

            ! approximate objective (for the true one I should recompute the plan with the new Z)
            obj_k(k)=-sum(f(1:n(k),k)*mu(1:n(k),k))-sum(g(:,k)*nu)+eps*sum(P(1:n(k),:,k))

        end do

        ! compute gradient step
        do j=1,m
            grad(j,:)=nu(j)*(Z(j,:)-Y(j,:))/(t-1)
        end do
        if (deb) then
            call gradient_deb_W2(m,d,nu,Z,fin_self,t,eps,itmax,tol,grad_self,fout_self,obj_self)
            fin_self=fout_self
        else
            grad_self=0.0
            obj_self=0.0
        end if
        grad=grad+grad_self
        obj_kk=0.0
        do k=1,Nm
            grad=grad+grad_k(:,:,k)
            errP(i)=errP(i)+errP_k(k)
            obj_kk=obj_kk+obj_k(k)
        end do
        Z=Z-tau*grad
        obj(i)=obj_kk+sum(nu*sum((Z-Y)**2,2))/(2*(t-1))+obj_self
        errG(i) = sqrt(sum(grad**2))

        if (verb==1) then
            write(*,'(A6,I6,2X,A13,E12.6,2X,A13,E12.6,2X,A11,E11.4)') &
            'Iter: ', i, 'Error marg.: ', errP(i), 'Error grad.: ', errG(i), 'Objective: ', obj(i)
        end if
        
        if (max(errP(i),errG(i))<tol) then
            niter=i
            exit
        end if

    end do

end subroutine extra_nto1_fortran

subroutine cost(n,m,d,X,Z,t,C)
    
    implicit none

    integer, intent(in) :: n, m, d
    real(wp), intent(in) :: X(n,d), Z(m,d), t
    real(wp), intent(out) :: C(n,m)

    integer :: i, j

    do i=1,n
        do j=1,m
            C(i,j)=0.5*sum((X(i,:)-Z(j,:))**2)/t
        end do
    end do

end subroutine

! subroutine mina_stab(n,m,C,f,mu,eps,g)

!     use para
!     implicit none

!     integer, intent(in) :: n, m
!     real(wp), intent(in) :: C(n,m), f(n), mu(n), eps

!     real(wp), intent(out) :: g(m)

!     integer :: j
!     real(wp) :: v(n), minv

!     do j=1,m
!         v=C(:,j)-f
!         minv=minval(v)
!         g(j)=-eps*log(sum(mu*exp(-(v-minv)/eps)))+minv
!     end do

! end subroutine

! subroutine minb_stab(n,m,C,g,nu,eps,f)

!     use para
!     implicit none

!     integer, intent(in) :: n, m
!     real(wp), intent(in) :: C(n,m), g(m), nu(m), eps

!     real(wp), intent(out) :: f(n)

!     integer :: i
!     real(wp) :: v(m), minv

!     do i=1,n
!         v=C(i,:)-g
!         minv=minval(v)
!         f(i)=-eps*log(sum(nu*exp(-(v-minv)/eps)))+minv
!     end do

! end subroutine

! subroutine get_P(n,m,f,g,C,mu,nu,eps,P)

!     use para
!     implicit none
    
!     integer, intent(in) :: n, m
!     real(wp), intent(in) :: f(n), g(m), C(n,m), mu(n), nu(m), eps

!     real(wp), intent(out) :: P(n,m)

!     integer :: i, j
    
!     do i=1,n
!         do j=1,m
!             P(i,j)=mu(i)*exp((f(i)+g(j)-C(i,j))/eps)*nu(j)
!         end do
!     end do

! end subroutine get_P

subroutine gradient(n,m,d,P,nu,X,Y,Z,t,grad)

    implicit none
    
    integer, intent(in) :: n, m, d
    real(wp), intent(in) :: P(n,m), nu(m), X(n,d), Y(m,d), Z(m,d), t

    real(wp), intent(out) :: grad(m,d)

    integer :: i, j
    real(wp) :: v(d)

    do j=1,m
        v=0.0
        do i=1,n
            v=v+(Z(j,:)-X(i,:))*P(i,j)
        end do
        grad(j,:)=nu(j)*(Z(j,:)-Y(j,:))/(t-1)-v/t
    end do

end subroutine gradient

subroutine gradient_k(n,m,d,P,X,Z,t,grad)

    implicit none
    
    integer, intent(in) :: n, m, d
    real(wp), intent(in) :: P(n,m), X(n,d), Z(m,d), t

    real(wp), intent(out) :: grad(m,d)

    integer :: i, j
    real(wp) :: v(d)

    do j=1,m
        v=0.0
        do i=1,n
            v=v+(Z(j,:)-X(i,:))*P(i,j)
        end do
        !grad(j,:)=nu(j)*(Z(j,:)-Y(j,:))/(t-1)-v/t
        grad(j,:)=-v/t
    end do

end subroutine gradient_k

end module extrapolation
