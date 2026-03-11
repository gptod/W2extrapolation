module ot_sinkhorn
    
    use para

    implicit none

    contains

subroutine ote(n,m,mu,nu,C,f0,eps,itmax,tol,verb,f,g,P,niter,errP,obj)
    
    implicit none

    ! inputs
    integer :: n, m, itmax, verb
    real(wp) :: mu(n), nu(m), C(n,m), f0(n)
    real(wp) :: eps, tol
    !f2py intent(in) :: n, m, itmax
    !f2py intent(in) :: mu, nu, C, f0
    !f2py intent(in) :: eps, tol, verb
    !f2py depend(n) :: mu, f0
    !f2py depend(m) :: nu
    !f2py depend(n,m) :: C

    ! outputs
    integer :: niter
    real(wp) :: f(n), g(m), P(n,m)
    real(wp) :: errP(itmax), obj(itmax)
    !f2py intent(out) :: f, g, P, niter, errP, obj

    ! locals
    integer :: i
    real(wp) :: norm_nu
    

    do i=1,n
        f(i)=f0(i)
    end do
    errP=0.0
    obj=0.0
    norm_nu=sqrt(sum(nu**2))
    niter=itmax

    do i=1,itmax
        
        call mina_stab(n,m,C,f,mu,eps,g)
        call minb_stab(n,m,C,g,nu,eps,f)
        
        call get_P(n,m,f,g,C,mu,nu,eps,P)
        errP(i)=sqrt(sum((sum(P,dim=1)-nu)**2))/norm_nu

        obj(i)=sum(f*mu)+sum(g*nu)-eps*sum(P)+eps

        if (verb==1) then
            print *, i, errP(i), obj(i)
        end if
        
        if (errP(i)<tol) then
            niter=i
            exit
        end if

    end do

end subroutine ote

subroutine cost_l2(n,m,d,X,Z,C)
    
    implicit none

    integer, intent(in) :: n, m, d
    real(wp), intent(in) :: X(n,d), Z(m,d)
    real(wp), intent(out) :: C(n,m)

    !f2py intent(in) :: n, m, d, X, Z
    !f2py depend(n,d) :: X
    !f2py depend(m,d) :: Z
    !f2py intent(out) :: C
    !f2py depend(n,m) :: C

    integer :: i, j

    do i=1,n
        do j=1,m
            C(i,j)=0.5*sum((X(i,:)-Z(j,:))**2)
        end do
    end do

end subroutine

subroutine cost_l1(n,m,d,X,Z,C)
    
    use para
    implicit none

    integer, intent(in) :: n, m, d
    real(wp), intent(in) :: X(n,d), Z(m,d)
    real(wp), intent(out) :: C(n,m)

    !f2py intent(in) :: n, m, d, X, Z
    !f2py depend(n,d) :: X
    !f2py depend(m,d) :: Z
    !f2py intent(out) :: C
    !f2py depend(n,m) :: C

    integer :: i, j

    do i=1,n
        do j=1,m
            C(i,j)=0.5*sqrt(sum((X(i,:)-Z(j,:))**2))
        end do
    end do

end subroutine

subroutine mina_stab(n,m,C,f,mu,eps,g)

    implicit none

    integer, intent(in) :: n, m
    real(wp), intent(in) :: C(n,m), f(n), mu(n), eps

    real(wp), intent(out) :: g(m)

    integer :: j
    real(wp) :: v(n), minv

    do j=1,m
        v=C(:,j)-f
        minv=minval(v)
        g(j)=-eps*log(sum(mu*exp(-(v-minv)/eps)))+minv
    end do

end subroutine

subroutine minb_stab(n,m,C,g,nu,eps,f)

    implicit none

    integer, intent(in) :: n, m
    real(wp), intent(in) :: C(n,m), g(m), nu(m), eps

    real(wp), intent(out) :: f(n)

    integer :: i
    real(wp) :: v(m), minv

    do i=1,n
        v=C(i,:)-g
        minv=minval(v)
        f(i)=-eps*log(sum(nu*exp(-(v-minv)/eps)))+minv
    end do

end subroutine

subroutine get_P(n,m,f,g,C,mu,nu,eps,P)

    implicit none
    
    integer, intent(in) :: n, m
    real(wp), intent(in) :: f(n), g(m), C(n,m), mu(n), nu(m), eps

    real(wp), intent(out) :: P(n,m)

    integer :: i, j
    
    do i=1,n
        do j=1,m
            P(i,j)=mu(i)*exp((f(i)+g(j)-C(i,j))/eps)*nu(j)
        end do
    end do

end subroutine get_P

subroutine gradient_deb_W2(m,d,nu,Z,fin_self,t,eps,itmax,tol,grad_self,fout_self,obj_self)

    implicit none
    
    integer, intent(in) :: m, d, itmax
    real(wp), intent(in) :: nu(m), Z(m,d), fin_self(m)
    real(wp), intent(in) :: t, eps, tol
    !f2py intent(in) :: m, d, itmax
    !f2py intent(in) :: nu, Z, fin_self
    !f2py intent(in) :: t, eps, tol
    !f2py depend(m) :: nu, fin_self
    !f2py depend(m,d) :: Z

    real(wp), intent(out) :: grad_self(m,d), obj_self
    real(wp), intent(inout) :: fout_self(m)
    !f2py intent(out) :: grad_self, obj_self
    !f2py intent(out) :: fout_self

    ! local variables
    integer :: i, j, niter_self
    real(wp) :: v(d), C_self(m,m), P_self(m,m), errP_self

    call cost_l2_self(m,d,Z,t,C_self)
    call ote_self(m,nu,C_self,fin_self,eps,itmax,tol,0,fout_self,P_self,niter_self,errP_self,obj_self)
    obj_self=0.5*obj_self

    do j=1,m
        v=0.0
        do i=1,m
            v=v+(Z(j,:)-Z(i,:))*P_self(i,j)
        end do
        grad_self(j,:)=v/t
    end do

end subroutine gradient_deb_W2

subroutine ote_self(n,mu,C,fin_self,eps,itmax,tol,verb,fout_self,P,niter_self,errP_self,obj_self)
    
    implicit none

    ! inputs
    integer, intent(in) :: n, itmax, verb
    real(wp), intent(in) :: mu(n), C(n,n), fin_self(n)
    real(wp), intent(in) :: eps, tol
    !f2py intent(in) :: n, itmax
    !f2py intent(in) :: mu, nu, C, fin_self
    !f2py intent(in) :: eps, tol, verb
    !f2py depend(n) :: mu, C, fin_self

    ! outputs
    integer, intent(out) :: niter_self
    real(wp), intent(out) :: fout_self(n), P(n,n)
    real(wp), intent(out) :: errP_self, obj_self
    !f2py intent(out) :: fout_self, P, niter_self, errp_self, obj_self

    ! locals
    integer :: i
    real(wp) :: fnew(n), norm_mu
    
    fout_self=fin_self
    errP_self=0.0
    obj_self=0.0
    norm_mu=sqrt(sum(mu**2))
    niter_self=itmax

    do i=1,itmax
        
        call mina_stab_self(n,C,fout_self,mu,eps,fnew)
        fout_self=0.5*fout_self+0.5*fnew

        call get_P_self(n,fout_self,C,mu,eps,P)
        errP_self=sqrt(sum((sum(P,dim=1)-mu)**2))/norm_mu

        obj_self=2*sum(fout_self*mu)-eps*sum(P)+eps

        if (verb==1) then
            print *, i, errP_self, obj_self
        end if
        
        if (errP_self<tol) then
            niter_self=i
            exit
        end if

    end do

end subroutine ote_self

subroutine mina_stab_self(n,C,f,mu,eps,fnew)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in) :: C(n,n), mu(n), eps
    real(wp), intent(in) :: f(n)

    real(wp), intent(out) :: fnew(n)

    integer :: i
    real(wp) :: v(n), minv

    do i=1,n
        v = C(:,i) - f
        minv = minval(v)
        fnew(i) = -eps * log(sum(mu * exp(-(v-minv)/eps))) + minv
    end do

end subroutine

subroutine cost_l2_self(n,d,X,t,C)
    
    implicit none

    integer, intent(in) :: n, d
    real(wp), intent(in) :: X(n,d), t
    !f2py intent(in) :: n, d
    !f2py intent(in) :: X
    !f2py depend(n,d) :: X

    real(wp), intent(out) :: C(n,n)
    !f2py intent(out) :: C

    integer :: i, j

    do i=1,n
        do j=1,n
            C(i,j)=0.5*sum((X(i,:)-X(j,:))**2)/t
        end do
    end do

end subroutine

subroutine get_P_self(n,f,C,mu,eps,P)

    implicit none

    integer, intent(in) :: n
    real(wp), intent(in) :: f(n), C(n,n), mu(n), eps
    real(wp), intent(out) :: P(n,n)

    integer :: i,j

    do i=1,n
        do j=1,n
            P(i,j) = mu(i) * mu(j) * &
                     exp((f(i)+f(j)-C(i,j))/eps)
        end do
    end do

end subroutine

end module ot_sinkhorn