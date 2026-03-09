subroutine extra_fortran(n,m,d,mu,nu,X,Y,f0,Z0,t,eps,deb,itmax,tol,tau,verb,f,g,Z,P,niter,errP,errG,obj)
    
    use para
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
    real(wp) :: grad_self(m,d), f_self(m), obj_self
    

    f=f0
    Z=Z0
    f_self=0.0
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
            call gradient_deb_W2(m,d,nu,Z,t,eps,itmax,tol,grad_self,f_self,obj_self)
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
    
    use para
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
    real(wp) :: grad_self(m,d), f_self(m), obj_self
    
    C(:,:,:)=0.0
    f=f0
    Z=Z0
    f_self=0.0
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
            call gradient_deb_W2(m,d,nu,Z,t,eps,itmax,tol,grad_self,f_self,obj_self)
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
    
    use para
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

subroutine mina_stab(n,m,C,f,mu,eps,g)

    use para
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

    use para
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

    use para
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

subroutine gradient(n,m,d,P,nu,X,Y,Z,t,grad)

    use para
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

    use para
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

subroutine gradient_deb_W2(m,d,nu,Z,t,eps,itmax,tol,grad_self,f_self,obj_self)

    use para
    implicit none
    
    integer, intent(in) :: m, d, itmax
    real(wp), intent(in) :: nu(m), Z(m,d)
    real(wp), intent(in) :: t, eps, tol
    !f2py intent(in) :: m, d, itmax
    !f2py intent(in) :: nu, Z
    !f2py intent(in) :: t, eps, tol
    !f2py depend(m) :: nu
    !f2py depend(m,d) :: Z

    real(wp), intent(out) :: grad_self(m,d), obj_self
    real(wp), intent(inout) :: f_self(m)
    !f2py intent(out) :: grad_self, obj_self
    !f2py intent(inout) :: f_self

    ! local variables
    integer :: i, j
    real(wp) :: v(d), C_self(m,m), P_self(m,m)

    call cost_l2_self(m,d,Z,t,C_self)
    call w2eps_self(m,nu,C_self,eps,itmax,tol,0,f_self,P_self,obj_self)
    obj_self=0.5*obj_self

    do j=1,m
        v=0.0
        do i=1,m
            v=v+(Z(j,:)-Z(i,:))*P_self(i,j)
        end do
        grad_self(j,:)=v/t
    end do

end subroutine gradient_deb_W2

subroutine w2eps_self(n,mu,C,eps,itmax,tol,verb,f,P,obj_self)
    
    use para
    implicit none

    ! inputs
    integer, intent(in) :: n, itmax, verb
    real(wp), intent(in) :: mu(n), C(n,n)
    real(wp), intent(in) :: eps, tol
    !f2py intent(in) :: n, itmax
    !f2py intent(in) :: mu, nu, C, f0
    !f2py intent(in) :: eps, tol, verb
    !f2py depend(n) :: mu, C

    ! outputs
    real(wp), intent(out) :: P(n,n)
    real(wp), intent(out) :: obj_self
    real(wp), intent(inout) :: f(n)
    !f2py intent(out) :: P, obj_self
    !f2py intent(out) :: f

    ! locals
    integer :: i
    real(wp) :: fnew(n), norm_mu
    real(wp) :: errP(itmax), obj(itmax)
    
    errP=0.0
    obj=0.0
    norm_mu=sqrt(sum(mu**2))

    do i=1,itmax
        
        call mina_stab_self(n,C,f,mu,eps,fnew)
        f=0.5*f+0.5*fnew

        call get_P_self(n,f,C,mu,eps,P)
        errP(i)=sqrt(sum((sum(P,dim=1)-mu)**2))/norm_mu

        obj(i)=2*sum(f*mu)-eps*sum(P)

        if (verb==1) then
            print *, i, errP(i), obj(i)
        end if
        
        if ((errP(i)<tol).or.(i==itmax)) then
            obj_self=obj(i)
            exit
        end if

    end do

end subroutine w2eps_self

subroutine mina_stab_self(n,C,f,mu,eps,fnew)

    use para
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
    
    use para
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

    use para
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