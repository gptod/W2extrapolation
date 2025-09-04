subroutine extra_fortran(n,m,d,mu,nu,X,Y,f0,Z0,t,eps,itmax,tol,tau,verb,f,g,Z,P,niter,errP,errG,obj)
    
    use para
    implicit none

    ! inputs
    integer :: n, m, d, itmax, verb
    real(wp) :: mu(n), nu(m), X(n,d), Y(m,d), f0(n), Z0(m,d)
    real(wp) :: t, eps, tol, tau 
    !f2py intent(in) :: n, m, d, itmax
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
    

    f=f0
    Z=Z0
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
        Z=Z-tau*grad
        errG(i) = sqrt(sum(grad**2))

        obj(i)=-sum(f*mu)-sum(g*nu)+eps*sum(P)+sum(nu*sum((Z-Y)**2,2))/(2*(t-1))
        
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
