c-----------------------------------------------------------------------
c     copied from my MATLAB code 
c     Solve the Advection Diffusion equation by Krylov subspace projections
      subroutine my_proj(x,rhs,wt,iter,tol,ifprint)
c     Input:
c        u:       initial guess
c        r:       rhs
c        wt:      mask .* mult
c        iter:    maxit
c        tol:     >0=abs, <0=rel
c        ifprint: print residual history
c     Output:
c        u:       solution 
c        iter:    number of iterations
      implicit none
      include 'SIZE'
      include 'TSTEP' ! istep

      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real x(lt), rhs(lt), wt(lt)
     $   , u(lt), r(lt), z(lt), w(lt), p(lt, lgmres), q(lt,lgmres)

      integer iter, maxit, j, n
      real tol, res0, res, rtol, ratio, alpha, beta, betai
      real glsc3, glamax

      real*8 etime1,etime2,etimep,dnekclock
      logical ifprint

      n = lx1*ly1*lz1*nelt

      etime1 = dnekclock()
      etimep = 0.
      maxit = iter
      rtol = abs(tol)
      if (maxit.gt.lgmres) 
     $   call exitt('my_proj needs lgmres>=maxit$',maxit)

      iter = 0
      call copy(u,x,n)
      call copy(r,rhs,n)
      if (glamax(u,n).gt.0) then
         call mypj_ax(w,u,n)     ! w = A z
         call sub2(r,w,n)        ! r = r - A u
      endif

      res = glsc3(r,r,wt,n)
      if (res.gt.0) res = sqrt(res)
      if (res.eq.0) goto 900     ! lucky conv.
      res0 = res
      if (tol.lt.0) rtol = abs(tol)*res0

      iter = 1
      do while (res.gt.rtol.AND.iter.le.maxit)

         etime2 = dnekclock()
         call mypj_prec(z,r,n)   ! z = M^{-1} r
         etimep = etimep + dnekclock()-etime2

         call mypj_ax(w,z,n)     ! w = A z

         do j=1,iter-1
            beta = glsc3(q(1,j), w, wt, n)
            call add2s2(z, p(1,j), -beta, n) ! z = z - beta * P(:,j)
            call add2s2(w, q(1,j), -beta, n) ! w = w - beta * Q(:,j)
         enddo

         beta = glsc3(w, w, wt, n)
         if (beta.gt.0) beta = sqrt(beta)
         if (beta.lt.1.1102e-12) goto 900 ! 5000*eps ~1.1102e-12, leftover iter is null space

         betai = 1./beta 
         call copy (p(1,iter), z, n)
         call cmult(p(1,iter), betai, n)  ! P(:,iter) = z/beta
         call copy (q(1,iter), w, n)
         call cmult(q(1,iter), betai, n)  ! Q(:,iter) = w/beta

         alpha = glsc3(Q(1,iter), r, wt, n)

         call add2s2(u, p(1,iter),  alpha, n) ! u = u + alpha * P(:,iter)
         call add2s2(r, q(1,iter), -alpha, n) ! r = r - alpha * Q(:,iter)

         res = glsc3(r,r,wt,n)
         if (res.gt.0) res = sqrt(res)
         ratio = res / res0
         if (ifprint.AND.nio.eq.0) 
     $      write(6,66) iter,rtol,res,res0,ratio,istep 
   66       format(i5,1p4e12.5,i8,' proj div')

         iter = iter + 1
      enddo
      iter = iter - 1
  900 continue

      call copy(x,u,n)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,8) istep,iter,res,res0,rtol,etimep,etime1
    8 format(4x,i7,'  my_gmres ',4x,i5,1p5e13.4)

      return
      end
c-----------------------------------------------------------------------
c     ax = A x 
      subroutine mypj_ax(ax,x,n)
      implicit none
      real ax(1), x(1)
      integer n
      
      call adfax(ax,x)
      
      return
      end
c-----------------------------------------------------------------------
c     z = M^{-1} v, where M might depends on r (Schwarz)
      subroutine mypj_prec(z,r,n)
      implicit none
      real z(1), r(1), v(1)
      integer n
            
      call copy(z,r,n) ! no-prec
            
      return
      end   
c-----------------------------------------------------------------------
   
      
