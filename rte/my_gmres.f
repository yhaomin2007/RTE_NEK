c-----------------------------------------------------------------------
c     copied from hmh_gmres (Nek5000/core/gmres.f)
c     Solve the Advection Diffusion equation by right-preconditioned 
c     GMRES iteration.
      subroutine my_gmres(x,res,wt,iter,tol,ifprint)
c     Input:
c        x:       initial guess
c        res:     rhs
c        wt:      mask .* mult
c        iter:    maxit
c        tol:     >0=abs, <0=rel
c        ifprint: print residual history
c     Output:
c        x:       solution
c        iter:    number of iterations

      implicit none
      include 'SIZE'
      include 'TSTEP' ! istep
      integer lt
      parameter (lt=lx1*ly1*lz1*lelt)

      real x(lt), res(lt), wt(lt)
     $   , x_gmres(lt), r_gmres(lt), w_gmres(lt)
     $   , v_gmres(lt,lgmres), z_gmres(lt,lgmres)
     $   , h_gmres(lgmres,lgmres), gamma_gmres(lgmres+1)
     $   , c_gmres(lgmres), s_gmres(lgmres)

      real rnorm0,rnorm,ratio,tol,tolpss

      real alpha, l, temp, wk1(lgmres), vlsc3, glsc3
      
      integer iconv, iter, maxit
      integer i,j,k,m,n

      real*8 etime1,etime2,etimep,dnekclock
      logical ifprint

      n = lx1*ly1*lz1*nelt ! TODO

      etime1 = dnekclock()
      etimep = 0.
      maxit = iter
      iter  = 0
      m     = lgmres
      tolpss= abs(tol)

      iconv = 0
c      call rzero(x_gmres,n)
      call copy(x_gmres,x,n)

      do while (iconv.eq.0)

         if(iter.eq.0) then 
            call copy(r_gmres,res,n)
         else
            !update residual
            call copy       (r_gmres,res,n)          ! r = res
            call mygmres_ax (w_gmres,x_gmres,n)      ! w = A x
            call sub3       (r_gmres,res,w_gmres,n)  ! r = r - w
         endif
                                                            !            ______
         gamma_gmres(1) = sqrt(glsc3(r_gmres,r_gmres,wt,n)) ! gamma  = \/ (r,r) 
         if(gamma_gmres(1) .eq. 0.) goto 9000     !check for lucky convergence

         if(iter.eq.0) then
            rnorm0 = gamma_gmres(1)
            if (tol.lt.0) tolpss=abs(tol)*gamma_gmres(1)
         endif

         temp = 1./gamma_gmres(1)
         call cmult2(v_gmres(1,1),r_gmres,temp,n) ! v  = r / gamma
                                                  !  1            1
         do j=1,m !  Krylov space
            iter = iter+1

            ! preconditioned
            etime2 = dnekclock()
            call mygmres_prec(z_gmres(1,j), v_gmres(1,j), n) ! z = M^{-1} v
            etimep = etimep + dnekclock()-etime2

            call mygmres_ax (w_gmres,z_gmres(1,j),n)  ! w = A z
                                                    !        j

            ! Gram-Schmidt
            do i=1,j
               h_gmres(i,j)=vlsc3(w_gmres,v_gmres(1,i),wt,n) ! h    = (w,v )
            enddo                                            !  i,j       i
            call gop(h_gmres(1,j),wk1,'+  ',j)               ! sum over P procs
            do i=1,j
               call add2s2(w_gmres,v_gmres(1,i),-h_gmres(i,j),n) ! w = w - h    v
            enddo                                                !          i,j  i

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h_gmres(i,j)
               h_gmres(i  ,j)=  c_gmres(i)*temp
     $                        + s_gmres(i)*h_gmres(i+1,j)
               h_gmres(i+1,j)= -s_gmres(i)*temp
     $                        + c_gmres(i)*h_gmres(i+1,j)
            enddo
                                                      !            ______
            alpha = sqrt(glsc3(w_gmres,w_gmres,wt,n)) ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged

            l = sqrt(h_gmres(j,j)*h_gmres(j,j)+alpha*alpha)
            temp = 1./l
            c_gmres(j) = h_gmres(j,j) * temp
            s_gmres(j) = alpha  * temp
            h_gmres(j,j) = l

            gamma_gmres(j+1) = -s_gmres(j) * gamma_gmres(j)
            gamma_gmres(j)   =  c_gmres(j) * gamma_gmres(j)

            rnorm = abs(gamma_gmres(j+1))
            ratio = rnorm/rnorm0
            if (ifprint.and.nio.eq.0)
     $         write (6,66) iter,tolpss,rnorm,rnorm0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')

            if (iter+1.gt.maxit) goto 900
            if (rnorm .lt.tolpss) goto 900  !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v_gmres(1,j+1),w_gmres,temp,n) ! v    = w / alpha
                                                       !  j+1            
         enddo
  900    iconv = 1
 1000    continue

         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma_gmres(k)
            do i=j,k+1,-1
               temp = temp - h_gmres(k,i)*c_gmres(i)
            enddo
            c_gmres(k) = temp/h_gmres(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x_gmres,z_gmres(1,i),c_gmres(i),n) ! x = x + c  z
         enddo                                             !          i  i
      enddo
 9000 continue

      call copy(x,x_gmres,n)

      etime1 = dnekclock()-etime1
      if (nio.eq.0) write(6,8) istep,iter,rnorm,rnorm0,tolpss,etimep,
     &                            etime1
    8 format(4x,i7,'  my_gmres ',4x,i5,1p5e13.4)

      return
      end
c-----------------------------------------------------------------------
c     ax = A x 
      subroutine mygmres_ax(ax,x,n)
      implicit none
      real ax(1), x(1)
      integer n

      call adfax(ax,x)

      return
      end
c-----------------------------------------------------------------------
c     z = M^{-1} v, where M might depends on r (Schwarz)
      subroutine mygmres_prec(z,v,n)
      implicit none
      real z(1), r(1), v(1)
      integer n

      call copy(z,v,n) ! no-prec

      return
      end
c-----------------------------------------------------------------------

