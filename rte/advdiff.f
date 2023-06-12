      subroutine set_advdiff_coef(ha,hb,hc,ux,uy,uz)
c     Set coefficient for advection diffusion and store them into common block
      implicit none
      include 'SIZE'
      include 'rte/ADVDIFF'

      real ux(1), uy(1), uz(1), ha(1), hb(1), hc(1)
      real glmin,glmax,tmn,tmx
      integer n, nd

      n = lx1*ly1*lz1*nelt
      nd = lxd*lyd*lzd*nelv

      call copy(cx,ux,n)
      call copy(cy,uy,n)
      if (ldim.eq.3) call copy(cz,uz,n)
      call set_convect_new(cxd,cyd,czd,cx,cy,cz) ! this is rotated into r-s-t

      call copy(coef_a,ha,n)
      call copy(coef_b,hb,n)
      call copy(coef_c,hc,n)

      ! print min/max of coef
      tmn=glmin(coef_a,n)
      tmx=glmax(coef_a,n)
      if (nio.eq.0) write(*,*)'set coef: a',tmn,tmx
      tmn=glmin(coef_b,n)
      tmx=glmax(coef_b,n)
      if (nio.eq.0) write(*,*)'set coef: b',tmn,tmx
      tmn=glmin(coef_c,n)
      tmx=glmax(coef_c,n)
      if (nio.eq.0) write(*,*)'set coef: c',tmn,tmx

      tmn=glmin(cx,n)
      tmx=glmax(cx,n)
      if (nio.eq.0) write(*,*)'set coef: cx',tmn,tmx
      tmn=glmin(cy,n)
      tmx=glmax(cy,n)
      if (nio.eq.0) write(*,*)'set coef: cy',tmn,tmx
      if (ldim.eq.3) then
         tmn=glmin(cz,n)
         tmx=glmax(cz,n)
         if (nio.eq.0) write(*,*)'set coef: cz',tmn,tmx
      endif

      tmn=glmin(cxd,nd)
      tmx=glmax(cxd,nd)
      if (nio.eq.0) write(*,*)'set coef: cxd',tmn,tmx
      tmn=glmin(cyd,nd)
      tmx=glmax(cyd,nd)
      if (nio.eq.0) write(*,*)'set coef: cyd',tmn,tmx
      if (ldim.eq.3) then
         tmn=glmin(czd,nd)
         tmx=glmax(czd,nd)
         if (nio.eq.0) write(*,*)'set coef: czd',tmn,tmx
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine adfax(ax,xx)
c     Apply advection diffusion operator, ax = H xx
c     H u = -div a grad u + c dot grad u + b u
c
c     Input: (see set_advdiff_coef)
c        (cxd, cyd, czd): dealiased flow
c        coef_a, coef_b, coef_c: coefficients
c        xx
c     Output:
c        ax

      implicit none
      include 'SIZE'
      include 'SOLN' ! mask
      include 'rte/ADVDIFF'

      real    ax(*), xx(*), dtmp(lx1*ly1*lz1*lelt)
      logical ifuf, ifcf
      integer n, imesh, isd

      n = lx1*ly1*lz1*nelt

      ifuf = .false.  ! is u already on the fine mesh
      ifcf = .true.   ! is c already on the fine mesh
      call rzero(dtmp,n) 
      call convect_new(dtmp,xx,ifuf,cxd,cyd,czd,ifcf)
      call col2(dtmp,coef_c,n)

      imesh = 2
      isd = 1
      call rzero(ax,n)
      call axhelm (ax,xx,coef_a,coef_b,imesh,isd)

      call add2 (ax,dtmp,n)
      call dssum(ax,lx1,ly1,lz1)
      call col2 (ax,tmask,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine adf_solve ! TODO: nek interface to set BC, rhs, etc
      include 'SIZE'
      include 'TOTAL'

      return
      end
c-----------------------------------------------------------------------
