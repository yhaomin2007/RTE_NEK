      include 'rte/my_gmres.f'
      include 'rte/my_proj.f'
      include 'rte/advdiff.f'

c--------------------------------------------------------------------
      subroutine RTE_setup()
c
c  setup for RTE solver
c
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"

c solid angle discretization
      call RTE_solid_angle_discretization()  

c bc setup for I fields
      call RTE_bc_setup()

      return
      end
c--------------------------------------------------------------------
      subroutine RTE_solve()
c solve RTE equations iteratively
c code borrowed from yuhsiang
c
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"
c 
c 
      integer n, lt, iter, itmp, ifield_bak
      parameter(lt=lx1*ly1*lz1*lelt)
      real ha(lt), hb(lt), hc(lt), rhs(lt), wt(lt), ub(lt), utmp(lt)
      real msk(lx1*ly1*lz1*lelt)

      real u_exact(lt), u_sol(lt)

      real tol, err, glamax

      integer i
      ntot = lx1*ly1*lz1*nelt
      n =  lx1*ly1*lz1*nelt

      niter = 1
      ifield_bak = ifield

      do iter_outer = 1,niter
  
      if(nid.eq.0) write(6,*) 'solving RTEs, iter_outer: ',iter_outer


c calculate incident radiation flux on boundaries
	  call RTE_incident_radiation_flux_on_bc()
c calculate source term (rhs) of RTE
      call RTE_source_term()

      if (ldim.eq.2) then  ! 2d
      
      do iphi = 1,nphi*4
      if(nid.eq.0) write(6,*) 'solving RTEs, iphi: ',iphi


      ipscalar = 1+iphi
      ifield = ipscalar+1

c solve RTE euqation for this angle.
c      call copy(u_exact,t(1,1,1,1,ipscalar),n) ! referenced sol from restart file

      call cfill(ha,sigma_d,n) ! TODO uservp
      call cfill(hb,sigma_t*sa_a(iphi,1),n) 
      call cfill(hc,1.0,n) 
	  
      call cfill(snxf,sn_x(iphi,1),n)
      call cfill(snyf,sn_y(iphi,1),n)
      call cfill(snzf,sn_z(iphi,1),n)

      call set_advdiff_coef(ha,hb,hc,snxf,snyf,snzf) ! send coef into common blocks

      !call rone(msk,n)
      !do ie = 1,nelt
      !do iface = 1,2*ldim
      !  if (cbc(iface,ie,ifield).eq.'t  ') then  ! for all other bcs
      !  call facev(msk,ie,iface,0.0,lx1,ly1,lz1) 
      !  endif 
      !enddo
      !enddo

      ! call col3  (wt,tmask,tmult,n)
      call rzero(wt,n)
      call col3(wt,tmask(1,1,1,1,ipscalar),tmult(1,1,1,1,ipscalar),n)

      call rzero(rhs,n)
c     call copy(rhs,Ir_src(1,1,1,1,iphi,1),n)
      call setqvol(rhs,utmp)
      call col2  (rhs,bm1,n)
      call dssum (rhs,lx1,ly1,lz1)

      ifield = ipscalar+1
      call rzero(ub,n)
      call bcdirsc(ub)
      call adfax(utmp,ub) ! dssum and mask included
      call sub2(rhs,utmp,n) ! split inhomogenuous Dirichlet BC

c      if (iphi.eq.1) then
c      call outpost(ub,tmask(1,1,1,1,ipscalar),
c     & tmult(1,1,1,1,ipscalar),wt,utmp,'ubb')
c      endif

c      ifield = 2 ! Neumann BC, need test
c      call bcneusc(utmp,1)
cc     call col2 (ta,h1,n)
c      call dssum(utmp,lx1,ly1,lz1)
c      call add2(rhs,utmp,n) ! add inhomogeneous Neumann 'f  ' contribution

      !call col2(rhs,tmask,n)
      call col2(rhs,tmask(1,1,1,1,ipscalar),n)

      ! linear solver
      iter = 10000
      tol = -1.e-5
      call rzero(u_sol,n) ! initial gues
c      call copy(u_sol,t(1,1,1,1,ipscalar),n) 

      call my_gmres(u_sol,rhs,wt,iter,tol,.true.)
c      iter = 300
c      call my_proj (u_sol,rhs,wt,iter,tol,.true.)
      call add2(u_sol,ub,n)

      call copy(t(1,1,1,1,ipscalar),u_sol,n) 

c      if (iphi.eq.1) then
c      call outpost(ub,tmask(1,1,1,1,ipscalar),
c     & tmult(1,1,1,1,ipscalar),wt,t(1,1,1,1,ipscalar),'ubb')
c      endif

      enddo

      else  ! 3d

      do itheta = 1,ntheta
      if(nid.eq.0) write(6,*) 'solving RTEs, itheta: ',itheta
      do iphi = 1,nphi*4
      if(nid.eq.0) write(6,*) 'solving RTEs, iphi: ',iphi
      ipscalar = 1+iphi+(itheta-1)*(nphi*4)
      ifield = ipscalar+1

c solve RTE euqation for this angle.
c      call copy(u_exact,t(1,1,1,1,ipscalar),n) ! referenced sol from restart file

      call cfill(ha,sigma_d,n) ! TODO uservp
      call cfill(hb,sigma_t,n) 
      call cfill(hc,1.0,n) 
	  
      call cfill(snxf,sn_x(iphi,itheta),n)
      call cfill(snyf,sn_y(iphi,itheta),n)
      call cfill(snzf,sn_z(iphi,itheta),n)
	  
      call set_advdiff_coef(ha,hb,hc,snxf,snyf,snzf) ! send coef into common blocks

      !call rone(msk,n)
      !do ie = 1,nelt
      !do iface = 1,2*ldim
      !  if (cbc(iface,ie,ifield).eq.'t  ') then  ! for all other bcs
      !  call facev(msk,ie,iface,0.0,lx1,ly1,lz1) 
      !  endif 
      !enddo
      !enddo
      !call col3  (wt,tmask,tmult,n)
      call col3(wt,tmask(1,1,1,1,ipscalar),tmult(1,1,1,1,ipscalar),n)

      call rzero(rhs,n)
c      call copy(rhs,Ir_src(1,1,1,1,iphi,itheta),n)
      call setqvol(rhs,utmp)
      call col2  (rhs,bm1,n)
      call dssum (rhs,lx1,ly1,lz1)

      ifield = ipscalar+1
      call rzero(ub,n)
      call bcdirsc(ub)
c      call adfax(utmp,ub) ! dssum and mask included
c      call sub2(rhs,utmp,n) ! split inhomogenuous Dirichlet BC

c      ifield = 2 ! Neumann BC, need test
c      call bcneusc(utmp,1)
cc     call col2 (ta,h1,n)
c      call dssum(utmp,lx1,ly1,lz1)
c      call add2(rhs,utmp,n) ! add inhomogeneous Neumann 'f  ' contribution

      !call col2(rhs,tmask,n)
      call col2(rhs,tmask(1,1,1,1,ipscalar),n)


      ! linear solver
      iter = 10000
      tol = -1.e-5
      call rzero(u_sol,n) ! initial gues
c      call copy(u_sol,t(1,1,1,1,ipscalar),n)
      call my_gmres(u_sol,rhs,wt,iter,tol,.true.)
c      iter = 300
c      call my_proj (u_sol,rhs,wt,iter,tol,.true.)
      call add2(u_sol,ub,n)

      call copy(t(1,1,1,1,ipscalar),u_sol,n) 

      enddo
      enddo

      endif

      enddo  ! do iter = 1,niter

      ifield = ifield_bak

	  
      return
      end
c--------------------------------------------------------------------
      subroutine RTE_solid_angle_discretization()
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"
      real dphi,dtheta

      ntot = lx1*ly1*lz1*nelt

      ! nphi = 2     ! change in RTE_DATA
      ! ntheta = 4   ! change in RTE_DATA

      emissivity = 0.5
      sigma_t= 1.0

      sn_pi = 3.1415926
      sigma_b = 5.67E-8   ! Boltzmann constant

      sigma_d = 1e-2 ! artificial viscosity
   
      dphi = 0.5*sn_pi/dble(nphi)
      dtheta = sn_pi/dble(ntheta)

      do iphi = 1,nphi*4
      sa_phi(iphi) = 0.5*dphi + dphi*dble(iphi-1)
      enddo

      do itheta = 1,ntheta
      sa_theta(itheta) = 0.5*dtheta + dtheta*dble(itheta-1)
      enddo

      if (ldim.eq.2) then  ! 2d
      
! https://www.openfoam.com/documentation/guides/latest/api/radiativeIntensityRay_8C_source.html
      do iphi = 1,nphi*4
      sn_x(iphi,1)= 2.0*sin(sa_phi(iphi))*sin(0.5*dphi) 
      sn_y(iphi,1)= 2.0*cos(sa_phi(iphi))*sin(0.5*dphi) 
      sn_z(iphi,1) = 0.0
      sa_a(iphi,1) = dphi
      enddo

      else  ! 3d

      sa_tot = 0.0
    
      do itheta = 1,ntheta
      do iphi = 1,nphi*4
      sn_x(iphi,itheta)=sin(sa_phi(iphi))*sin(0.5*dphi)*
     & (dtheta-cos(2.0*sa_theta(itheta))*sin(dtheta))
      sn_y(iphi,itheta)=cos(sa_phi(iphi))*sin(0.5*dphi)*
     & (dtheta-cos(2.0*sa_theta(itheta))*sin(dtheta))
      sn_z(iphi,itheta)=0.5*dphi*sin(2*sa_theta(itheta))*sin(dtheta)     
      
      sa_a(iphi,itheta) = 2.0*sin(sa_theta(itheta))*sin(0.5*dtheta)*dphi
      enddo

      endif
 
      return
      end
c--------------------------------------------------------------------
      subroutine RTE_bc_setup()
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"

      integer ipscalar
      real angle,ssnn

      real sn2(3),sint2,sarea2
      real sn3(3)
c calculate normal vector of each element face

c judge element face normal with each sold angle, 
c if angle less than 90 degree, then this face should be t for this angle.
c if angle more than 90 degree, then this face should be f for this angle

      ntot = lx1*ly1*lz1*nelt

      call rzero(snfnx,ntot)
      call rzero(snfny,ntot)
      call rzero(snfnz,ntot)

      do ie = 1,nelt
      do iface = 1,2*ldim
        if ((cbc(iface,ie,1).ne.'   ')
     & .and.(cbc(iface,ie,1).ne.'E  ')
     & .and.(cbc(iface,ie,1).ne.'P  ')) then  ! for all other bcs

        call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1

        call getSnormal(sn2,i,j,k,iface,ie)

        snfnx(i,j,k,ie) = sn2(1)
        snfny(i,j,k,ie) = sn2(2)
        snfnz(i,j,k,ie) = sn2(3)

        enddo
        enddo
        enddo
        endif
      enddo
      enddo


      do ie = 1,nelt
      do iface = 1,2*ldim

        if ((cbc(iface,ie,1).ne.'   ')
     & .and.(cbc(iface,ie,1).ne.'E  ')
     & .and.(cbc(iface,ie,1).ne.'P  ')) then

        call surface_int(sint2,sarea2,snfnx,ie,iface)  
        sn2(1) = - sint2/sarea2
        call surface_int(sint2,sarea2,snfny,ie,iface)  
        sn2(2) = - sint2/sarea2
        call surface_int(sint2,sarea2,snfnz,ie,iface)  
        sn2(3) = - sint2/sarea2

        efn(1,iface,ie) = sn2(1)   ! element face normal, point to inside
        efn(2,iface,ie) = sn2(2) 
        efn(3,iface,ie) = sn2(3) 

        if (ldim.eq.2) then

          do iphi = 1,nphi*4
          ipscalar = 1+iphi

          sn3(1) = sin(sa_phi(iphi))
          sn3(2) = cos(sa_phi(iphi))

          ssnn = (sn2(1)*sn3(1)+sn2(2)*sn3(2))

            if (ssnn.ge.0.0) then
            cbc(iface,ie,ipscalar+1) = 't  '
            else 
            cbc(iface,ie,ipscalar+1) = 'f  '
            endif
          enddo
      
        else

          do itheta = 1,ntheta
          do iphi = 1,nphi*4
          ipscalar = 1+iphi+(itheta-1)*(nphi*4)

          sn3(1) = sin(sa_theta(itheta))*sin(sa_phi(iphi))
          sn3(2) = sin(sa_theta(itheta))*cos(sa_phi(iphi))
          sn3(3) = cos(sa_theta(itheta))

          ssnn = (sn2(1)*sn3(1)+sn2(2)*sn3(2)+sn2(3)*sn3(3))

            if  (ssnn.ge.0.0) then
            cbc(iface,ie,ipscalar+1) = 't  '
            else 
            cbc(iface,ie,ipscalar+1) = 'f  '
            endif

          enddo
          enddo

        endif
      
        endif
      
      enddo      
      enddo

c     call rzero(snfnx,ntot)
c     call rzero(snfny,ntot)
c     call rzero(snfnz,ntot)
c	  
c     do ie = 1,nelt
c     do iface = 1,2*ldim
c         iphi = 1
c         ipscalar = 1+iphi
c        if (cbc(iface,ie,ipscalar+1).eq.'t  ') then
c		 
c       call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
c       do k=k0,k1
c       do j=j0,j1
c       do i=i0,i1
c
c       snfnx(i,j,k,ie) = 1
c
c       enddo
c       enddo
c       enddo
c
c        endif
c		 
c        if (cbc(iface,ie,ipscalar+1).eq.'f  ') then
c		 
c       call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
c       do k=k0,k1
c       do j=j0,j1
c       do i=i0,i1
c
c       snfnx(i,j,k,ie) = 2
c
c       enddo
c       enddo
c       enddo
c
c        endif
c	 
c     enddo
c     enddo
c
c     call outpost(snfnx,snfnx,snfnx,snfnx,snfnx,'bcf')
c
      return
      end
c--------------------------------------------------------------------
      subroutine RTE_incident_radiation_flux_on_bc()
c
c sum up all incident fglux on surface...
c
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"


      ntot = lx1*ly1*lz1*nelt

      call rzero(Ir_ibc,ntot)

      do ie = 1,nelt
      do iface = 1,2*ldim

        if ((cbc(iface,ie,1).ne.'   ')
     & .and.(cbc(iface,ie,1).ne.'E  ')
     & .and.(cbc(iface,ie,1).ne.'P  ')) then


        call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1

        if (ldim.eq.2) then

          do iphi = 1,nphi*4          
          ipscalar = 1+iphi
          if (cbc(iface,ie,ipscalar+1).eq.'f  ') then
       Ir_ibc(i,j,k,ie) = Ir_ibc(i,j,k,ie)
     & + t(i,j,k,ie,ipscalar)
     & *(sn_x(iphi,1)*efn(1,iface,ie)
     & +sn_y(iphi,1)*efn(2,iface,ie))

          endif
          enddo
      
        else

          do itheta = 1,ntheta
          do iphi = 1,nphi*4
          ipscalar = 1+iphi+(itheta-1)*(nphi*4)

         if (cbc(iface,ie,ipscalar+1).eq.'f  ') then
      Ir_ibc(i,j,k,ie) = Ir_ibc(i,j,k,ie)
     & + t(i,j,k,ie,ipscalar)
     & *(sn_x(iphi,itheta)*efn(1,iface,ie)
     & +sn_y(iphi,itheta)*efn(2,iface,ie)
     & +sn_z(iphi,itheta)*efn(3,iface,ie))

          endif

          enddo
          enddo

        endif

        enddo
        enddo
        enddo

        endif
      
      enddo      
      enddo

      
      end
c--------------------------------------------------------------------
c--------------------------------------------------------------------
      subroutine RTE_bc(ix,iy,iz,ie,ray)
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"
      
      integer ix,iy,iz,ie
      real ray,t_loc

      t_loc = t(ix,iy,iz,ie,1)

      ray =  (emissivity*sigma_b*t_loc**4.0 
     & + (1.0-emissivity)*Ir_ibc(ix,iy,iz,ie))/sn_pi

      return
      end
c--------------------------------------------------------------------
      subroutine RTE_source_term()
c
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "rte/RTE_DATA"

      real dtdx(lx1,ly1,lz1,lelt)
      real dtdy(lx1,ly1,lz1,lelt)
      real dtdz(lx1,ly1,lz1,lelt)

      integer ipscalar

      ntot = lx1*ly1*lz1*nelt

c Ir_src = sn  dot grad ir

       if (ldim.eq.2) then

          do iphi = 1,nphi*4          
           ipscalar = 1+iphi
c           call gradm1(t(1,1,1,1,ipscalar),dtdx,dtdy,dtdz)
           do i = 1,ntot
c      Ir_src(i,1,1,1,iphi,1) =
c     & - sn_x(iphi,1)*dtdx(i,1,1,1)
c     & - sn_y(iphi,1)*dtdy(i,1,1,1) 
c     & - sigma_t*t(i,1,1,1,ipscalar)
c     & + sigma_t*sigma_b*t(i,1,1,1,1)**4.0/sn_pi
	 
      Ir_src(i,1,1,1,iphi,1) =
     & sigma_t*sa_a(iphi,1)*sigma_b*t(i,1,1,1,1)**4.0/sn_pi
           enddo         
          enddo

        else

          do itheta = 1,ntheta
          do iphi = 1,nphi*4
          ipscalar = 1+iphi+(itheta-1)*(nphi*4)

c           call gradm1(t(1,1,1,1,ipscalar),dtdx,dtdy,dtdz)

           do i = 1,ntot
c      Ir_src(i,1,1,1,iphi,itheta) = 
c     & - sn_x(iphi,itheta)*dtdx(i,1,1,1)
c     & - sn_y(iphi,itheta)*dtdy(i,1,1,1) 
c     & - sn_z(iphi,itheta)*dtdz(i,1,1,1)
c     & - sigma_t*t(i,1,1,1,ipscalar)
c     & + sigma_t*sigma_b*t(i,1,1,1,1)**4.0/sn_pi

      Ir_src(i,1,1,1,iphi,itheta) = 
     & sigma_t*sa_a(iphi,itheta)*sigma_b*t(i,1,1,1,1)**4.0/sn_pi
           enddo         

          enddo
          enddo

        endif


      return
      end
c--------------------------------------------------------------------
