c--------------------------------------------------------------------
      subroutine RTE_solid_angle_discretization()
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "RTE_DATA"

      real dphi,dtheta

      ntot = lx1*ly1*lz1*nelt

      ! nphi = 2     ! change in RTE_DATA
      ! ntheta = 4   ! change in RTE_DATA

      emissivity = 0.5
      sigma_t= 1.0

      sn_pi = 3.1415926
      sigma_b = 5.67E-8   ! Boltzmann constant
   
      dphi = 0.5*sn_pi/dble(nphi)
      dtheta = sn_pi/dble(ntheta)

      do iphi = 1,nphi*4
      sa_phi(iphi) = 0.5*dphi + dphi*dble(iphi-1)
      enddo

      do itheta = 1,ntheta
      sa_theta(itheta) = 0.5*dtheta + dtheta*dble(itheta-1)
      enddo
 

      if (ldim.eq.2) then  ! 2d
      itheta = 1
      do iphi = 1,nphi*4
      sn_x(iphi,itheta)= 2.0*sin(sa_phi(iphi))*sin(0.5*dphi)
      sn_y(iphi,itheta)= 2.0*cos(sa_phi(iphi))*sin(0.5*dphi)
      sa_w(iphi,1) = dphi/(2.0*sn_pi)
      enddo

      else  ! 3d

      sa_tot = 0.0
    
      do itheta = 1,ntheta
      do iphi = 1,nphi*4
      sn_x(iphi,itheta)=sin(sa_phi(iphi))*sin(0.5*dphi)*
     & (dtheta-cos(2.0*sa_theta(itheta))*sin(dtheta))
      sn_y(iphi,itheta)=cos(sa_phi(iphi))*sin(0.5*dphi)*
     & (dtheta-cos(2.0*sa_theta(itheta))*sin(dtheta))
      sn_z(iphi,itheta)=0.5*dphi*sin(2.0*sa_theta(itheta))*sin(dtheta)     
      
      sa_w(iphi,itheta) = sin(sa_theta(itheta))*dtheta*dphi
      sa_tot = sa_tot + sa_w(iphi,itheta)
      enddo
      enddo

      ! solid angle weight
      do itheta = 1,ntheta
      do iphi = 1,nphi*4
      sa_w(iphi,itheta)  =  sa_w(iphi,itheta)/sa_tot
      enddo
      enddo

      endif
 
      return
      end
c--------------------------------------------------------------------
      subroutine RTE_bc_setup()
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "RTE_DATA"

      integer ipscalar
      real angle

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
     & .or.(cbc(iface,ie,1).ne.'E  ')
     & .or.(cbc(iface,ie,1).ne.'P  ')) then

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
     & .or.(cbc(iface,ie,1).ne.'E  ')
     & .or.(cbc(iface,ie,1).ne.'P  ')) then

        call surface_int(sint2,sarea2,snfnx,ie,iface)  
        sn2(1) = - sint2/sarea2
        call surface_int(sint2,sarea2,snfny,ie,iface)  
        sn2(2) = - sint2/sarea2
        call surface_int(sint2,sarea2,snfnz,ie,iface)  
        sn2(3) = - sint2/sarea2

        efn(1,iface,ie) =   sn2(1)   ! element face normal, point to inside
        efn(2,iface,ie) =   sn2(2) 
        efn(3,iface,ie) =   sn2(3) 

        if (ldim.eq.2) then

          do iphi = 1,nphi*4
          ipscalar = 1+iphi

          sn3(1) = cos(sa_phi(iphi))
          sn3(2) = sin(sa_phi(iphi))

          angle = acos(sn2(1)*sn3(1)+sn2(2)*sn3(2))

            if (angle.le.(sn_pi/2.0)) then
            cbc(iface,ie,ipscalar+1) = 't  '
            else 
            cbc(iface,ie,ipscalar+1) = 'f  '
            endif
          enddo
      
        else

          do itheta = 1,ntheta
          do iphi = 1,nphi*4
          ipscalar = 1+iphi+(itheta-1)*(nphi*4)

          sn3(1) = cos(sa_phi(iphi))*sin(sa_theta(itheta))
          sn3(2) = sin(sa_phi(iphi))*sin(sa_theta(itheta))
          sn3(3) = cos(sa_theta(itheta))

          angle = acos(sn2(1)*sn3(1)+sn2(2)*sn3(2)+sn2(3)*sn3(3))
          
            if (angle.le.(sn_pi/2.0)) then
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
      include "RTE_DATA"


      ntot = lx1*ly1*lz1*nelt

      call rzero(Ir_ibc,ntot)

      do ie = 1,nelt
      do iface = 1,2*ldim

        if ((cbc(iface,ie,1).ne.'   ')
     & .or.(cbc(iface,ie,1).ne.'E  ')
     & .or.(cbc(iface,ie,1).ne.'P  ')) then


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

      varmin = glmin(Ir_ibc,ntot)
      varmax = glmax(Ir_ibc,ntot)
      if (nid.eq.0) write(6,*) 'Ir_ibc min/max:',varmin,'-',varmax   
   
      end
c--------------------------------------------------------------------
c--------------------------------------------------------------------
      subroutine RTE_bc(ix,iy,iz,iside,ie,temp)
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "RTE_DATA"
      
      integer ix,iy,iz,iside,ie
      real temp,t_loc

      t_loc = t(ix,iy,iz,ie,1)

      temp =  (emissivity*sigma_b*t_loc**4.0 
     & + (1.0-emissivity)*Ir_ibc(ix,iy,iz,ie))/sn_pi

      return
      end
c--------------------------------------------------------------------
      subroutine RTE_source_term()
c
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "RTE_DATA"

      real dtdx(lx1,ly1,lz1,lelt)
      real dtdy(lx1,ly1,lz1,lelt)
      real dtdz(lx1,ly1,lz1,lelt)

      integer ipscalar

      ntot = lx1*ly1*lz1*nelt

c Ir_src = sn  dot grad ir

       if (ldim.eq.2) then

          do iphi = 1,nphi*4          
           ipscalar = 1+iphi

      call gradm1(dtdx,dtdy,dtdz,t(1,1,1,1,ipscalar))

      !call opgrad  (dtdx,dtdy,dtdz,t(1,1,1,1,ipscalar))

      ! ensure continuity at element boundaries
      call opcolv (dtdx,dtdy,dtdz,bm1)
      call opdssum (dtdx,dtdy,dtdz)
      call opcolv  (dtdx,dtdy,dtdz,binvm1)


           do i = 1,ntot
      Ir_src(i,1,1,1,iphi,1) =
     & - sn_x(iphi,1)*dtdx(i,1,1,1)
     & - sn_y(iphi,1)*dtdy(i,1,1,1) 
     & - sigma_t*t(i,1,1,1,ipscalar)
     & + sigma_t*sigma_b*t(i,1,1,1,1)**4.0/sn_pi
           enddo         
          enddo

        else

          do itheta = 1,ntheta
          do iphi = 1,nphi*4
           ipscalar = 1+iphi+(itheta-1)*(nphi*4)
      call gradm1(dtdx,dtdy,dtdz,t(1,1,1,1,ipscalar))

	  !call opgrad  (dtdx,dtdy,dtdz,t(1,1,1,1,ipscalar))

      ! ensure continuity at element boundaries
      call opcolv (dtdx,dtdy,dtdz,bm1)
      call opdssum (dtdx,dtdy,dtdz)
      call opcolv  (dtdx,dtdy,dtdz,binvm1)
		   
           do i = 1,ntot
      Ir_src(i,1,1,1,iphi,itheta) = 
     & - sn_x(iphi,itheta)*dtdx(i,1,1,1)
     & - sn_y(iphi,itheta)*dtdy(i,1,1,1) 
     & - sn_z(iphi,itheta)*dtdz(i,1,1,1)
     & - sigma_t*t(i,1,1,1,ipscalar)
     & + sigma_t*sigma_b*t(i,1,1,1,1)**4.0/sn_pi

           enddo         

          enddo
          enddo

        endif


       if (ldim.eq.2) then

          do iphi = 1,nphi*4          
          varmin = glmin(Ir_src(1,1,1,1,iphi,1),ntot)
          varmax = glmax(Ir_src(1,1,1,1,iphi,1),ntot)
          if (nid.eq.0) then
       write(6,*) 'Ir_src iphi ',iphi,' min/max:',varmin,'-',varmax
          endif
          enddo

        else

          do itheta = 1,ntheta
          do iphi = 1,nphi*4
          varmin = glmin(Ir_src(1,1,1,1,iphi,itheta),ntot)
          varmax = glmax(Ir_src(1,1,1,1,iphi,itheta),ntot)
          if (nid.eq.0) then
       write(6,*) 'Ir_src iphi ',iphi,' itheta ',itheta,
     & ' min/max:',varmin,'-',varmax
          endif

          enddo
          enddo

        endif

      return
      end
c--------------------------------------------------------------------
