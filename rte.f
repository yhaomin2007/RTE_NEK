c--------------------------------------------------------------------
      subroutine solid_angle_discretization()
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "RTE_DATA"

      ntot = lx1*ly1*lz1*nelt

      nphi = 2
      ntheta = 4

      sn_pi = 3.1415926
   
      dphi = 0.5*sn_pi/dble(nphi)
      dtheta = sn_pi/dble(ntheta)


      do iphi = 1,nphi*4
      sa_phi(iphi) = 0.5*dphi + dphi*dble(iphi-1)
      enddo

      do itheta = 1,ntheta
      sa_theta(itheta) = 0.5*dtheta + dtheta*dble(itheta-1)
      enddo
 

      if (ldim.eq.2) then  ! 2d
      
      do iphi = 1,nphi*4
      sn_x(itheta,iphi)= 2*sin(sa_phi(iphi))*sin(0.5*dphi)
      sn_y(itheta,iphi)= 2*cos(sa_phi(iphi))*sin(0.5*dphi)
      sa_w(iphi,1) = dphi/(2.0*sn_pi)
      enddo

      else  ! 3d

      sa_tot = 0.0
    
      do ithete = 1,ntheta
      do iphi = 1,nphi*4
      sn_x(itheta,iphi)=sin(sa_phi(iphi))*sin(0.5*dphi)*
     & (dtheta-cos(2*sa_theta(itheta))*sin(dtheta))
      sn_y(itheta,iphi)=cos(sa_phi(iphi))*sin(0.5*dphi)*
     & (dtheta-cos(2*sa_theta(itheta))*sin(dtheta))
      sn_z(itheta,iphi)=0.5*dphi*sin(2*sa_theta(itheta))*sin(dtheta)     
      
      sa_w(iphi,itheta) = sin(sa_theta(itheta))*dtheta*dphi
      sa_tot = sa_tot + sa_w(iphi,itheta)
      enddo
      enddo

      ! solid angle weight
      do ithete = 1,ntheta
      do iphi = 1,nphi*4
      sa_w(iphi,itheta)  =  sa_w(iphi,itheta)/sa_tot
      enddo
      enddo

      endif
 
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
     & .or.(cbc(iface,ie,1).ne.'E  ')) then
        call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1

        call getSnormal(sn2,ix,iy,iz,iface,ie)

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
     & .or.(cbc(iface,ie,1).ne.'E  ')) then

        call surface_int(sint2,sarea2,snfnx,ie,iface)  
        sn2(1) = sint2/sarea2
        call surface_int(sint2,sarea2,snfny,ie,iface)  
        sn2(2) = sint2/sarea2
        call surface_int(sint2,sarea2,snfnz,ie,iface)  
        sn2(3) = sint2/sarea2

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

          do ithete = 1,ntheta
          do iphi = 1,nphi*4
          ipscalar = 1+iphi+(ithete-1)*(nphi*4)

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

      end
c--------------------------------------------------------------------
      subroutine incident_flux_on_bc()
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
     & .or.(cbc(iface,ie,1).ne.'E  ')) then

        call facind(i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,iface)
        do k=k0,k1
        do j=j0,j1
        do i=i0,i1

        if (ldim.eq.2) then

          do iphi = 1,nphi*4          
          ipscalar = 1+iphi
          if (cbc(iface,ie,ipscalar+1).eq.'f  ') then
          Ir_ibc(i,j,k,ie) = Ir_ibc(i,j,k,ie) + t(i,j,k,ie,ipscalar)
          endif
          enddo
      
        else

          do ithete = 1,ntheta
          do iphi = 1,nphi*4
          ipscalar = 1+iphi+(ithete-1)*(nphi*4)

          if (cbc(iface,ie,ipscalar+1).eq.'f  ') then
          Ir_ibc(i,j,k,ie) = Ir_ibc(i,j,k,ie) + t(i,j,k,ie,ipscalar)
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
      subroutine RTE_bc(ix,iy,iz,iside,ie,temp)
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include "RTE_DATA"
      
      integer ix,iy,iz,iside,ie


      t_loc = t(ix,iy,iz,ie,1)

      temp =  emmissivity*sigma*t_loc**4.0/sn_pi + Ir_ibc(ix,iy,iz,ie)

      
      end
c--------------------------------------------------------------------
