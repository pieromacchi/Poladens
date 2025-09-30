        Module variables
        Real*8, allocatable:: poldens(:,:,:),dens(:,:),xyz(:,:)
		Real*8, allocatable:: totpol(:),totatpol(:)
		Real*8, allocatable:: rzeta(:),xyzatom(:,:),atpol(:,:,:)
		Real*8, allocatable:: atpro(:),atch(:,:),atdip(:,:,:)
		Real*8, allocatable:: wh(:),atpolint(:,:,:)
		Real*8, allocatable:: atvol(:)
		Integer*4, allocatable:: iZ(:)
        Real*8:: iamden(500,92)
        Real*8::r(3),step(3,3),xyz0(3),efield(3,-3:3),field,molpol(3,3)
		Real*8:: Totel,Totdens,Stepvol,MolVol,DipMol(3,-3:3)
        Integer*4::i,j,k,l,m,n,npoints,natom,ndim,nstep(3)
        Character*120:: line,cubename(-3:3),alphaname(3,3)
        Character*120:: totalphaname,totatalphaname
		logical::lcalcpoldens,lhirsh,lvor

        End Module Variables


		program Poladens
		Use Variables

	
		field=0.005
		lcalcpoldens=.false.
        lhirsh=.false.
        lvor=.false.
		Open (3,file='poladens.inp',form='formatted')
		Open (2,file='poladens.out',form='formatted')
        call readinput
		Print *, ' Input read'
		call setvariables
		Print *, 'variables set'
        call setproden
        Print *, 'Proden set'
		call readcubefiles
		Print *, 'cubefiles read'
		call calcpolarizability
        Print *, 'Polarizabilities calculated'
        if (lcalcpoldens) then 	
	 	 call writealphafiles
 	 	 Print *, 'cube files written'
		end if
		call writeoutput
		Print *, 'Output files written'
		close(2)
		End


		Subroutine SetVariables
		Use Variables
	
		molpol=0.0d0
		cubename=' '
		cubename(0)='zero.cube'
		cubename(1)='xp.cube'
		cubename(-1)='xm.cube'
		cubename(2)='yp.cube'
		cubename(-2)='ym.cube'
		cubename(3)='zp.cube'
		cubename(-3)='zm.cube'
		alphaname='alphaij.cube'
		do i=1,3
		do j=1,3
		 write(alphaname(i,j)(6:7),'(2i1)') i,j
		end do
		end do
		totalphaname='alpha.cube'
        totatalphaname='intrinsicalpha.cube'
		efield=0.0d0
		efield(1,1)=-field
		efield(1,-1)=field
		efield(2,2)=-field
		efield(2,-2)=field
		efield(3,3)=-field
		efield(3,-3)=field

		Return
		End Subroutine SetVariables

		Subroutine AllocateVariables
		Use Variables
		allocate(xyzatom(3,natom),iZ(natom),rzeta(natom))
        allocate(poldens(3,3,npoints),dens(npoints,-3:3)
     . ,xyz(3,npoints))
        allocate(atvol(natom))
		allocate(totpol(npoints),totatpol(npoints))
		allocate(atpol(3,3,natom),atpro(natom),atch(natom,-3:3),
     . atdip(3,natom,-3:3),atpolint(3,3,natom),wh(natom))
		totpol=0.0d0
        totatpol=0.0d0
		atpol=0.0d0
		atpolint=0.0d0
		atch=0.0d0
		atdip=0.0d0
        poldens=0.0d0
        dens=0.0d0
        xyz=0.0d0
		Return
		End Subroutine AllocateVariables

		Subroutine readcubefiles
		Use Variables
	
		open(10,file=cubename(0),form='formatted')
		read(10,'(a)') line
		read(10,'(a)') line
		read(10,*) natom, xyz0(1:3)	
		do i=1,3
		 read(10,*) nstep(i),step(i,1:3)
		end do
		stepvol=step(1,1)*step(2,2)*step(3,3)
		Print *, 'grid point volume', stepvol
        npoints=nstep(1)*nstep(2)*nstep(3)
        Print *, 'Number of steps', npoints
        
        call allocatevariables
		do i=1,natom
		 read(10,*) iz(i),rzeta(i),xyzatom(1:3,i)
		end do
		read(10,*) dens(1:npoints,0)
		close(10)

		do j=-3,3
         if(j==0) cycle
	 	 open(j+10,file=cubename(j),form='formatted')
	 	 do i=1,natom+6
	 	  read(j+10,'(a)') line
	 	 end do
	 	 read(j+10,*) dens(1:npoints,j)
	 	 close(j+10)
		end do

		Return
		End Subroutine readcubefiles

		Subroutine calcpolarizability
		Use Variables
        real*8:: ra(3)
		integer*4:: j1,j2,j3,jpoint

		TotEl=0.0d0
		Totdens=0.0d0
		jpoint=0
		Molvol=0.0d0
		atvol=0.0d0
		DipMol=0.0d0
		do j1=1,nstep(1)
	 	 do j2=1,nstep(2)
		  do j3=1,nstep(3)
		    jpoint=jpoint+1	
		    xyz(1,jpoint)=xyz0(1)+(j1-1)*step(1,1)+(j2-1)*step(1,2)
     . +(j3-1)*step(1,3)
		    xyz(2,jpoint)=xyz0(2)+(j1-1)*step(2,1)+(j2-1)*step(2,2)
     . +(j3-1)*step(2,3)
		    xyz(3,jpoint)=xyz0(3)+(j1-1)*step(3,1)+(j2-1)*step(3,2)
     . +(j3-1)*step(3,3)
	    r(1:3)=xyz(1:3,jpoint)
c        write(22,'(3f12.6)') r
c	    if ((dens(jpoint,0) < 1.0e-6)) cycle
	    Totel=totel+dens(jpoint,0)*stepvol
	    Totdens=totdens+dens(jpoint,0)
	    MolVol=MolVol+stepvol
	    do i=-3,3
             dipMol(1:3,i)=
     . dipMol(1:3,i)+dens(jpoint,i)*r(1:3)*stepvol
	    end do
	    call weight
	    do k=1,natom
             ra(1:3)=xyzatom(1:3,k)
             do i=-3,3
	      atch(k,i)=atch(k,i)+dens(jpoint,i)*stepvol*wh(k)
	      atdip(1:3,k,i)=atdip(1:3,k,i)+
     . wh(k)*dens(jpoint,i)*(r(1:3)-ra(1:3))
             end do
	    end do
            do j=1,3
	     do i=1,3
              poldens(i,j,jpoint)=
     . r(i)*(dens(jpoint,j)-dens(jpoint,-j))/(2*field)
c   	      if(r(1)==r(2)) write(22,*)
c     . jpoint,i,j,r,dens(jpoint,j),dens(jpoint,-j)
		      molpol(i,j)=molpol(i,j)+poldens(i,j,jpoint)
		      do k=1,natom
		       atpol(i,j,k)=
     . atpol(i,j,k)+(r(i)-xyzatom(i,k))*wh(k)
     . *(dens(jpoint,j)-dens(jpoint,-j))/(2*field)
               if(i==j) totatpol(jpoint)=totatpol(jpoint)+
     . (r(i)-xyzatom(i,k))*wh(k)
     . *(dens(jpoint,j)-dens(jpoint,-j))/(2*field)
		      end do
		      if(i==j) totpol(jpoint)=totpol(jpoint)+poldens(i,j,jpoint)
		     end do
		    end do
		   end do
		 end do
		end do
        atdip=atdip*stepvol
        atpol=atpol*stepvol
        molpol=molpol*stepvol
		totpol=totpol/3
        totatpol=totatpol/3
		do k=1,natom
         do j=1,3
          do i=1,3
           atpolint(i,j,k)=
     . (atdip(i,k,j)-atdip(i,k,-j))/(2*field)
		  end do
		 end do
		end do
		Return
		End Subroutine calcpolarizability

		Subroutine writealphafiles
		Use variables

		do i=1,3
		 do j=1,3
		  open(1,file=alphaname(i,j),form='formatted')
          Print '(a,a)', 'Writing ', alphaname(i,j) 
		  write(1,'(a)') alphaname(i,j)(1:7)
		  write(1,'(a)') 'Polarizability density'
		  write(1,'(i5,3f12.6)') natom,xyz0
		  do k=1,3
		   write(1,'(i5,3f12.6)') nstep(k),step(k,1:3)
		  end do
		  do k=1,natom
		   write(1,'(i5,4f12.6)') iZ(k),rzeta(k),xyzatom(1:3,k)
		  end do
		  do k=1,nstep(1)
		   do l=1,nstep(2)
		    jpoint=(k-1)*nstep(2)*nstep(3)+(l-1)*nstep(3)
		    write(1,'(6e13.5)') poldens(i,j,jpoint+1:jpoint+nstep(3))
		   end do
		  end do
		  close(1)
		 end do
		end do
		open(1,file=totalphaname,form='formatted')
		write(1,'(a)') 'Total Polarizability'
        write(1,'(a)') 'Polarizability density'
        write(1,'(i5,3f12.6)') natom,xyz0
        do k=1,3
          write(1,'(i5,3f12.6)') nstep(k),step(k,1:3)
        end do
        do k=1,natom
         write(1,'(i5,4f12.6)') iZ(k),rzeta(k),xyzatom(1:3,k)
        end do
        do k=1,nstep(1)
         do l=1,nstep(2)
          jpoint=(k-1)*nstep(2)*nstep(3)+(l-1)*nstep(3)
          write(1,'(6e13.5)') totpol(jpoint+1:jpoint+nstep(3))
         end do
        end do
        close(1)
        open(1,file=totatalphaname,form='formatted')
        write(1,'(a)') 'Intrinsic Polarizability'
        write(1,'(a)') 'Polarizability density'
        write(1,'(i5,3f12.6)') natom,xyz0
        do k=1,3
          write(1,'(i5,3f12.6)') nstep(k),step(k,1:3)
        end do
        do k=1,natom
         write(1,'(i5,4f12.6)') iZ(k),rzeta(k),xyzatom(1:3,k)
        end do
        do k=1,nstep(1)
         do l=1,nstep(2)
          jpoint=(k-1)*nstep(2)*nstep(3)+(l-1)*nstep(3)
          write(1,'(6e13.5)') totatpol(jpoint+1:jpoint+nstep(3))
         end do
        end do
        close(1)

        Print *, 'Alpha.cube files written'

		Return
		End subroutine writealphafiles

		Subroutine Readinput
		Use Variables
		Integer*4::ierr	

		ierr=0
		do 
		 read(3,'(a)',iostat=ierr)  line
		 if (ierr/=0) exit
   	      if(line(1:1)=='!') cycle
		 j1=index(line,'field')
		 if(j1>0) read(line(j1+5:),*) field
		 j1=index(line,'write poldens')
		 if(j1>0) lcalcpoldens=.true.
  	     j1=index(line,'hirsh')
         if(j1>0) lhirsh=.true.
         j1=index(line,'voronoi')
         if(j1>0) lvor=.true.
         j1=index(line,'end')
         if(j1>0) return
		end do

		Return
		End Subroutine Readinput

		Subroutine Writeoutput
        Use Variables
		Real*8:: Tot(4),TotInt(3,3)	

        Write(2,'(a,f12.4)')
     . 'TOTAL Molecular Volume', MolVol	
        Write(2,'(a)') 'FIELD AND DIPOLE MOMENTS'
		do i=-3,3
		 Write(2,'(6f12.6)') efield(1:3,i),dipmol(1:3,i)*2.541746
		end do
		Write(2,'(a,f12.4)') 
     . 'TOTAL NUMBER OF INTEGRATED ELECTRONS', Totel
		Write(2,'(a,f12.4)') 
     . 'TOTAL DENS', Totdens
     
		Write(2,'(a)') ' ATOMIC CHARGES, ATOMIC DIPOLES and VOLUMES'
		do j=-3,3
		 tot=0.0d0
         Write (2,*)
         write (2,'(a,3f12.6)') 'FIELD: ',efield(1:3,j)
		 do i=1,natom
		  Write(2,'(i4,f12.4,4f12.6)') iZ(i),atch(i,j),atdip(1:3,i,j),
     .atvol(i)
		  tot(1)=tot(1)+atch(i,j)
		  tot(2:4)=tot(2:4)+atdip(1:3,i,j)
 		 end do
		 Write(2,'(a4,f12.4,3f12.6)') 'Tot ',tot
		end do
		Write (2,'(a)') 'MOLECULAR POLARIZABILITIES'
		write(2,*)
 		Write(2,'(3f12.5)') molpol
        Write(2,'(a)') ' INTRINSIC ATOMIC POLARIABILITIES'
        totint=0.0d0
        do i=1,natom
         Write(2,'(i4)') iZ(i)
          totint=totint+atpol(1:3,1:3,i)
		  do j=1,3
		   Write(2,'(6f12.6)') atpol(j,1:3,i), atpolint(j,1:3,i)
		  end do
        end do 
        write(2,*)    
        Write(2,'(a)') ' INTRINSIC ATOMIC POLARIABILITIES'
        do j=1,3
         Write(2,'(3f12.6)') totint(j,1:3)
        end do
        Return
        End Subroutine Writeoutput

		Subroutine weight
		Use Variables
		Real*8 :: proden,rdist,ra(3),rd(3),rav(3),norm
c Hirshfeld ---> lhirsh
c Z/r ----> default
c Voronoi ---> lVor

		proden=0.0d0

		do i=1,natom
         norm=(rzeta(i)**4)/6.0
         atpro(i)=0.0d0
         ra(1:3)=xyzatom(1:3,i)
         rd=ra-r
		 rdist=sqrt(rd(1)**2+rd(2)**2+rd(3)**2)
		 if(rdist<1.0d-10) then 
		  rdist = 1.0d-10
		 end if
         if(lhirsh) then 
          j=int(rzeta(i))
          jr=int(100*rdist)
          if(jr<=500) atpro(i)=iamden(jr,j)
         elseif(lvor) then
          wh(i)=1.0
          do j=1,natom
           if(i==j) cycle
           rav(1:3)=xyzatom(1:3,j)
           rd=rav-r
           rvor=sqrt(rd(1)**2+rd(2)**2+rd(3)**2)
           if(rvor<rdist) wh(i)=0.0d0
           if(rvor==rdist.and.wh(i).gt.0) wh(i)=0.5
          end do
          atvol(i)=atvol(i)+stepvol*wh(i)
         else
c          atpro(i)=rzeta(i)/rdist
          atpro(i)=norm*exp(-rzeta(i)*rdist/0.529177)
         end if
         if(lvor) cycle
         proden=proden+atpro(i)
		end do
        if(lvor) return
        if (proden < 1.0d-8) then
         proden=1.0d-8
        end if
 		wh(1:natom)=atpro(1:natom)/proden
		Return
		End    

        Subroutine SetProden
        Use Variables
        character*120 poladenspath
        integer*4 jzeta
c specify the path here
		poladenspath='./iam.dat'
        open(22,file=poladenspath,form='formatted')

        do
         read(22,'(a)',iostat=ierr)  line
         if (ierr/=0) exit
         if (index(line,'DATAFILE')>0) then 
          read(line(9:),*) jzeta
          read(22,*)
          read(22,*)
          read(22,*)
          do j=1,500    
           read(22,*) iamden(j,jzeta)
          end do
         end if
        end do
        
        Return 
        End Subroutine setproden        



