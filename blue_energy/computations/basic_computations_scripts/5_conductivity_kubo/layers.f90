program layers

implicit none

integer :: i,j,k,l
integer :: nsteps,nspecies,nspecwork,nlayers,ncount,ncounttot
integer, allocatable, dimension(:) :: nat
integer, allocatable, dimension(:,:) :: zbin

double precision :: du,z
double precision, allocatable, dimension(:) :: zmin,zmax
double precision, allocatable, dimension(:,:) :: xdisp,ydisp,zdisp

character*500 :: dispfile,posfile,outfile1,outfile2, nome_file_input

DO i = 1, iargc()
  CALL getarg(i,nome_file_input)
END DO

open(10,file=nome_file_input)

read(10,*)nsteps
read(10,*)nspecies
read(10,*)nspecwork
allocate(nat(nspecies))
do i=1,nspecies
   read(10,*)nat(i)
   if(i.eq.nspecwork)then
      allocate(zbin(nat(i),nsteps))
      allocate(xdisp(nat(i),nsteps),ydisp(nat(i),nsteps),zdisp(nat(i),nsteps))
   endif
enddo
!print *, 'nat', nat
read(10,*)nlayers
allocate(zmin(nlayers),zmax(nlayers))
do i=1,nlayers
   read(10,*)zmin(i),zmax(i)
enddo
read(10,'(a)')posfile
open(11,file=posfile)
read(10,'(a)')dispfile
open(12,file=dispfile)
read(10,'(a)')outfile1
open(21,file=outfile1)
read(10,'(a)')outfile2
open(22,file=outfile2)!keydisp
close(10)

!print*, 'nsteps', nsteps
do i=1,nsteps
   do j=1,nspecies
      if(j.ne.nspecwork)then
         do k=1,nat(j)
            read(11,*)du,du,du
            read(12,*)du,du,du
         enddo
      else
         do k=1,nat(j)
            read(11,*)du,du,z
            do l=1,nlayers
               if((z.gt.zmin(l)).and.(z.le.zmax(l)))then
                  zbin(k,i)=l
               endif
            enddo
            read(12,*)xdisp(k,i),ydisp(k,i),zdisp(k,i)
            if(zbin(k,i).eq.0)then
               write(6,*)'atom',k,'step',i,'z',z
               stop
            endif
         enddo
      endif
   enddo
enddo
ncounttot=0
do i=2,nsteps
   do k=1,nat(nspecwork)
      if((zbin(k,i).ne.zbin(k,i-1)).and.(zbin(k,i).ne.0))then
         ncount=0
         do l=1,i-1
            if(zbin(k,l).ne.0)then
               write(21,*)xdisp(k,l),ydisp(k,l),zdisp(k,l)
               ncount=ncount+1
               ncounttot=ncounttot+1
            endif
         enddo
         write(21,*)
         ncounttot=ncounttot+1
         write(22,*)zbin(k,i-1),ncount
         do l=1,i-1
            zbin(k,l)=0
         enddo
      endif
   enddo
enddo
do k=1,nat(nspecwork)
   ncount=0
   do l=1,nsteps
      if(zbin(k,l).ne.0)then
         write(21,*)xdisp(k,l),ydisp(k,l),zdisp(k,l)
         ncount=ncount+1
         ncounttot=ncounttot+1
      endif
   enddo
   write(21,*)
   ncounttot=ncounttot+1
   write(22,*)zbin(k,nsteps),ncount
!   write(6,*)ncounttot
enddo

end program
