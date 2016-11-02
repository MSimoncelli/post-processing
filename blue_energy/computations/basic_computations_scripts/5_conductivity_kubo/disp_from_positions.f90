program dispfrompositions

implicit none

integer :: i,j,k,l,icount,nspecies,nposfiles,natmax
integer, allocatable, dimension(:) :: nat,nstep

double precision :: lbox,halfbox,dx,dy,dz
double precision, allocatable, dimension(:,:) :: dxtot,dytot,dztot
double precision, allocatable, dimension(:,:) :: xold,xnew
double precision, allocatable, dimension(:,:) :: yold,ynew
double precision, allocatable, dimension(:,:) :: zold,znew

logical, allocatable, dimension(:) :: extractdisp

character*500 ::outfile,outposfile,nome_file_input,out_all_disp_file
character*200, allocatable, dimension(:) :: posfile

natmax=0

DO i = 1, iargc()
  CALL getarg(i,nome_file_input)
END DO



open(10,file=nome_file_input)
read(10,*)nspecies
allocate(nat(nspecies),extractdisp(nspecies))
do i=1,nspecies
   read(10,*)nat(i),extractdisp(i)
   if(nat(i).gt.natmax)natmax=nat(i)
enddo
allocate(xnew(nspecies,natmax),ynew(nspecies,natmax),znew(nspecies,natmax))
allocate(xold(nspecies,natmax),yold(nspecies,natmax),zold(nspecies,natmax))
allocate(dxtot(nspecies,natmax),dytot(nspecies,natmax),dztot(nspecies,natmax))
read(10,*)nposfiles
allocate(posfile(nposfiles),nstep(nposfiles))
do i=1,nposfiles
   read(10,'(a)')posfile(i)
   read(10,*)nstep(i)
enddo
read(10,*)lbox
halfbox=lbox/2.0d0
read(10,"(A)") outfile   
read(10,"(A)") outposfile   
read(10,"(A)") out_all_disp_file
close(10)

dxtot=0.d0
dytot=0.d0
dztot=0.d0

!print*, 'tot_num_steps=',sum(nstep)
icount=0
open(12,file=outfile)
open(13,file=outposfile)
do i=1,nposfiles
   open(11,file=posfile(i))
   do j=1,nstep(i)
      icount=icount+1
      do k=1,nspecies
         if(.not.extractdisp(k))then
            do l=1,nat(k)
               read(11,*)
            enddo
         else
            do l=1,nat(k)
               if(icount.gt.1)then
                  xold(k,l)=xnew(k,l)
                  yold(k,l)=ynew(k,l)
                  zold(k,l)=znew(k,l)
               endif
               read(11,*)xnew(k,l),ynew(k,l),znew(k,l)
               write(13,*)xnew(k,l),ynew(k,l),znew(k,l)
               if(icount.gt.1)then
                  dx=xnew(k,l)-xold(k,l)
                  dy=ynew(k,l)-yold(k,l)
                  dz=znew(k,l)-zold(k,l)
                  if(dx.gt.halfbox)dx=dx-lbox
                  if(dy.gt.halfbox)dy=dy-lbox
                  if(dx.lt.(-halfbox))dx=dx+lbox
                  if(dy.lt.(-halfbox))dy=dy+lbox
                  write(12,*)dx,dy,dz
               endif
               dxtot(k,l)=dxtot(k,l)+dx
               dytot(k,l)=dytot(k,l)+dy
               dztot(k,l)=dztot(k,l)+dz
            enddo
         endif
      enddo
   enddo
   close(11)
enddo
close(12)
close(13)
open(14,file=out_all_disp_file)
icount=0
do k=1,nspecies
   if(extractdisp(k))then
   do l=1,nat(k)
      icount=icount+1
      write(14,*)icount,dxtot(k,l),dytot(k,l),dztot(k,l)
   enddo
   endif
enddo
close(14)

end program

