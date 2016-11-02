program pm2xyz

implicit none 

integer :: nspecmax
parameter (nspecmax=10)
integer :: nconfigs,i,j,nspec,nionstot,nions(nspecmax),k,nskip
double precision :: x,y,z,conv
character*3 :: spc(nspecmax)
character*80 :: filein, fileout

write(6,*)'How many configs?'
read(5,*)nconfigs
print *, "nconfigs=", nconfigs
write(6,*)'How many steps do you skip?'
read(5,*)nskip
print *, "nskip=", nskip
write(6,*)'How many species?'
read(5,*)nspec
print *, "nspec=", nspec
if(nspec.gt.nspecmax)then
   write(6,*)'Too many species -- modify your program'
endif

nionstot=0
do i=1,nspec
     write(6,*)'specie',i
     write(6,*)'Chemical element?'
     read(5,*)spc(i)
     print'(A i1 A A)', 'specie ',i, '=',spc(i) 
     write(6,*)'How many ions?'
     read(5,*)nions(i)
     print'(A i1 A i5 A)', 'specie ',i, 'has ',nions(i), 'ions' 
     nionstot=nionstot+nions(i)
enddo
write(6,*)'Input file?'
read(5,*)filein
print *, "file input=", filein
open(10,file=filein)
write(6,*)'Output file?'
read(5,*)fileout
print *, "file output=", fileout
write(6,*)'Conversion? 1.0 si non'
read(5,*)conv
print *, "conversion=", conv

open(20,file=fileout)

do i=1,nconfigs
   if(mod(i,nskip).eq.0) write(20,*)nionstot
   if(mod(i,nskip).eq.0) write(20,*)'pas',i
   do k=1,nspec
      do j=1,nions(k)
         read(10,*) x,y,z
         if(mod(i,nskip).eq.0)write(20,*)spc(k),x*conv,y*conv,z*conv
      enddo
   enddo
enddo

close(10)
close(20)
print *, "DONE!"

end program
