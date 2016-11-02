! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  real(8), intent(in out), dimension(:,:) :: A
  integer :: iq

  if(size(A(:,3)) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1,:))
     call QsortC(A(iq:,:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(8), intent(in out), dimension(:,:) :: A
  integer, intent(out) :: marker
  real(8), dimension(3) :: temp
  integer :: i, j
  real(8) :: x      ! pivot point
  x = A(1,3)
  i= 0
  j= size(A(:,3)) + 1

  do
     j = j-1
     do
        if (A(j,3) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i,3) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp(1) = A(i,1)
        temp(2) = A(i,2)
        temp(3) = A(i,3)
        A(i,:) = A(j,:)
        A(j,1)=temp(1)
        A(j,2)=temp(2)
        A(j,3)=temp(3)
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module

integer function num_lines_file(num_file) result(res)
implicit none 
integer, intent(in) ::num_file
integer :: io
character*200 :: inputline
integer :: n_lin
n_lin=0
do
read(num_file,*,IOSTAT=io)  inputline
if (io > 0) then
write(*,*) 'Check input.  Something was wrong'
exit
else if (io < 0) then
exit
else
n_lin = n_lin + 1
end if
end do
rewind(num_file)
res=n_lin
end function num_lines_file

subroutine read_sim_box_infos(memo_array_coord) 
REAL(8), DIMENSION(3), INTENT(INOUT) :: memo_array_coord
integer, EXTERNAL :: num_lines_file

integer :: length, i
character*60 :: crap
open (unit = 19, file = "restart.dat")
length=num_lines_file(19)
do i=1,(length-3)
 read(19,*) crap
end do
do i=1,3
read(19,'(f100.10)',IOSTAT=io)  memo_array_coord(i)
end do
close(19)
end subroutine read_sim_box_infos

REAL(8) FUNCTION dist_w_PBC_opt(coor_a,coor_b,cell_size) RESULT(res)
IMPLICIT NONE
REAL(8), DIMENSION(3), INTENT(IN) :: coor_a
REAL(8), DIMENSION(3), INTENT(IN) :: coor_b
REAL(8), DIMENSION(3), INTENT(IN) :: cell_size

real(8) :: dxcf, dycf,dzcf, halfboxxrec, halfboxyrec

halfboxxrec=2.0/cell_size(1)
halfboxyrec=2.0/cell_size(2)

dxcf=coor_a(1)-coor_b(1)
dycf=coor_a(2)-coor_b(2)
dzcf=coor_a(3)-coor_b(3)

! minimal distance convenction
dxcf=dxcf-cell_size(1)*int(dxcf*halfboxxrec)
dycf=dycf-cell_size(2)*int(dycf*halfboxyrec)


res=sqrt(dxcf**2+dycf**2+dzcf**2)
END FUNCTION dist_w_PBC_opt

program water_RDF
use qsort_c_module
implicit none 
integer :: w,nspecmax, num_line_pos, nspec,io 
integer :: num_configs, n_part_tot, n_ions,k
parameter (nspecmax=10)
integer, EXTERNAL :: num_lines_file
real(8), EXTERNAL :: dist_w_PBC_opt
integer, dimension(64) :: buffer_num
integer :: idx_z, update_1,update_2,update_3,update_4,update_5,update_6
character buffer_chr(64)*3
character*60 :: crap
character*80 ::  fileout
character*200 :: junk
real(8) :: conv
integer :: Nbins, i,j, atoms_per_conf, r, histo_pos
integer :: n_dist, n_spec_ions,num_regions,current_run
character(len=3) current_run_chr
!real(8), ALLOCATABLE, DIMENSION(3,100) :: histo
real(8), DIMENSION(3,100) :: histo
real(8), ALLOCATABLE, DIMENSION(:,:) :: array_conf
real(8), ALLOCATABLE, DIMENSION(:) :: num_density
real(8), ALLOCATABLE, DIMENSION(:) :: num_balls
real(8), DIMENSION(0:2,0:1) :: regions
integer, DIMENSION(6) :: indexes_regions
integer, dimension(nspecmax) :: n_part
integer, dimension(nspecmax) :: array_index
character*3 spc(nspecmax)
integer, allocatable :: n_specie_dist(:)
real(8), dimension(3) :: box_size, coor_a, coor_b
real(8) :: big_cut_off, buf_dist, vol_electrode_1, vol_electrode_2, vol_bulk, volume_cut_off
character :: tb 

tb=char(9)
conv= 0.529177 
Nbins=100
n_dist=4
n_spec_ions=2
num_regions=3 
!print *, '1st check'

allocate(n_specie_dist(1:n_dist))
!electrode 1, bulk, electrode 2

CALL GetArg(1, current_run_chr)
!print *, current_run_chr
read( current_run_chr, '(i10)' ) current_run


!allocate(histo(1:num_regions,1:Nbins))
allocate(num_density(1:num_regions))
allocate(num_balls(1:num_regions))

!initialize to zero.
do i=1,num_regions
num_density(i)=0.0
num_balls(i)=0.0
do j=1,Nbins
histo(i,j)=0.0
end do
end do

open (unit = 17, file = "positions.out",status='old', iostat=io)
open (unit = 18, file = "runtime.inpt",status='old', iostat=io)
  
! read the number of lines of the file
num_line_pos = num_lines_file(17)
do k=1,64
read(18,*,IOSTAT=io)  buffer_num(k)
end do
rewind(18)
do k=1,64
  read(18,'(A A)',IOSTAT=io)  buffer_chr(k), junk
end do

array_index = (/7,12,17,22,27,32,37,42,0,0/)
atoms_per_conf=0
n_part_tot=0
n_ions=0
nspec=5

do k=1,8
spc(k)=buffer_chr(array_index(k))
!!print *, 'spc(', k, ')=', spc(k)
n_part(k)=buffer_num(array_index(k)+1)
!print *, 'n_part(', k, ')=', n_part(k)
!if (((spc(k)=='Na') .or. (spc(k)=='Cl')) .or. (spc(k)=='K') ) then
if ((k==4) .or. (k==5)) then
    n_ions=n_ions+n_part(k)
  end if
atoms_per_conf=atoms_per_conf+n_part(k)
!print *, n_part(k)
end do 
num_configs=num_line_pos/atoms_per_conf

n_specie_dist(1)=n_part(1) !Oxigen
n_specie_dist(2)=n_part(4) !Na
n_specie_dist(3)=n_part(5) !Cl
n_specie_dist(4)=n_part(6)+n_part(7)+n_part(8) !C1+C2+P

allocate(array_conf(1:n_part(1),1:3))

regions(0,0)=(17.000-5)/conv
regions(0,1)=(40.000+5)/conv

!regions(1,0)=55.000/conv
regions(1,0)=70.000/conv

call read_sim_box_infos(box_size)
!symmetric reasoning
regions(1,1)=box_size(3)-regions(1,0)
regions(2,0)=box_size(3)-regions(0,1)
regions(2,1)=box_size(3)-regions(0,0)

!big_cut_off= 0.5*(regions(0,1)-regions(0,0))
big_cut_off=10.0/conv
!print *, 'OK_bef_LOOP W'
!num_configs=1
do w=1,num_configs
!print '(A i2 A i5)', 'frame ', w, ' out of ', num_configs
  do i=1, n_part(1)
  read(17, *,IOSTAT=io)  array_conf(i,1), array_conf(i,2), array_conf(i,3)  
  end do
  do i=1,atoms_per_conf-n_part(1)
  read(17, *,IOSTAT=io) crap ! the reading continues
  end do

  call QsortC(array_conf)

  update_1=1
  update_2=1
  update_3=1
  update_4=1
  update_5=1
  update_6=1
  do idx_z=1,n_part(1)
    if(update_1.eq.1) then
    indexes_regions(1)=idx_z  
    end if    

    if(update_2.eq.1) then
    indexes_regions(2)=idx_z
    end if    

    if(update_3.eq.1) then
    indexes_regions(3)=idx_z
    end if    

    if(update_4.eq.1) then
    indexes_regions(4)=idx_z
    end if    

    if(update_5.eq.1) then
    indexes_regions(5)=idx_z
    end if    

    if(update_6.eq.1) then
    indexes_regions(6)=idx_z
    end if

    if((array_conf(idx_z,3).ge.regions(0,0)) .and. (update_1==1)) then
    update_1=0
    end if
    if (array_conf(idx_z,3).ge.regions(0,1) .and. (update_2==1)) then
    update_2=0
    end if
    if (array_conf(idx_z,3).ge.regions(1,0) .and. (update_3==1)) then
    update_3=0
    end if
    if (array_conf(idx_z,3).ge.regions(1,1) .and. (update_4==1)) then
    update_4=0
    end if
    if (array_conf(idx_z,3).ge.regions(2,0) .and. (update_5==1)) then
    update_5=0
    end if
    if (array_conf(idx_z,3).ge.regions(2,1) .and. (update_6==1)) then
    update_6=0
    end if  
  end do

  do r=1,3
!print *, r
    do i=indexes_regions(r*2-1),indexes_regions(r*2)
      coor_a(1)=array_conf(i,1)
      coor_a(2)=array_conf(i,2)
      coor_a(3)=array_conf(i,3)
      if ((abs(coor_a(3)-regions(r-1,0)).ge.big_cut_off) .and. (abs(coor_a(3)-regions(r-1,1)).ge.big_cut_off)) then
         !num_density(r)=num_density(r)+1!take into accoutn the central particle
         do j=indexes_regions(r*2-1),i-1
           coor_b(1)=array_conf(j,1)
           coor_b(2)=array_conf(j,2)
           coor_b(3)=array_conf(j,3)
           buf_dist=dist_w_PBC_opt(coor_a,coor_b,box_size)
           if (buf_dist.le.big_cut_off) then
              num_density(r)=num_density(r)+1
              histo_pos=ceiling((buf_dist/big_cut_off)*(Nbins)) 
              histo(r,histo_pos)=histo(r,histo_pos)+1
           end if
         end do
         do j=i+1,indexes_regions(r*2)
            coor_b(1)=array_conf(j,1)
            coor_b(2)=array_conf(j,2)
            coor_b(3)=array_conf(j,3)
            buf_dist=dist_w_PBC_opt(coor_a,coor_b,box_size)
            if (buf_dist.le.big_cut_off) then
               num_density(r)=num_density(r)+1
               histo_pos=ceiling((buf_dist/big_cut_off)*(Nbins)) 
               histo(r,histo_pos)=histo(r,histo_pos)+1
            end if
         end do
       end if
    end do
end do

end do



do i=1,Nbins
histo(:,i)=histo(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
end do
!vol_electrode_1=box_size(1)*conv*box_size(2)*conv*(regions(0,1)-regions(0,0))*conv
!vol_bulk=box_size(1)*conv*box_size(2)*conv*(regions(1,1)-regions(1,0))*conv
!vol_electrode_2=box_size(1)*conv*box_size(2)*conv*(regions(2,1)-regions(2,0))*conv

!print *, vol_electrode_1
!print *, vol_bulk
!print *, vol_electrode_2

!print *, num_density(1)
!print *, num_density(2)
!print *, num_density(3)

volume_cut_off=((4.0/3.0)*3.141592*(big_cut_off*conv)**3)
num_density(1)=num_density(1)*1.0/(volume_cut_off)
num_density(2)=num_density(2)*1.0/(volume_cut_off)
num_density(3)=num_density(3)*1.0/(volume_cut_off)


do r=1,3
histo(r,:)=histo(r,:)*(1.0/(num_density(r)))
end do 

write(fileout,'(AAA)'),'run',current_run_chr,'RDF_W_W_for.dat'
open(unit=22,file=fileout,iostat=io)
write(22,'(A i3)')'#num_configs=',num_configs
do r=1,3
write(22,'(A i1)')'# RDF region',r
do i=1,Nbins
write(22,'(f10.6 A f10.6)') (i*1.0/(1.0*(Nbins))*big_cut_off*conv),tb, histo(r,i)
enddo
enddo
!print*, 'FINE'
close(22)
close(18)
close(17)
end program


