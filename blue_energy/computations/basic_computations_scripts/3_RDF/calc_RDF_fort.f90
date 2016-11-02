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

module all_my_subroutines
use qsort_c_module
implicit none
contains

subroutine read_sim_box_infos(memo_array_coord) 
REAL(8), DIMENSION(3), INTENT(INOUT) :: memo_array_coord
integer, EXTERNAL :: num_lines_file

integer :: length, i,io
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

subroutine find_indexes_regions(ar_c_SPC,indexes_regions_SPC,regions_SPC,len_array) 
implicit none 
 real(8), intent(in out), dimension(:,:) :: ar_c_SPC
  integer, intent(in) :: len_array
  integer, intent(in out), dimension(6) :: indexes_regions_SPC
  real(8), intent(in),DIMENSION(0:2,0:1) :: regions_SPC
  integer :: idx_z, update_1,update_2,update_3,update_4,update_5,update_6
  call QsortC(ar_c_SPC)

  update_1=1
  update_2=1
  update_3=1
  update_4=1
  update_5=1
  update_6=1
  do idx_z=1,len_array
    if(update_1.eq.1) then
    indexes_regions_SPC(1)=idx_z  
    end if    

    if(update_2.eq.1) then
    indexes_regions_SPC(2)=idx_z
    end if    

    if(update_3.eq.1) then
    indexes_regions_SPC(3)=idx_z
    end if    

    if(update_4.eq.1) then
    indexes_regions_SPC(4)=idx_z
    end if    

    if(update_5.eq.1) then
    indexes_regions_SPC(5)=idx_z
    end if    

    if(update_6.eq.1) then
    indexes_regions_SPC(6)=idx_z
    end if

    if((ar_c_SPC(idx_z,3).ge.regions_SPC(0,0)) .and. (update_1==1)) then
    update_1=0
    end if
    if (ar_c_SPC(idx_z,3).ge.regions_SPC(0,1) .and. (update_2==1)) then
    update_2=0
    end if
    if (ar_c_SPC(idx_z,3).ge.regions_SPC(1,0) .and. (update_3==1)) then
    update_3=0
    end if
    if (ar_c_SPC(idx_z,3).ge.regions_SPC(1,1) .and. (update_4==1)) then
    update_4=0
    end if
    if (ar_c_SPC(idx_z,3).ge.regions_SPC(2,0) .and. (update_5==1)) then
    update_5=0
    end if
    if (ar_c_SPC(idx_z,3).ge.regions_SPC(2,1) .and. (update_6==1)) then
    update_6=0
    end if
  end do
end subroutine find_indexes_regions

subroutine build_rdf(ar_S_1,idx_1,ar_S_2,idx_2,histo,num_dens_1_2,Nbins,regions,box_size,big_cut_off,flag_equal) 
  implicit none
  real(8), intent(in out), dimension(:,:) :: ar_S_1
  integer, intent(in out), dimension(6) :: idx_1
  real(8), intent(in out), dimension(:,:) :: ar_S_2
  integer, intent(in out), dimension(6) :: idx_2
  real(8), intent(in out), dimension(:,:) :: histo
  real(8), intent(in out), dimension(3) :: num_dens_1_2
  real(8), intent(in),DIMENSION(0:2,0:1) :: regions
  real(8), intent(in), dimension(3) :: box_size
  real(8), intent(in) :: big_cut_off
  integer, intent(in) :: Nbins
  integer, intent(in) :: flag_equal
  real(8), dimension(3) :: coor_a, coor_b
  real(8) :: buf_dist, conv
  real(8), EXTERNAL :: dist_w_PBC_opt
  integer :: r, histo_pos,i,j

  conv= 0.529177
  !big_cut_off= 0.5*(regions(0,1)-regions(0,0))
  
do r=1,3
!print *, idx_1(r*2-1),idx_1(r*2)
!print *, idx_2(r*2-1),idx_2(r*2)
    do i=idx_1(r*2-1),idx_1(r*2)
      !select first atom
      coor_a(1)=ar_S_1(i,1) 
      coor_a(2)=ar_S_1(i,2)
      coor_a(3)=ar_S_1(i,3)
      if ((abs(coor_a(3)-regions(r-1,0)).ge.big_cut_off) .and. (abs(coor_a(3)-regions(r-1,1)).ge.big_cut_off)) then
         if (flag_equal.eq.0) then
              !print*, 'index Na=', i
             !num_dens_1_2(r)=num_dens_1_2(r)+1!take into accoutn the central particle
             do j=idx_2(r*2-1),idx_2(r*2)
               !print*, 'ind Cl=', j
               coor_b(1)=ar_S_2(j,1)
               coor_b(2)=ar_S_2(j,2)
               coor_b(3)=ar_S_2(j,3)
               buf_dist=dist_w_PBC_opt(coor_a,coor_b,box_size)
               if (buf_dist.le.big_cut_off) then
                  num_dens_1_2(r)=num_dens_1_2(r)+1
                  histo_pos=ceiling((buf_dist/big_cut_off)*(Nbins)) 
                  histo(r,histo_pos)=histo(r,histo_pos)+1
               end if
             end do
         end if
         if (flag_equal.eq.1) then
             !print*, 'ERROR'
             !num_dens_1_2(r)=num_dens_1_2(r)+1
              do j=idx_2(r*2-1),i-1
              coor_b(1)=ar_S_2(j,1)
              coor_b(2)=ar_S_2(j,2)
              coor_b(3)=ar_S_2(j,3)
              buf_dist=dist_w_PBC_opt(coor_a,coor_b,box_size)
               if (buf_dist.le.big_cut_off) then
                 num_dens_1_2(r)=num_dens_1_2(r)+1
                 histo_pos=ceiling((buf_dist/big_cut_off)*(Nbins)) 
                 histo(r,histo_pos)=histo(r,histo_pos)+1
               end if
             end do

             do j=i+1,idx_2(r*2)
                coor_b(1)=ar_S_2(j,1)
                coor_b(2)=ar_S_2(j,2)
                coor_b(3)=ar_S_2(j,3)
                buf_dist=dist_w_PBC_opt(coor_a,coor_b,box_size)
                if (buf_dist.le.big_cut_off) then
                   num_dens_1_2(r)=num_dens_1_2(r)+1
                   histo_pos=ceiling((buf_dist/big_cut_off)*(Nbins)) 
                   histo(r,histo_pos)=histo(r,histo_pos)+1
                end if
              end do
         end if  
       end if
    end do
end do
end subroutine build_rdf


subroutine write_to_file(histo,name_file,current_run_chr,num_configs,Nbins, big_cut_off) 
implicit none
real(8), intent(in out), dimension(:,:) :: histo
character*5, intent(in) :: name_file
character*3, intent(in) :: current_run_chr
integer, intent(in) :: num_configs
integer, intent(in) :: Nbins
real(8), intent(in) :: big_cut_off
real(8) :: conv
character :: tb 
integer :: r,i,io
character*80 ::  fileout

tb=char(9)
conv= 0.529177
!print*, current_run_chr
!print*, name_file
write(fileout,'(AAAAA)'),'run',current_run_chr,'RDF_',name_file,'fin.dat'
!print*, fileout
open(unit=22,file=fileout,iostat=io)
write(22,'(A i3)')'#num_configs=',num_configs
do r=1,3
write(22,'(A i1)')'#RDF region',r
do i=1,Nbins
write(22,'(f10.6 A f10.6)') (i*1.0/(1.0*(Nbins))*big_cut_off*conv),tb, histo(r,i)
enddo
enddo
!print*, 'FINE'
close(22)
end subroutine write_to_file

end module all_my_subroutines

program water_RDF
use all_my_subroutines
implicit none 
integer :: w,nspecmax, num_line_pos, nspec,io, Nbins
integer :: num_configs, n_part_tot, n_ions,k
parameter (nspecmax=10)
parameter (Nbins=100)
integer, EXTERNAL :: num_lines_file
real(8), EXTERNAL :: dist_w_PBC_opt
integer, dimension(64) :: buffer_num

character buffer_chr(64)*3
character*60 :: crap
character*200 :: junk
real(8) :: conv
integer ::  i,j, atoms_per_conf, r, histo_pos, total
integer :: n_dist, n_spec_ions,num_regions,current_run
character(len=3) current_run_chr
!real(8), ALLOCATABLE, DIMENSION(3,100) :: histo

real(8), ALLOCATABLE, DIMENSION(:,:) :: ar_c_O
real(8), ALLOCATABLE, DIMENSION(:,:) :: ar_c_Na
real(8), ALLOCATABLE, DIMENSION(:,:) :: ar_c_Cl
real(8), ALLOCATABLE, DIMENSION(:,:) :: ar_c_C

integer, DIMENSION(6) :: ind_reg_O
integer, DIMENSION(6) :: ind_reg_Na
integer, DIMENSION(6) :: ind_reg_Cl
integer, DIMENSION(6) :: ind_reg_C

real(8), ALLOCATABLE, DIMENSION(:) :: num_density_C_Cl
real(8), ALLOCATABLE, DIMENSION(:) :: num_density_C_Na
real(8), ALLOCATABLE, DIMENSION(:) :: num_density_O_Cl
real(8), ALLOCATABLE, DIMENSION(:) :: num_density_O_Na
real(8), ALLOCATABLE, DIMENSION(:) :: num_density_Na_Cl
real(8), ALLOCATABLE, DIMENSION(:) :: num_density_Cl_Cl
real(8), ALLOCATABLE, DIMENSION(:) :: num_density_Na_Na

real(8), DIMENSION(3,Nbins) :: histo_C_Cl
real(8), DIMENSION(3,Nbins) :: histo_C_Na
real(8), DIMENSION(3,Nbins) :: histo_O_Cl
real(8), DIMENSION(3,Nbins) :: histo_O_Na
real(8), DIMENSION(3,Nbins) :: histo_Na_Cl
real(8), DIMENSION(3,Nbins) :: histo_Cl_Cl
real(8), DIMENSION(3,Nbins) :: histo_Na_Na

real(8), DIMENSION(0:2,0:1) :: regions

integer, dimension(nspecmax) :: n_part
integer, dimension(nspecmax) :: array_index
character*3 spc(nspecmax)
integer, allocatable :: n_specie_dist(:)
real(8), dimension(3) :: box_size
character :: tb 
real(8) :: big_cut_off, volume_cut_off
tb=char(9)

conv=0.529177
big_cut_off=10.0/conv !more appropriate (empirical value!)

call read_sim_box_infos(box_size)


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
allocate(num_density_C_Cl(1:num_regions))
allocate(num_density_C_Na(1:num_regions))
allocate(num_density_O_Cl(1:num_regions))
allocate(num_density_O_Na(1:num_regions))
allocate(num_density_Na_Cl(1:num_regions))
allocate(num_density_Cl_Cl(1:num_regions))
allocate(num_density_Na_Na(1:num_regions))





!initialize to zero.
do i=1,num_regions
   num_density_C_Cl(i)=0.0
   num_density_C_Na(i)=0.0
   num_density_O_Cl(i)=0.0
   num_density_O_Na(i)=0.0
   num_density_Na_Cl(i)=0.0
   num_density_Cl_Cl(i)=0.0
   num_density_Na_Na(i)=0.0
   do j=1,Nbins
      histo_C_Cl(i,j)=0.0
      histo_C_Na(i,j)=0.0
      histo_O_Cl(i,j)=0.0
      histo_O_Na(i,j)=0.0
      histo_Na_Cl(i,j)=0.0
      histo_Cl_Cl(i,j)=0.0
      histo_Na_Na(i,j)=0.0
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
!print*, num_line_pos
n_specie_dist(1)=n_part(1) !Oxigen
n_specie_dist(2)=n_part(4) !Na
n_specie_dist(3)=n_part(5) !Cl
n_specie_dist(4)=n_part(6)+n_part(7)+n_part(8) !C1+C2+P

allocate(ar_c_O(1:n_specie_dist(1),1:3))
allocate(ar_c_Na(1:n_specie_dist(2),1:3))
allocate(ar_c_Cl(1:n_specie_dist(3),1:3))
allocate(ar_c_C(1:n_specie_dist(4),1:3))



regions(0,0)=(17.000-10)/conv
regions(0,1)=(40.000+10)/conv
regions(1,0)=55.000/conv
!regions(1,0)=55.000/conv

!print*, box_size
!symmetric reasoning
regions(1,1)=box_size(3)-regions(1,0)
regions(2,0)=box_size(3)-regions(0,1)
regions(2,1)=box_size(3)-regions(0,0)



!print *, 'OK_bef_LOOP W'
!num_configs=3
!total=0
!print*, num_configs
do w=1,num_configs
!print*, w
!print '(A i2 A i5)', 'frame ', w, ' out of ', num_configs
  do i=1, n_part(1)
  !total=total+1
  read(17, *,IOSTAT=io)  ar_c_O(i,1), ar_c_O(i,2), ar_c_O(i,3)  
  end do
  do i=1,(n_part(2)+n_part(3)) !do not consider the Hidrogen!
  read(17, *,IOSTAT=io) crap ! the reading continues
  !total=total+1
  end do
  do i=1, n_part(4)
  read(17, *,IOSTAT=io)  ar_c_Na(i,1), ar_c_Na(i,2), ar_c_Na(i,3)  
  !total=total+1
  end do
  do i=1, n_part(5)
  read(17, *,IOSTAT=io)  ar_c_Cl(i,1), ar_c_Cl(i,2), ar_c_Cl(i,3)  
  !total=total+1  
  end do
  do i=1, n_part(6)+n_part(7)+n_part(8)
  read(17, *,IOSTAT=io)  ar_c_C(i,1), ar_c_C(i,2), ar_c_C(i,3)  
  !total=total+1
  end do

!print*, (1.0*total)/(1.0*atoms_per_conf)

  call find_indexes_regions(ar_c_O,ind_reg_O,regions,n_specie_dist(1))
  call find_indexes_regions(ar_c_Na,ind_reg_Na,regions,n_specie_dist(2))
  call find_indexes_regions(ar_c_Cl,ind_reg_Cl,regions,n_specie_dist(3))
  call find_indexes_regions(ar_c_C,ind_reg_C,regions,n_specie_dist(4))
!print*,'ind_reg_Na=', ind_reg_Na
!print*,'ind_reg_Cl=', ind_reg_Cl
!print*,'ind_reg_C=', ind_reg_C
!print*,'ind_reg_Na=', ind_reg_Na

 !(ar_S_1,idx_1,ar_S_2,idx_2,histo,num_dens_1_2,Nbins,regions,box_size,big_cut_off,flag_equal) 

  !C_Cl
  !C_Na
  !O_Cl
  !O_Na
  !Na_Cl
  !Cl_Cl
  !Na_Na
  !print*, 'C-CL'
  call build_rdf(ar_c_C,ind_reg_C,ar_c_Cl,ind_reg_Cl,histo_C_Cl,num_density_C_Cl,Nbins,regions,box_size,big_cut_off,0) 
  !print*, 'C-Na'
  call build_rdf(ar_c_C,ind_reg_C,ar_c_Na,ind_reg_Na,histo_C_Na,num_density_C_Na,Nbins,regions,box_size,big_cut_off,0) 
  !print*, 'O-Cl'
  call build_rdf(ar_c_O,ind_reg_O,ar_c_Cl,ind_reg_Cl,histo_O_Cl,num_density_O_Cl,Nbins,regions,box_size,big_cut_off,0) 
  !print*, 'O-Na'
  call build_rdf(ar_c_O,ind_reg_O,ar_c_Na,ind_reg_Na,histo_O_Na,num_density_O_Na,Nbins,regions,box_size,big_cut_off,0) 
  !print*, 'Na-Cl'
  call build_rdf(ar_c_Na,ind_reg_Na,ar_c_Cl,ind_reg_Cl,histo_Na_Cl,num_density_Na_Cl,Nbins,regions,box_size,big_cut_off,0) 
  !print*, 'Cl-Cl'
  call build_rdf(ar_c_Cl,ind_reg_Cl,ar_c_Cl,ind_reg_Cl,histo_Cl_Cl,num_density_Cl_Cl,Nbins,regions,box_size,big_cut_off,1) 
  !print*, 'Na-Na'
  call build_rdf(ar_c_Na,ind_reg_Na,ar_c_Na,ind_reg_Na,histo_Na_Na,num_density_Na_Na,Nbins,regions,box_size,big_cut_off,1) 

end do

!print*, 'out_sub'
!print*, num_density_O_Cl
do i=1,Nbins
  histo_C_Cl(:,i)=histo_C_Cl(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
  histo_C_Na(:,i)=histo_C_Na(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
  histo_O_Cl(:,i)=histo_O_Cl(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
  histo_O_Na(:,i)=histo_O_Na(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
  histo_Na_Cl(:,i)=histo_Na_Cl(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
  histo_Cl_Cl(:,i)=histo_Cl_Cl(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
  histo_Na_Na(:,i)=histo_Na_Na(:,i)/(4.0*3.141592*((conv*big_cut_off*(i)*1.0/Nbins)**2)*conv*big_cut_off*(1.0/Nbins))
end do

volume_cut_off=((4.0/3.0)*3.141592*(big_cut_off*conv)**3)

num_density_C_Cl=num_density_C_Cl*(1.0/volume_cut_off)
num_density_C_Na=num_density_C_Na*(1.0/volume_cut_off)
num_density_O_Cl=num_density_O_Cl*(1.0/volume_cut_off)
num_density_O_Na=num_density_O_Na*(1.0/volume_cut_off)
num_density_Na_Cl=num_density_Na_Cl*(1.0/volume_cut_off)
num_density_Cl_Cl=num_density_Cl_Cl*(1.0/volume_cut_off)
num_density_Na_Na=num_density_Na_Na*(1.0/volume_cut_off)

!print*, num_density_O_Cl
do r=1,3
  histo_C_Cl(r,:)=histo_C_Cl(r,:)*(1.0/num_density_C_Cl(r))
  histo_C_Na(r,:)=histo_C_Na(r,:)*(1.0/num_density_C_Na(r))
  histo_O_Cl(r,:)=histo_O_Cl(r,:)*(1.0/num_density_O_Cl(r))
  histo_O_Na(r,:)=histo_O_Na(r,:)*(1.0/num_density_O_Na(r))
  histo_Na_Cl(r,:)=histo_Na_Cl(r,:)*(1.0/num_density_Na_Cl(r))
  histo_Cl_Cl(r,:)=histo_Cl_Cl(r,:)*(1.0/num_density_Cl_Cl(r))
  histo_Na_Na(r,:)=histo_Na_Na(r,:)*(1.0/num_density_Na_Na(r))
end do

call write_to_file(histo_C_Cl,'C_Cl_',current_run_chr,num_configs,Nbins, big_cut_off)
call write_to_file(histo_C_Na,'C_Na_',current_run_chr,num_configs,Nbins, big_cut_off)
call write_to_file(histo_O_Cl,'O_Cl_',current_run_chr,num_configs,Nbins, big_cut_off)
call write_to_file(histo_O_Na,'O_Na_',current_run_chr,num_configs,Nbins, big_cut_off)
call write_to_file(histo_Na_Cl,'Na_Cl',current_run_chr,num_configs,Nbins, big_cut_off)
call write_to_file(histo_Cl_Cl,'Cl_Cl',current_run_chr,num_configs,Nbins, big_cut_off)
call write_to_file(histo_Na_Na,'Na_Na',current_run_chr,num_configs,Nbins, big_cut_off)

close(18)
close(17)
end program

