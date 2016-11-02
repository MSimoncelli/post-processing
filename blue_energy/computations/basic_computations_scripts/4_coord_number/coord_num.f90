REAL(8) FUNCTION dist_w_PBC_opt(coor_a,coor_b,cell_size) RESULT(res)
	IMPLICIT NONE
	REAL(8), DIMENSION(3), INTENT(IN) :: coor_a
	REAL(8), DIMENSION(3), INTENT(IN) :: coor_b
	REAL(8), DIMENSION(3), INTENT(IN) :: cell_size

	integer :: i
	REAL(8), DIMENSION(9,3) :: coor_b_rep
	REAL(8), DIMENSION(9) :: dists
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


program coord_num

implicit none 

integer :: nspecmax, num_line_pos, io,m,n,nspec, base_shift_ion, AllocateStatus
parameter (nspecmax=10)
integer :: nconfigs,i,j,n_part_tot,k,nskip, n_ions, offset, counter
integer, dimension(64) :: buffer_num
character buffer_chr(64)*3
integer :: unit_1, unit_2,g,w, ios
double precision :: x,y,z,conv
character*80 :: filein, fileout
character*200 :: inputline
character*200 :: junk
real(8) :: bohr
real(8), allocatable :: r_cut_off(:,:)
real(8), allocatable :: theta(:)
real(8), dimension(16818, 3) :: array_conf !with allocate it does not work!!!
!real(8), allocatable :: array_conf(:,:)
integer, allocatable :: ar_out(:,:)
integer, dimension(4) :: n_specie_dist
!real(8), allocatable :: ion_p(:,:)
real(8), dimension(320, 3) :: ion_p
integer, dimension(nspecmax) :: n_part
integer, dimension(nspecmax) :: array_index
character*3 spc(nspecmax)
integer :: atoms_per_conf, num_configs, rel_atoms_per_dist, n_dist, specie
integer, EXTERNAL :: num_lines_file
real(8), EXTERNAL :: dist_w_PBC_opt
real(8), dimension(3) :: box_size
real(8), dimension(3) :: coor_1, coor_2
character*60 :: crap
real(8) :: buf_dist, PI, dcc, confdegrenorm
character :: tb 
CHARACTER(len=64) ::nome_file_cutoffs

PI=4.0*atan(1.0)
! For graphite, in bohr.
! Approximate carbon-carbon distance for a structure close to a graphitic structure
! Used to calculate the solid angle 
dcc=2.70 != 1.43/conv
! Confinement degre equal : total solid angle/maximum solid angle
! with maximum solid angle being 4*PI*1.65 because carbon surface does not cover all hexagonal surface
! See Merlet et al., Nature Commun., 4, 2701 (2013) for more details, in particular SI.
confdegrenorm=100*3.0*sqrt(3.0)/(4.0*PI*PI)

tb=char(9)

open (unit = 17, file = "positions.out",status='old', iostat=io)
open (unit = 18, file = "runtime.inpt",status='old', iostat=io)
          
DO i = 1, iargc()
  CALL getarg(i,nome_file_cutoffs)
END DO

num_line_pos=0
call read_sim_box_infos(box_size)
!print *, box_size

! read the number of lines of the file
num_line_pos = num_lines_file(17)

!print *, "num_line_pos=", num_line_pos

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
	!print *, 'spc(', k, ')=', spc(k)
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
!print *, 'num_atoms_pr_conf=', atoms_per_conf
!print *, 'num_configs=', num_configs
!print *, 'ntot_ions=',n_ions

! read relevant data in a configuration and store in an array
! for water you consider only oxigen!
rel_atoms_per_dist=atoms_per_conf-n_part(3)-n_part(2)

n_dist=4

fileout='out_coord_num.dat'
open(unit=22,file=fileout,iostat=io)

write(22, *)tb,'#label',tb,'coor x',tb,'  coor y',tb,'coor z',tb,'N(O)',tb,'N(Na)',tb,'N(Cl)',tb,'N(C)'

allocate(r_cut_off(2,n_dist))
allocate(ar_out(n_ions,n_dist))
allocate(theta(n_ions))

do i=1,n_ions
!print *, 'init_theta', i
theta(i)=0.0
end do
!print *, 'init ok'

!allocate(ion_p(n_ions,3))
!print *, rel_atoms_per_dist
!allocate(array_conf(rel_atoms_per_dist,3))
   
ar_out(:,:)=0
! 1st index:
! 1= Na
! 2= Cl

! 2nd index
! 1= O
! 2= Na
! 3= Cl
! 4= C

! to convert from angstrom to atomic units
bohr= 0.529177

open(15,file=nome_file_cutoffs)
read(15,*)crap
!print *, crap
read(15,*)crap
!print *, crap
read(15,*)crap
!print *, crap
!Na-O
read(15,*)r_cut_off(1,1)
!print *, r_cut_off(1,1)
r_cut_off(1,1)=r_cut_off(1,1)/bohr
read(15,*)crap
!Na-Na
read(15,*)crap
read(15,*)r_cut_off(1,2)
!print *, r_cut_off(1,2)
r_cut_off(1,2)=r_cut_off(1,2)/bohr
read(15,*)crap
!Na-Cl
read(15,*)crap
read(15,*)r_cut_off(1,3)
!print *, r_cut_off(1,3)
r_cut_off(1,3)=r_cut_off(1,3)/bohr
read(15,*)crap
!Na-C
read(15,*)crap
read(15,*)r_cut_off(1,4)
!print *,r_cut_off(1,4)
r_cut_off(1,4)=r_cut_off(1,4)/bohr
read(15,*)crap
!Cl-O
read(15,*)crap
read(15,*)r_cut_off(2,1)
!print *,r_cut_off(2,1)
r_cut_off(2,1)=r_cut_off(2,1)/bohr
read(15,*)crap
!Cl-Na
read(15,*)crap
read(15,*)r_cut_off(2,2)
!print *,r_cut_off(2,2)
r_cut_off(2,2)=r_cut_off(2,2)/bohr
read(15,*)crap
!Cl-Cl
read(15,*)crap
read(15,*)r_cut_off(2,3)
!print *,r_cut_off(2,3)
r_cut_off(2,3)=r_cut_off(2,3)/bohr
read(15,*)crap
!Cl-C
read(15,*)crap
read(15,*)r_cut_off(2,4)
!print *,r_cut_off(2,4)
r_cut_off(2,4)=r_cut_off(2,4)/bohr


!print *,'done'
!print *, 'r_cut_off=',r_cut_off
!print *, 'space'
!print *, r_cut_off

n_specie_dist(1)=n_part(1) !Oxigen
n_specie_dist(2)=n_part(4) !Na
n_specie_dist(3)=n_part(5) !Cl
n_specie_dist(4)=n_part(6)+n_part(7)+n_part(8) !C1+C2+P	


!index: O, H1, H2, Na, Cl, C1, C2, P
base_shift_ion=n_part(1)+n_part(2)+n_part(3)
!print *, base_shift_ion


do w=1,num_configs
	!print '(A i2 A i5)', 'frame ', w, ' out of ', num_configs
	g=1
	do i=1, atoms_per_conf
		if ((i.le.n_part(1)).or.(i.gt.base_shift_ion)) then
			read(17, *,IOSTAT=io)  array_conf(g,1), array_conf(g,2), array_conf(g,3)  
			g=g+1
		else
			read(17, *,IOSTAT=io) crap ! the reading continues
		end if
	end do

	specie=1
	do i=1,n_ions
		offset=0
		ion_p(i,1)=array_conf(n_part(1)+i,1) !select one ion
		ion_p(i,2)=array_conf(n_part(1)+i,2)
		ion_p(i,3)=array_conf(n_part(1)+i,3)
		coor_1(1)=ion_p(i,1)
		coor_1(2)=ion_p(i,2)
		coor_1(3)=ion_p(i,3)
		if (i<=n_ions/2) then
			specie=1
			! remove the overcounting (the reference ion must not contribute!)
			ar_out(i,2)=ar_out(i,2)-1
		else 
			specie=2
			! remove the overcounting (the reference ion must not contribute!)
			ar_out(i,3)=ar_out(i,3)-1
		end if 
		do k=1,n_dist	
			do j=1,n_specie_dist(k)
				coor_2(1)=array_conf(offset+j,1) !select one ion
				coor_2(2)=array_conf(offset+j,2)
				coor_2(3)=array_conf(offset+j,3)
				buf_dist=dist_w_PBC_opt(coor_1,coor_2,box_size)
				!counter=counter+1
				!buf_dist=dist_w_PBC(coor_1,coor_2,box_size)
				if (buf_dist.lt.r_cut_off(specie,k)) then
					ar_out(i,k)=ar_out(i,k)+1
					if (k.eq.4) then
					theta(i)=theta(i)+2*PI*(1-(buf_dist/sqrt(buf_dist*buf_dist+(dcc/2.0)*(dcc/2.0))))*100.0/(4.0*PI*0.6046)
					!print *,'cycle_in', i 
					end if
				end if
			end do
			offset=offset+n_specie_dist(k)
		end do
	end do

	write(22,'(A i3)' ) '#Run ', w
	do i=1,n_ions
		if (i.le.n_ions/2) then
			!write(22,'(A A i4 A i4 A i4 A i4 A i4)')'Na',tb,i,tb, ar_out(i,1),tb,ar_out(i,2),tb, ar_out(i,3),tb,ar_out(i,4)
			write(22,'(A A i4 f10.3 f10.3 f10.3 A i4 A i4 A i4 A i4 A f10.3)')'Na',tb,i,ion_p(i,:), &
			tb, ar_out(i,1),tb,ar_out(i,2),tb, ar_out(i,3),tb,ar_out(i,4),tb, theta(i)
		else
			write(22,'(A A i4 f10.3 f10.3 f10.3 A i4 A i4 A i4 A i4 A f10.3)')'Cl',tb,i-160,ion_p(i,:), &
			tb, ar_out(i,1),tb,ar_out(i,2),tb, ar_out(i,3),tb,ar_out(i,4),tb, theta(i)
		end if
		theta(i)=0
	end do
	ar_out(:,:)=0
end do
!print *, 'Done!'

!deallocate(array_conf)
!deallocate(ion_p)
!deallocate(n_specie_dist)
deallocate(ar_out)
deallocate(r_cut_off)
close(22)
close(18)
close(17)
end program
