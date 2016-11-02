program msdlayers

implicit none

integer :: i,j,k
integer :: ntrajpieces_1,nbin,npts_1,ibin_1,nmsdcalltime
integer :: ntrajpieces_2,npts_2,ibin_2,nmsdlen
integer :: nt,mcorrtime

double precision :: du,dtime,time
double precision :: xdisp_1,ydisp_1,zdisp_1
double precision :: xdisp_2,ydisp_2,zdisp_2
double precision, allocatable, dimension(:) :: xdispstore,ydispstore,zdispstore
double precision, allocatable, dimension(:) :: xmsd,ymsd,zmsd,msdtot,msdxy,norm
logical :: overflow

character*500:: dispfile_1,keyfile_1,outfile_1
character*500:: dispfile_2,keyfile_2,outfile_2
character*500:: nome_file_input_1,nome_file_input_2


CALL getarg(1,nome_file_input_1)
CALL getarg(2,nome_file_input_2)

!print*, nome_file_input_1
!print*, '---'
!print*, nome_file_input_2
open(10,file=nome_file_input_1)

read(10,*)ntrajpieces_1
!print*, ntrajpieces_1
read(10,*)nbin
read(10,*)nmsdlen
read(10,*)nmsdcalltime,dtime
read(10,'(a)')dispfile_1
open(11,file=dispfile_1)
read(10,'(a)')keyfile_1
open(12,file=keyfile_1)
read(10,'(a)')outfile_1
open(21,file=outfile_1)

close(10)

open(30,file=nome_file_input_2)

read(30,*)ntrajpieces_2
!print*, ntrajpieces_2
read(30,*)nbin
read(30,*)nmsdlen
read(30,*)nmsdcalltime,dtime
read(30,'(a)')dispfile_2
open(31,file=dispfile_2)
read(30,'(a)')keyfile_2
open(32,file=keyfile_2)
read(30,'(a)')outfile_2
!open(41,file=outfile_2)

close(30)

!print*, 'nmsdlen=',nmsdlen

allocate(xmsd(0:nmsdlen),ymsd(0:nmsdlen),zmsd(0:nmsdlen))
allocate(norm(0:nmsdlen),msdtot(0:nmsdlen),msdxy(0:nmsdlen))
allocate(xdispstore(nmsdlen),ydispstore(nmsdlen),zdispstore(nmsdlen))

xdispstore=0.d0
ydispstore=0.d0
zdispstore=0.d0
xmsd=0.d0
ymsd=0.d0
zmsd=0.d0
norm=0.d0


do i=1,max(ntrajpieces_1,ntrajpieces_2)
   !write(6,*)i,'over',max(ntrajpieces_1,ntrajpieces_2),'done'
     if (i.lt.ntrajpieces_1) then
          read(12,*)ibin_1,npts_1
     end if
     if (i.lt.ntrajpieces_2) then
          read(32,*)ibin_2,npts_2
     end if
     if ((ibin_1.ne.nbin).or.(ibin_2.ne.nbin)) then
            if ((i.lt.ntrajpieces_1).and.(ibin_1.ne.nbin)) then
                      do j=1,npts_1
                        read(11,*)du,du,du
                      enddo
            end if
            if  ((i.lt.ntrajpieces_2).and.(ibin_2.ne.nbin)) then
                  do j=1,npts_2
                      read(31,*)du,du,du
                  enddo 
            end if
    else
      mcorrtime=1
      overflow=.false.
            do j=1,max(npts_1,npts_2)
                    xdisp_1=0
                    ydisp_1=0
                    zdisp_1=0
                    xdisp_2=0
                    ydisp_2=0
                    zdisp_2=0
                    if ((j.le.npts_1).and.(i.lt.ntrajpieces_1)) then 
                       read(11,*)xdisp_1,ydisp_1,zdisp_1
                    end if
                    if ((j.le.npts_2).and.(i.lt.ntrajpieces_2)) then 
                     read(31,*)xdisp_2,ydisp_2,zdisp_2
                    end if

                     xdispstore(mcorrtime)=0.0d0             
                     ydispstore(mcorrtime)=0.0d0             
                     zdispstore(mcorrtime)=0.0d0  
           
                     do k=1,mcorrtime
                        xdispstore(k)=xdispstore(k)+xdisp_1-xdisp_2
                        ydispstore(k)=ydispstore(k)+ydisp_1-ydisp_2
                        zdispstore(k)=zdispstore(k)+zdisp_1-zdisp_2
                     enddo
                     if(overflow)then
                           if (i.lt.ntrajpieces_1) then
                              do k=mcorrtime+1,nmsdlen
                                 xdispstore(k)=xdispstore(k)+xdisp_1
                                 ydispstore(k)=ydispstore(k)+ydisp_1
                                 zdispstore(k)=zdispstore(k)+zdisp_1
                              enddo
                            end if
                            if (i.lt.ntrajpieces_2) then
                              do k=mcorrtime+1,nmsdlen
                                 xdispstore(k)=xdispstore(k)-xdisp_2
                                 ydispstore(k)=ydispstore(k)-ydisp_2
                                 zdispstore(k)=zdispstore(k)-zdisp_2
                              enddo
                            end if
                     endif
             
                     do k=1,mcorrtime
                        nt=mcorrtime-k
                        xmsd(nt)=xmsd(nt)+xdispstore(k)**2.0d0
                        ymsd(nt)=ymsd(nt)+ydispstore(k)**2.0d0
                        zmsd(nt)=zmsd(nt)+zdispstore(k)**2.0d0
                        norm(nt)=norm(nt)+1.0d0
                     enddo
                     if(overflow)then
                        do k=mcorrtime+1,nmsdlen
                           nt=mcorrtime-k+nmsdlen
                           xmsd(nt)=xmsd(nt)+xdispstore(k)**2.0d0
                           ymsd(nt)=ymsd(nt)+ydispstore(k)**2.0d0
                           zmsd(nt)=zmsd(nt)+zdispstore(k)**2.0d0
                           norm(nt)=norm(nt)+1.0d0
                        enddo
                     endif

                     if(mod(float(mcorrtime),float(nmsdlen)).eq.0)then
                        overflow=.true.
                     endif
                     mcorrtime=int(mod(float(mcorrtime),float(nmsdlen)))
                     mcorrtime=mcorrtime+1
                enddo 
              endif
   read(11,*)
   read(31,*)
enddo


do i=0,nmsdlen-1
   !xmsd(i)=xmsd(i)/norm(i)
   !ymsd(i)=ymsd(i)/norm(i)
   !zmsd(i)=zmsd(i)/norm(i)
   !msdtot(i)=xmsd(i)+ymsd(i)+zmsd(i)
   !msdxy(i)=xmsd(i)+ymsd(i)
   time=(dble(i)+1)*dble(nmsdcalltime)*dtime*2.418d-5
   write(21,*)time,xmsd(i),ymsd(i)
   !write(21,*)time,msdxy_1(i),msdtot_1(i)
enddo
end program

