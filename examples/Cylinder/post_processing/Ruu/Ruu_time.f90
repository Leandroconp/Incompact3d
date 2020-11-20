  program autocorrelation

  implicit none
 
  character(7) :: chits
  integer :: dig1, dig2, dig3, dig4, dig5, dig6, dig7
  integer :: nfiles, icrfile, file1, filen, ifile, nts,nss, nz,nzb, cont, i,j,k,xi,xf,nsonda,var
  real(4) , allocatable :: sondex(:,:,:),sondey(:,:,:),sondez(:,:,:), u(:,:),v(:,:), ul(:,:), vl(:,:), R(:,:,:), Nlength(:)
  real(4) , allocatable :: hist(:,:), Rmed(:,:), ul2t(:), vl2t(:)
  real(4) :: Lz, Rtime, ul2m, vl2m, eps, vmach, epmach, length, ulmed
  
  ! Size of the file:
  nts=80000; nss=63; nzb=2; nsonda=38; nz=128
  xi=35000; xf=60000; Lz=6.
  nfiles=64; icrfile=1; file1=1; filen=64
  var=2 !number of variables
   

  allocate ( sondex(nts,nss,nzb),sondey(nts,nss,nzb),sondez(nts,nss,nzb))
  allocate (u(nts,nz), v(nts,nz), ul(nts,nz), vl(nts,nz), R(nts,nz/2,var))
  allocate ( Nlength(nts), ul2t(nts), vl2t(nts) )
  allocate ( hist(2,nz/2) , Rmed(2,nz/2))

  u(:,:)=0.0; v(:,:)=0.0; ul(:,:)=0.0; vl(:,:)=0.0; R(:,:,:)=0.0
  sondex(:,:,:)=0.0; sondey(:,:,:)=0.0; sondez(:,:,:)=0.0; Rmed(:,:)=0.0

  cont=0
!  do ifile = file1,filen,icrfile
!
!    dig1 =   ifile/1000000 + 48
!    dig2 = ( ifile - 1000000*( ifile/1000000 ) )/100000 + 48
!    dig3 = ( ifile - 100000*( ifile/100000 ) )/10000 + 48
!    dig4 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
!    dig5 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
!    dig6 = ( ifile - 100*( ifile/100 ) )/10 + 48
!    dig7 = ( ifile - 10*( ifile/10 ) )/1 + 48
!    chits(1:7) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//char(dig5)//char(dig6)//char(dig7)
!    write (*,*) 'file ','sondef',chits
!
!
!    OPEN(10,FILE='../../dados/sondas3/'//'sondef'//chits,form='unformatted', STATUS='unknown')
!    READ(10) sondex,sondey,sondez
!    CLOSE(10)
!
!    do i=1,nzb
!       u(:,cont+i)=sondex(:,nsonda,i)
!       v(:,cont+i)=sondey(:,nsonda,i)
!    enddo
!    cont=cont+nzb
!  enddo
 
!   open(12,file='u',form='unformatted',status='unknown')
!   write(12) u
!   close(12)
!   open(12,file='v',form='unformatted',status='unknown')
!   write(12) v
!   close(12)
 
  OPEN(10,FILE='u',form='unformatted', STATUS='unknown')
  READ(10) u
  CLOSE(10)
  OPEN(10,FILE='v',form='unformatted', STATUS='unknown')
  READ(10) v
  CLOSE(10)
 

  ! fluctuation ul and vl 
  do j=xi,xf  
     ul(j,1:nz/2)=u(j,1:nz/2) - sum(u(j,1:nz/2))/(nz/2)
     vl(j,1:nz/2)=v(j,1:nz/2) - sum(v(j,1:nz/2))/(nz/2)
  enddo

  ! autocorelation of ul
  ul2m=0.0; vl2m=0.0; cont=0; Rtime=0.0
  do k=1,nz/2
     do j=xi,xf
        do i=1,nz/2-1
           Rtime = Rtime + ul(j,i)*ul(j,i+cont)
           ul2m = ul2m + ul(j,i)*ul(j,i)
        enddo
        R(j,k,1) = R(j,k,1)+Rtime/ul2m
        Rtime=0.0
        ul2m=0.0
     enddo
     cont=cont+1
  enddo

  ! autocorelation of vl
  cont=0; Rtime=0.0; 
  do k=1,nz/2
     do j=xi,xf
        do i=1,nz/2-1
           Rtime = Rtime + vl(j,i)*vl(j,i+cont)
           vl2m = vl2m + vl(j,i)*vl(j,i)
        enddo
        R(j,k,2) = R(j,k,2)+Rtime/vl2m
        Rtime=0.0
        vl2m=0.0 
     enddo
     cont=cont+1
  enddo

   ! Error of machine
   eps=1.
1  eps=eps/2.
   vmach=eps+1.
   do while(vmach > 1.)
      go to 1
   enddo
   epmach=abs(eps) ! Error machine


  ! Wavelength estimation by ul, if want vl change all R(j,i,1) to R(j,i,2)
  Nlength(:)=0.0
  do j=xi,xf
     i=3
     do while (Nlength(j).eq.0.0)
        if (R(j,i,1).gt.R(j,i-1,1)+epmach.and.R(j,i,1).gt.R(j,i+1,1)+epmach) then
           Nlength(j)=(i-1)*Lz/(nz-1)
        else if (i.eq.nz/2-1) then 
           Nlength(j)=Lz/2 + Lz/(nz-1)
        end if
        i=i+1
     end do
  enddo

  ! averaged fluctuation for ul
  Rtime=0.0; ul2m=0.0; ul2t=0.0
  do j=xi,xf
     do i=1,nz/2
        Rtime = Rtime + ul(j,i)*ul(j,i)
     enddo
     ul2t(j)=sqrt(Rtime/(nz/2))
     Rtime=0.
  enddo
  ul2m=maxval(ul2t(xi:xf))/2.

  ! averaged fluctuation for vl
  Rtime=0.0; vl2m=0.0; vl2t=0.0
  do j=xi,xf
     do i=1,nz/2
        Rtime = Rtime + vl(j,i)*vl(j,i)
     enddo
     vl2m=sqrt(Rtime/(nz/2))
     Rtime=0.
  enddo
  vl2m=maxval(vl2t(:))/2 
  print *, ul2m,vl2m

  ! Histogram of wavelength
  hist(:,:)=0.;cont=0; length=0.
  do i=1,nz/2
     hist(1,i)=(i)*Lz/(nz-1)
     do j=xi,xf
        if (sqrt(sum(ul(j,1:nz/2)*ul(j,1:nz/2))/(nz/2)).gt.0.5*ul2m) then 
           if (Nlength(j).le.(i)*Lz/(nz-1).and.Nlength(j).gt.(i-1)*Lz/(nz-1)) then
              hist(2,i)=hist(2,i)+1.
              cont=cont+1
              length=length+Nlength(j)
           endif
        end if
     enddo
  enddo
  hist(2,:)=hist(2,:)*100/cont
  length=length/cont
  print *,'wavelength=',length

  
  open (10,file='hist' ,form='formatted',status='unknown')
  write (10,FMT='(2ES16.7)') hist
 ! close(10)
 ! open (11,file='Nlength' ,form='formatted',status='unknown')
 ! write (11,FMT='(1ES16.7)') Nlength
 ! close(11)
  do i=1,nz/2
   Rmed(1,i)=(i-1)*Lz/(nz-1)
  enddo
   Rmed(2,:)=R(48500,:,1)
  open (11,file='Ruut97' ,form='formatted',status='unknown')
  write (11,FMT='(2ES16.7)') Rmed
  close(11)
  open(12,file='Ruutime',form='unformatted',status='unknown')
  write(12) R(:,:,1)
  close(12)
  open(11,file='ul',form='unformatted',status='unknown')
  write(11) ul
  close(11)
 ! open(12,file='v',form='unformatted',status='unknown')
 ! write(12) v
 ! close(12)
  

  end program autocorrelation
