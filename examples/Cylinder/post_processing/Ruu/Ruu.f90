  program autocorrelation

  implicit none
 
  character(7) :: chits
  integer :: dig1, dig2, dig3, dig4, dig5, dig6, dig7
  integer :: nfiles, icrfile, file1, filen, ifile, nts,nss, nz,nzb, cont, i,j,k,xi,xf,nsonda
  real(4), allocatable :: sondex(:,:,:),sondey(:,:,:),sondez(:,:,:), u(:,:),v(:,:), ul(:,:), vl(:,:), R(:,:)
  real(4) :: Lz, Rtime, ul2m, vl2m
  
  ! Size of the file:
  nts=100000; nss=63; nzb=2; nsonda=46; nz=128
  xi=18001; xf=69000; Lz=6.
  nfiles=64; icrfile=1; file1=1; filen=64
   

  allocate ( sondex(nts,nss,nzb),sondey(nts,nss,nzb),sondez(nts,nss,nzb))
  allocate (u(nts,nz), v(nts,nz), ul(nts,nz), vl(nts,nz), R(1:3,1:nz/2))
  u(:,:)=0.0; v(:,:)=0.0; ul(:,:)=0.0; vl(:,:)=0.0; R(:,:)=0.0
  sondex(:,:,:)=0.0; sondey(:,:,:)=0.0; sondez(:,:,:)=0.0

  cont=0
  do ifile = file1,filen,icrfile

    dig1 =   ifile/1000000 + 48
    dig2 = ( ifile - 1000000*( ifile/1000000 ) )/100000 + 48
    dig3 = ( ifile - 100000*( ifile/100000 ) )/10000 + 48
    dig4 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
    dig5 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
    dig6 = ( ifile - 100*( ifile/100 ) )/10 + 48
    dig7 = ( ifile - 10*( ifile/10 ) )/1 + 48
    chits(1:7) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//char(dig5)//char(dig6)//char(dig7)
    write (*,*) 'file ','sondef',chits


    OPEN(10,FILE='../../dados/sondas/'//'sondef'//chits,form='unformatted', STATUS='unknown')
    READ(10) sondex,sondey,sondez
    CLOSE(10)

    do i=1,nzb
       u(:,cont+i)=sondex(:,nsonda,i)
       v(:,cont+i)=sondey(:,nsonda,i)
    enddo
    cont=cont+nzb
  enddo  
    
  ul(xi:xf,:)=u(xi:xf,:) - sum(u(xi:xf,:))/((xf-xi+1)*nz)
  vl(xi:xf,:)=v(xi:xf,:) - sum(v(xi:xf,:))/((xf-xi+1)*nz)

  print *,'umed,vmed=',sum(u(xi:xf,:))/((xf-xi+1)*nz),sum(v(xi:xf,:))/((xf-xi+1)*nz) 

  Rtime=0.0; ul2m=0.0
  do j=xi,xf
     do i=1,nz/2-1
        Rtime = Rtime + ul(j,i)*ul(j,i)
     enddo
     ul2m=ul2m+Rtime/(nz/2-1)
     Rtime=0.
  enddo
  ul2m=ul2m/(xf-xi+1)

  Rtime=0.0; vl2m=0.0
  do j=xi,xf
     do i=1,nz/2-1
        Rtime = Rtime + vl(j,i)*vl(j,i)
     enddo
     vl2m=vl2m+Rtime/(nz/2-1)
     Rtime=0.
  enddo
  vl2m=vl2m/(xf-xi+1)  

  cont=0; Rtime=0.0
  do k=1,nz/2
     R(1,k) = (k-1)*Lz/(nz-1)
     do j=xi,xf
        do i=1,nz/2-1
           Rtime = Rtime + ul(j,i)*ul(j,i+cont)
        enddo
        R(2,k) = R(2,k)+Rtime/(nz/2-1)
        Rtime=0.0
     enddo
     R(2,k) = R(2,k)/((xf-xi+1)*ul2m)
     cont=cont+1
  enddo

  cont=0; Rtime=0.0
  do k=1,nz/2
     do j=xi,xf
        do i=1,nz/2-1
           Rtime = Rtime + vl(j,i)*vl(j,i+cont)
        enddo
        R(3,k) = R(3,k)+Rtime/(nz/2-1)
        Rtime=0.0
     enddo
     R(3,k) = R(3,k)/((xf-xi+1)*vl2m)
     cont=cont+1
  enddo

  open(11,file='u',form='unformatted',status='unknown')
  write(11) u
  close(11)
  open(12,file='v',form='unformatted',status='unknown')
  write(12) v
  close(12)
  open(11,file='ul',form='unformatted',status='unknown')
  write(11) ul
  close(11)
  open(12,file='vl',form='unformatted',status='unknown')
  write(12) vl
  close(12)

  open (11,file='Ruu' ,form='formatted',status='unknown')
  write (11,FMT='(3ES16.7)') R
  close(11)

  end program autocorrelation
