     !***************************************************************
     !*  Calculate the median value of an array with the Heapsort   *
     !*  method                                                     *
     !* ----------------------------------------------------------- *
     !* REFERENCE:                                                  *
     !*      "NUMERICAL RECIPES By W.H. Press, B.P. Flannery,       *
     !*       S.A. Teukolsky and W.T. Vetterling, Cambridge         *
     !*       University Press, 1986" [BIBLI 08].                   *
     !* ----------------------------------------------------------- *
     !* SAMPLE RUN:                                                 *
     !*                                                             *
     !* 407.8 192.8 851.1 604.4 932.3 799.4 914.5 965.8 453.7 295.1 *
     !* 154.5 977.4 410.2 916.2 934.7 504.8 823.2 225.2 456.6  49.0 *
     !* 933.5 663.0 335.3 346.6 568.7 956.1 654.7 300.7 379.6 591.9 *
     !* 992.9 689.6 644.7 305.4 148.2 257.2 664.6 612.1 713.0  99.7 *
     !*  46.5 167.6 984.6 847.2  55.4  82.7 999.0  10.7 877.7 929.4 *
     !* 398.1 972.8 874.1 755.1 472.1 122.8 671.4  35.5 128.8  76.8 *
     !* 454.2 959.2 510.1 791.3 122.8 176.6 237.9 995.8 548.3 309.8 *
     !* 162.6 996.5 750.0 250.6 577.7 761.1 101.9 797.1 539.0 723.5 *
     !*                                                             *
     !*  Median value:    558.5251                                  *
     !*                                                             *
     !*  Sorted table (Heapsort method):                            *
     !*                                                             *
     !*  10.7  35.5  46.5  49.0  55.4  76.8  82.7  99.7 101.9 122.8 *
     !* 122.8 128.8 148.2 154.5 162.6 167.6 176.6 192.8 225.2 237.9 *
     !* 250.6 257.2 295.1 300.7 305.4 309.8 335.3 346.6 379.6 398.1 *
     !* 407.8 410.2 453.7 454.2 456.6 472.1 504.8 510.1 539.0 548.3 *
     !* 568.7 577.7 591.9 604.4 612.1 644.7 654.7 663.0 664.6 671.4 *
     !* 689.6 713.0 723.5 750.0 755.1 761.1 791.3 797.1 799.4 823.2 *
     !* 847.2 851.1 874.1 877.7 914.5 916.2 929.4 932.3 933.5 934.7 *
     !* 956.1 959.2 965.8 972.8 977.4 984.6 992.9 995.8 996.5 999.0 *
     !*                                                             *
     !*                                                             *
     !*                         F90 Release By J-P Moreau, Paris.   *
     !*                                (www.jpmoreau.fr)            *
     !***************************************************************
     !*******************************************************
     !* Given an array X of N numbers, returns their median *
     !* value XMED. The array X is modified and returned    *
     !* sorted in ascending order.                          *
     !*******************************************************
     SUBROUTINE MDIAN(X,N,XMED)
       real X(N)
       call hpsort(N,X)
       N2=N/2
       if (2*N2.eq.N) then
         XMED = 0.5*(X(N2)+X(N2+1))
       else
         XMED = X(N2+1)
       endif
       return
     END
     
     
     !*****************************************************
     !*  Sorts an array RA of length N in ascending order *
     !*                by the Heapsort method             *
     !* ------------------------------------------------- *
     !* INPUTS:                                           *
     !*	    N	  size of table RA                       *
     !*      RA	  table to be sorted                     *
     !* OUTPUT:                                           *
     !*	    RA    table sorted in ascending order        *
     !*                                                   *
     !* NOTE: The Heapsort method is a N Log2 N routine,  *
     !*       and can be used for very large arrays.      *
     !*****************************************************         
     SUBROUTINE HPSORT(N,RA)
       real RA(N)
       L=N/2+1
       IR=N
       !The index L will be decremented from its initial value during the
       !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
       !will be decremented from its initial value down to 1 during the
       !"retirement-and-promotion" (heap selection) phase.
     10 continue
       if(L > 1)then
         L=L-1
         RRA=RA(L)
       else
         RRA=RA(IR)
         RA(IR)=RA(1)
         IR=IR-1
         if(IR.eq.1)then
           RA(1)=RRA
           return
         end if
       end if
       I=L
       J=L+L
     20 if(J.le.IR)then
       if(J < IR)then
         if(RA(J) < RA(J+1))  J=J+1
       end if
       if(RRA < RA(J))then
         RA(I)=RA(J)
         I=J; J=J+J
       else
         J=IR+1
       end if
       goto 20
       end if
       RA(I)=RRA
       goto 10
     END
     
     !write table of size N to standard output
     SUBROUTINE TWRIT(N,ARR)
     real ARR(N)
       print *,' '
       WRITE(*,10) (ARR(I),I=1,N)
       return
     10 FORMAT(10F6.1)
     END
     
     !end of file tmdian.f90
     
