program generateInputGeometryFiles
	implicit double precision (a-h,o-z)
        PARAMETER (NN=90000)
        DIMENSION LNODE(NN,4),NC1(NN),NC2(NN),NC3(NN),NC4(NN),X(NN),Y(NN)
        open(unit=21,file='tape7.neu')
        open(unit=22,file='grid.dat')
        open(unit=23,file='bctemp.dat')

        read(21,*)nodes,neles
        write(22,*)nodes,neles
		z=0.0d0
		xmax=-10000.0
		xmin=10000.0
		ymax=-10000.0
		ymin=10000.0
        do i = 1, nodes
	        read(21,*)n,x(i),y(i)
	        write(22,100)n,x(i),y(i),z
	        write(27,100)n,x(i),y(i),z
100       	format(i8,3f18.10)
			xmax=amax1(xmax,x(i))
			xmin=amin1(xmin,x(i))
			ymax=amax1(ymax,y(i))
			ymin=amin1(ymin,y(i))
        enddo
		print *,xmax,xmin,ymin,ymax
		print*,nodes,neles

        do i = 1, neles
	        read(21,*)nel,n1,n2,n3,n4
    	    lnode(nel,1)=n1
    	    lnode(nel,2)=n2
    	    lnode(nel,3)=n3
    	    lnode(nel,4)=n4
        enddo

        NGHOST=NELES
        DO I = 1, NELES
	        N1=LNODE(I,1)
    	    N2=LNODE(I,2)
    	    N3=LNODE(I,3)
    	    N4=LNODE(I,4)
            NC1(I)=0
            NC2(I)=0
            NC3(I)=0
            NC4(I)=0
			itest=0
            DO J = 1, NELES
    	        IF(I.NE.J)THEN
                M1=LNODE(J,1)
                M2=LNODE(J,2)
                M3=LNODE(J,3)
                M4=LNODE(J,4)
      if(n1.eq.m2.and.n2.eq.m1.or.n1.eq.m1.and.n2.eq.m2.or.		&
          n1.eq.m3.and.n2.eq.m2.or.n1.eq.m2.and.n2.eq.m3.or.	&
          n1.eq.m4.and.n2.eq.m3.or.n1.eq.m3.and.n2.eq.m4.or.	&
          n1.eq.m1.and.n2.eq.m4.or.n1.eq.m4.and.n2.eq.m1)then	
      nc1(i)=j
      itest=itest+1
      elseif(n2.eq.m2.and.n3.eq.m1.or.n2.eq.m1.and.n3.eq.m2.or.		&
              n2.eq.m3.and.n3.eq.m2.or.n2.eq.m2.and.n3.eq.m3.or.	&
              n2.eq.m4.and.n3.eq.m3.or.n2.eq.m3.and.n3.eq.m4.or.	&
              n2.eq.m1.and.n3.eq.m4.or.n2.eq.m4.and.n3.eq.m1)then
      nc2(i)=j
      itest=itest+1
      elseif(n3.eq.m2.and.n4.eq.m1.or.n3.eq.m1.and.n4.eq.m2.or.		&
              n3.eq.m3.and.n4.eq.m2.or.n3.eq.m2.and.n4.eq.m3.or.	&
              n3.eq.m4.and.n4.eq.m3.or.n3.eq.m3.and.n4.eq.m4.or.	&
              n3.eq.m1.and.n4.eq.m4.or.n3.eq.m4.and.n4.eq.m1)then
      nc3(i)=j
      itest=itest+1
      elseif(n4.eq.m2.and.n1.eq.m1.or.n4.eq.m1.and.n1.eq.m2.or.		&
              n4.eq.m3.and.n1.eq.m2.or.n4.eq.m2.and.n1.eq.m3.or.	&
              n4.eq.m4.and.n1.eq.m3.or.n4.eq.m3.and.n1.eq.m4.or.	&
              n4.eq.m1.and.n1.eq.m4.or.n4.eq.m4.and.n1.eq.m1)then
      nc4(i)=j
      itest=itest+1
      endif
                ENDIF
      if(itest.ge.4)go to 99
                ENDDO
99    continue
          IF(NC1(I).EQ.0)THEN
          NGHOST=NGHOST+1
          NC1(I)=NGHOST
          NTYPE=5
	    NSIDE=1
      	!if(x(n1).lt. 0.00 .and. x(n2) .lt. 0.00)ntype=6
      	if(abs(x(n1)-xmin).lt. 0.0001 .and. abs(x(n2)-xmin) .lt. 0.0001)ntype=2
      	if(abs(x(n1)-xmax).lt. 0.0001 .and. abs(x(n2)-xmax) .lt. 0.0001)ntype=3
       	if(abs(y(n1)-ymax).lt. 0.0001 .and.abs(y(n2)-ymax) .lt. 0.0001)ntype=2
      	if(abs(y(n1)) .lt. 0.0001 .and. abs(y(n2)) .lt. 0.0001)ntype=4
          WRITE(23,*)I ,NGHOST,NTYPE,NSIDE
          ENDIF
          
          IF(NC2(I).EQ.0)THEN
          NGHOST=NGHOST+1
          NC2(I)=NGHOST
          NTYPE=5
	    NSIDE=2
      	!if(x(n2).lt. 0.00 .and. x(n3) .lt. 0.00)ntype=6
      	if(abs(x(n2)-xmin).lt. 0.0001 .and. abs(x(n3)-xmin) .lt. 0.0001)ntype=2
      	if(abs(x(n2)-xmax).lt. 0.0001 .and. abs(x(n3)-xmax) .lt. 0.0001)ntype=3
       	if(abs(y(n2)-ymax).lt. 0.0001 .and.abs(y(n3)-ymax) .lt. 0.0001)ntype=2
      	if(abs(y(n2)) .lt. 0.0001 .and. abs(y(n3)) .lt. 0.0001)ntype=4
          WRITE(23,*)I ,NGHOST,NTYPE,NSIDE
          ENDIF
          
          IF(NC3(I).EQ.0)THEN
          NGHOST=NGHOST+1
          NC3(I)=NGHOST
          NTYPE=5
	    NSIDE=3
      	!if(x(n3).lt. 0.00 .and. x(n3) .lt. 0.00)ntype=6
      	if(abs(x(n3)-xmin).lt. 0.0001 .and. abs(x(n4)-xmin) .lt. 0.0001)ntype=2
      	if(abs(x(n3)-xmax).lt. 0.0001 .and. abs(x(n4)-xmax) .lt. 0.0001)ntype=3
       	if(abs(y(n3)-ymax).lt. 0.0001 .and.abs(y(n4)-ymax) .lt. 0.0001)ntype=2
      	if(abs(y(n3)) .lt. 0.0001 .and. abs(y(n4)) .lt. 0.0001)ntype=4
          WRITE(23,*)I ,NGHOST,NTYPE,NSIDE
          ENDIF
          
          IF(NC4(I).EQ.0)THEN
          NGHOST=NGHOST+1
          NC4(I)=NGHOST
          NTYPE=5
	    NSIDE=4
      	!if(x(n4).lt. 0.00 .and. x(n1) .lt. 0.00)ntype=6
      	if(abs(x(n4)-xmin).lt. 0.0001 .and. abs(x(n1)-xmin) .lt. 0.0001)ntype=2
      	if(abs(x(n4)-xmax).lt. 0.0001 .and. abs(x(n1)-xmax) .lt. 0.0001)ntype=3
       	if(abs(y(n4)-ymax).lt. 0.0001 .and.abs(y(n1)-ymax) .lt. 0.0001)ntype=2
      	if(abs(y(n4)) .lt. 0.0001 .and. abs(y(n1)) .lt. 0.0001)ntype=4
          WRITE(23,*)I ,NGHOST,NTYPE,NSIDE
          ENDIF
         
        WRITE(22,200)I,N1,N2,N3,N4,NC1(I),NC2(I),NC3(I),NC4(I)
        ENDDO
200       FORMAT(9I8)
        close(21)
        close(22)
		close(23)
		
		
		
		open(26,file='bctemp.dat')
		open(25,file='bc.in')
		write(25,*)nghost-neles
		print*,nghost-neles
        do i=1,nghost-neles
        	read(26,*) indx,ng,nt,ns
        	write(25,*) indx,ng,nt,ns
        enddo
        close(26)
        close(25)
        close(27)
        CALL REMOVE(file='bctemp.dat')
endprogram generateInputGeometryFiles
