MODULE AXITURB

	IMPLICIT none
	INTEGER, PARAMETER :: NN=60000,NB=6000
	CHARACTER*15, PARAMETER :: FNAMEGRID='grid1.dat',FNAMEBC='bc.in'
	CHARACTER*50 :: CHOICE,tempchar
	
	INTEGER	:: NODES,NELES,NGHOSTS,LOCAL,NITER,ISTART,IAXIS,IVISCOUS, &
	NCC1,NCC2,NCC3,NCC4,N1,N2,N3,N4,NG,NP,NT,NS,NCOPP,I,IE,ITER,K,J,L,&
	N,NITEROLD,NOLD,nel,npseudo,num,na,nbb
	
	DOUBLE PRECISION :: ALPHA1,ALPHA2,ALPHA3,ALPHA4,ALPHA,BETA,GAMA,  &
	AMOL,PRANDTL,PRANDTLT,CFL,P0,T0,PAMB,AMACHINIT,CP,CV,RR,RUNIV,S0, &
	CONST1,CONST2,CONST3,DIV23,C02,C04,RTOLERANCE,ANR,DNR,DUX,DUY,DVX,&
	DVY,DIVVEL,CMU,TXX,TXY,TYY,QXX,QYY,DQX,DQY,DRX,DRY,AMU,AK,AMUINIT,&
	RHOINIT,VINIT,UINIT,TINIT,EINIT,RINIT,AMUT,AKT,AMUEFF,AKEFF,SX,SY,&
	XA,XB,XC,XD,YA,YB,YC,YD,UA,UB,UC,UD,VA,VB,VC,VD,QA,QB,QC,QD,RA,RB,&
	RC,RD,CHI,TTT,RHO0,DYNRATIO,RHO,E,DTMIN,VEL,DTBYVOL,PHI,E2NC1, 	  &
	E2NC2,E2NC3,E2NC4,E4NC1,E4NC2,E4NC3,E4NC4,D2NC1,D2NC2,D2NC3,D2NC4,&
	D24NC1,D24NC2,D24NC3,D24NC4,RI,FR2,FMU,PK,DEE,UAVG,VAVG,TAVG,DL12,&
	DL23,DL34,DL41,XCCFACE1,XCCFACE2,XCCFACE3,XCCFACE4,YCCFACE1,	  &
	YCCFACE2,YCCFACE3,YCCFACE4,RT,ENX,ENY,TFLAME,ANCOMB,ACOMB,CONVKSC,&
	RHOPROP,D21,D22,D23,D2F,D25,D41,D42,D43,D4F,D45,ANUN,ANUD,YAVG,	  &
	VBYY,CDIFF1,CDIFF2,CDIFF3,CDIFF4,CDIFF5,PNN,TNN,AMN,UNN,VNN,RON,  &
	DTX,DTY,TA,TB,TC,TD,C5I,PI,D4NC1,D4NC2,D4NC3,D4NC4,X1,X2,X3,X4,Y1,&
	Y2,Y3,Y4,PINIT,THRUST,UI,VI,TI,RHOI,UNG,RNG,TNG,VNG,DFLOATITEST,  &
	RAVG,RHOAVG,RHOB,RHOD,RHONG,VOLBYDT,ANUI,ANURTI,PIE,DELTABL,PINF, &
	TINF,AMINF,TWALL,F1I,F2I,F3I,F4I,F5I,G1I,G2I,G3I,G4I,G5I,SX1,SX2, &
	SX3,SX4,SY1,SY2,SY3,SY4, CFLTURB,RESNN,CSOUND,RHOA,RHOC,RTRED,    &
	BURNRATE,C1I,C2I,C3I,C4I,D2Q1I,D2Q2I,D2Q3I,D2Q4I,D2Q5I,rdot,gth,g,&
	rdoteb,Gport,dhyd,DELY
		
	DOUBLE PRECISION, DIMENSION(NN) :: X,Y,XCC,YCC,VOL,DL,DT,SC1X,    &
	SC2X,SC3X,SC4X,SC1Y,SC2Y,SC3Y,SC4Y,C1,C2,C3,C4,C5,F1,F2,F3,F4,F5, &
	G1,G2,G3,G4,G5,R1,R2,R3,R4,R5,COLD1,COLD2,COLD3,COLD4,COLD5,D2Q1, &
	D2Q2,D2Q3,D2Q4,D2Q5,D241,D242,D243,D244,D245,U,V,P,T,VISU,VISV,   &
	VIST,VISR,UN,VN,TN,RN,RHON,SOURCE3,SOURCE5,ANU,ANURT
	
	INTEGER,DIMENSION(NN,4) :: NOD
	
	INTEGER, DIMENSION(NN) :: NC1,NC2,NC3,NC4,ITEST,NPARENT,NGHOST,	  &
	NTYPE,NSIDE,NWALL,NWALLNODES
	
	CONTAINS
	
	SUBROUTINE 	GEOMETRY
		IMPLICIT NONE
		NC1=-1;	NC2=-1;	NC3=-1;	NC4=-1
		OPEN(21, FILE=FNAMEGRID)
        READ(21, *) NODES, NELES
        DO I = 1, NODES
        READ(21, *) N, X(N), Y(N)
	    X(N) = X(N) * 0.001D0
	    Y(N) = Y(N) * 0.001D0
        END DO
        ITEST = 0
        DO I = 1, NELES
            READ(21, *)NEL,NOD(I,1),NOD(I,2),NOD(I,3),NOD(I,4),NC1(I),&
            NC2(I),NC3(I),NC4(I)
        END DO
       ! do i=1,10*neles
        !	tempchar='boundaryoff'
        !	if(i.le.neles) tempchar='compdomain'
        !	print*,i,nc1(i),tempchar
        !enddo
        !pause
        CLOSE(21)
        OPEN(22, FILE=FNAMEBC)
        READ(22, *)NGHOSTS
        DO I = 1, NGHOSTS
            READ(22, *)NPARENT(I), NGHOST(I), NTYPE(I), NSIDE(I)
        END DO
        CLOSE(22)
        DO I = 1, NELES
            N1 = NOD(I, 1)
            N2 = NOD(I, 2)
            N3 = NOD(I, 3)
            N4 = NOD(I, 4)            
            ITEST(N1) = ITEST(N1) +1
            ITEST(N2) = ITEST(N2) +1
            ITEST(N3) = ITEST(N3) +1
            ITEST(N4) = ITEST(N4) +1
            X1 = X(N1)
            X2 = X(N2)
            X3 = X(N3)
            X4 = X(N4)            
            Y1 = Y(N1)
            Y2 = Y(N2)
            Y3 = Y(N3)
            Y4 = Y(N4)
                        
            XCC(I) = (X1 + X2 + X3 + X4)*0.25D0
            YCC(I) = (Y1 + Y2 + Y3 + Y4)*0.25D0            
            
            DL12 = DSQRT((X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2))
            DL23 = DSQRT((X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3))
            DL34 = DSQRT((X3-X4)*(X3-X4) + (Y3-Y4)*(Y3-Y4))
            DL41 = DSQRT((X4-X1)*(X4-X1) + (Y4-Y1)*(Y4-Y1))           
            DL(I) = DMIN1(DL12, DL23, DL34, DL41)
            
            SC1X(I) = (Y2-Y1)
            SC1Y(I) = -(X2-X1)
            SC2X(I) = (Y3-Y2)
            SC2Y(I) = -(X3-X2)
            SC3X(I) = (Y4-Y3)
            SC3Y(I) = -(X4-X3)
            SC4X(I) = (Y1-Y4)
            SC4Y(I) = -(X1-X4)

            VOL(I) = 0.5D0*((X3-X1)*(Y4-Y2)-(X4-X2)*(Y3-Y1))

            IF (VOL(I) .LE. 0.0D0) THEN
              VOL(I) = -VOL(I)
            END IF
            
        END DO
        DO I = 1, NGHOSTS
            NP = NPARENT(I)
            NG = NGHOST(I)
            NS = NSIDE(I)            
            N1 = NOD(NP, 1)
            N2 = NOD(NP, 2)
            N3 = NOD(NP, 3)
            N4 = NOD(NP, 4)            
            X1 = X(N1)
            X2 = X(N2)
            X3 = X(N3)
            X4 = X(N4)            
            Y1 = Y(N1)
            Y2 = Y(N2)
            Y3 = Y(N3)
            Y4 = Y(N4)
            
            IF (NS .EQ. 1) THEN
               XCCFACE1 = 0.50D0*(X1 + X2)
               YCCFACE1 = 0.50D0*(Y1 + Y2)               
               XCC(NG) = 2.0D0*XCCFACE1 - XCC(NP)
               YCC(NG) = 2.0D0*YCCFACE1 - YCC(NP)               
               ITEST(N1) = ITEST(N1) + 1
               ITEST(N2) = ITEST(N2) + 1
            
         ELSE IF (NS .EQ. 2) THEN
               XCCFACE2 = 0.50D0*(X2 + X3)
               YCCFACE2 = 0.50D0*(Y2 + Y3)               
               XCC(NG) = 2.0D0*XCCFACE2 - XCC(NP)
               YCC(NG) = 2.0D0*YCCFACE2 - YCC(NP)               
               ITEST(N2) = ITEST(N2) + 1
               ITEST(N3) = ITEST(N3) + 1
               
         ELSE IF (NS .EQ. 3) THEN
               XCCFACE3 = 0.50D0*(X3 + X4)
               YCCFACE3 = 0.50D0*(Y3 + Y4)               
               XCC(NG) = 2.0D0*XCCFACE3 - XCC(NP)
               YCC(NG) = 2.0D0*YCCFACE3 - YCC(NP)               
               ITEST(N3) = ITEST(N3) + 1
               ITEST(N4) = ITEST(N4) + 1
               
          ELSE IF (NS .EQ. 4) THEN
               XCCFACE4 = 0.50D0*(X4 + X1)
               YCCFACE4 = 0.50D0*(Y4 + Y1)               
               XCC(NG) = 2.0D0*XCCFACE4 - XCC(NP)
               YCC(NG) = 2.0D0*YCCFACE4- YCC(NP)               
               ITEST(N1) = ITEST(N1) + 1
               ITEST(N4) = ITEST(N4) + 1
               
          END IF
          
          VOL(NG) = VOL(NP)
          
          NC1(NG)=NG
          NC2(NG)=NG
          NC3(NG)=NG
          NC4(NG)=NG
          
        END DO
       RETURN
		
	ENDSUBROUTINE GEOMETRY

	
	SUBROUTINE 	DISSIPATION
		IMPLICIT none
		C02=0.25D0
    	C04=0.00390625D0	!1/256		
		
    	DO I=1,NGHOSTS
    		NG=NGHOST(I)
    		NP=NPARENT(I)
      		P(NG)=P(NP)
      		DT(NG)=DT(NP)
      		C1(NG)=C1(NP)
      		C2(NG)=C2(NP)
      		C3(NG)=C3(NP)
      		C4(NG)=C4(NP)
      		C5(NG)=C5(NP)
    	ENDDO		

    	DO I=1,NELES
    		NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
    		PI=P(I)
    		C5I=C5(i)
    		ANUN=DABS(P(NCC1)-PI)+DABS(P(NCC2)-PI)+	&
          		 DABS(P(NCC3)-PI)+DABS(P(NCC4)-PI)
      		ANUD=DABS(P(NCC1)+PI)+DABS(P(NCC2)+PI)+	&
            	 DABS(P(NCC3)+PI)+DABS(P(NCC4)+PI)
		
			ANU(I)= ANUN/ANUD
       	
       		!ANUN=DABS(C5(NCC1)-C5I)+DABS(C5(NCC2)-C5I)+	&
          	!	 DABS(C5(NCC3)-C5I)+DABS(C5(NCC4)-C5I)
      		!ANUD=DABS(C5(NCC1)+C5I)+DABS(C5(NCC2)+C5I)+	&
            !	 DABS(C5(NCC3)+C5I)+DABS(C5(NCC4)+C5I)
            
            !ANUrt(I)= anurt(i)+ ANUN/(ANUD + 0.000001D0)
            
       	
			D2Q1(I)=C1(NCC1)+C1(NCC2)+C1(NCC3)+C1(NCC4)-4.0D0*C1(I)
      		D2Q2(I)=C2(NCC1)+C2(NCC2)+C2(NCC3)+C2(NCC4)-4.0D0*C2(I)
      		D2Q3(I)=C3(NCC1)+C3(NCC2)+C3(NCC3)+C3(NCC4)-4.0D0*C3(I)
      		D2Q4(I)=C4(NCC1)+C4(NCC2)+C4(NCC3)+C4(NCC4)-4.0D0*C4(I)
      		D2Q5(I)=C5(NCC1)+C5(NCC2)+C5(NCC3)+C5(NCC4)-4.0D0*C5(I)
    	ENDDO   

   		DO I=1,NGHOSTS
   			NG=NGHOST(I)
			NP=NPARENT(I)
      		ANU(NG)=ANU(NP)
      		D2Q1(NG)=D2Q1(NP)
      		D2Q2(NG)=D2Q2(NP)
      		D2Q3(NG)=D2Q3(NP)
      		D2Q4(NG)=D2Q4(NP)
      		D2Q5(NG)=D2Q5(NP)
    	ENDDO
	
		DO I=1,NELES
			VOLBYDT=VOL(I)/DT(I)
			ANUI=ANU(I);	ANURTI=ANURT(I)
			C1I=C1(I);	C2I=C2(I);	C3I=C3(I);	C4I=C4(I);	C5I=C5(I)
			D2Q1I=D2Q1(I);	D2Q2I=D2Q2(I);	D2Q3I=D2Q3(I);	
			D2Q4I=D2Q4(I);	D2Q5I=D2Q5(I)
			
			NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
			
	        E2NC1=C02*DMAX1(ANUI,ANU(NCC1))
	        E2NC2=C02*DMAX1(ANUI,ANU(NCC2))
	        E2NC3=C02*DMAX1(ANUI,ANU(NCC3))
	        E2NC4=C02*DMAX1(ANUI,ANU(NCC4))
	
	        E4NC1=DMAX1(0.0D0,C04-E2NC1)
	        E4NC2=DMAX1(0.0D0,C04-E2NC2)
	        E4NC3=DMAX1(0.0D0,C04-E2NC3)
	        E4NC4=DMAX1(0.0D0,C04-E2NC4)
			
	        D24NC1=0.5D0*(VOLBYDT+VOL(NCC1)/DT(NCC1))
	        D24NC2=0.5D0*(VOLBYDT+VOL(NCC2)/DT(NCC2))
	        D24NC3=0.5D0*(VOLBYDT+VOL(NCC3)/DT(NCC3))
	        D24NC4=0.5D0*(VOLBYDT+VOL(NCC4)/DT(NCC4))
	
	        D2NC1=E2NC1*D24NC1
	        D2NC2=E2NC2*D24NC2        
        	D2NC3=E2NC3*D24NC3
        	D2NC4=E2NC4*D24NC4
        
			D4NC1=E4NC1*D24NC1
        	D4NC2=E4NC2*D24NC2        
        	D4NC3=E4NC3*D24NC3
        	D4NC4=E4NC4*D24NC4
					    
			D21=D2NC1*(C1(NCC1)-C1I)+		&
       	    	D2NC2*(C1(NCC2)-C1I)+			&
       	    	D2NC3*(C1(NCC3)-C1I)+			&
      	    	D2NC4*(C1(NCC4)-C1I)      
			D41=D4NC1*(D2Q1(NCC1)-D2Q1I)+	&
      	 		D4NC2*(D2Q1(NCC2)-D2Q1I)+		&
       	    	D4NC3*(D2Q1(NCC3)-D2Q1I)+		&
       	    	D4NC4*(D2Q1(NCC4)-D2Q1I)
    		D241(I)=D21-D41
	
        	D22=D2NC1*(C2(NCC1)-C2I)+		&
       	 		D2NC2*(C2(NCC2)-C2I)+			&
       	    	D2NC3*(C2(NCC3)-C2I)+			&
       	    	D2NC4*(C2(NCC4)-C2I)			
			D42=D4NC1*(D2Q2(NCC1)-D2Q2I)+	&
       	 		D4NC2*(D2Q2(NCC2)-D2Q2I)+		&
       	    	D4NC3*(D2Q2(NCC3)-D2Q2I)+		&
       	    	D4NC4*(D2Q2(NCC4)-D2Q2I)	    		
			D242(I)=D22-D42
	
        	D23=D2NC1*(C3(NCC1)-C3I)+		&
       	 		D2NC2*(C3(NCC2)-C3I)+			&
       	    	D2NC3*(C3(NCC3)-C3I)+			&
       	    	D2NC4*(C3(NCC4)-C3I)       
			D43=D4NC1*(D2Q3(NCC1)-D2Q3I)+	&
       	  		D4NC2*(D2Q3(NCC2)-D2Q3I)+		&
      	    	D4NC3*(D2Q3(NCC3)-D2Q3I)+		&
       	    	D4NC4*(D2Q3(NCC4)-D2Q3I)    			 
        	D243(I)=D23-D43
	
        	D2F=D2NC1*(C4(NCC1)-C4I)+		&
       	   		D2NC2*(C4(NCC2)-C4I)+			&
       	    	D2NC3*(C4(NCC3)-C4I)+			&
       	    	D2NC4*(C4(NCC4)-C4I)
			D4F=D4NC1*(D2Q4(NCC1)-D2Q4I)+	&
       	   		D4NC2*(D2Q4(NCC2)-D2Q4I)+		&
       	    	D4NC3*(D2Q4(NCC3)-D2Q4I)+		&
       	    	D4NC4*(D2Q4(NCC4)-D2Q4I)    			 
        	 D244(I)=D2F-D4F 
        	 
        	!E2NC1=C02*DMAX1(ANUrt(I),ANUrt(NCC1))
	        !E2NC2=C02*DMAX1(ANUrt(I),ANUrt(NCC2))
	        !E2NC3=C02*DMAX1(ANUrt(I),ANUrt(NCC3))
	        !E2NC4=C02*DMAX1(ANUrt(I),ANUrt(NCC4))
	
	        !E4NC1=DMAX1(0.0D0,C04-E2NC1)
	        !E4NC2=DMAX1(0.0D0,C04-E2NC2)
	        !E4NC3=DMAX1(0.0D0,C04-E2NC3)
	        !E4NC4=DMAX1(0.0D0,C04-E2NC4)
	
	        !D24NC1=0.5D0*(VOLBYDT+VOL(NCC1)/DT(NCC1))
	        !D24NC2=0.5D0*(VOLBYDT+VOL(NCC2)/DT(NCC2))
	        !D24NC3=0.5D0*(VOLBYDT+VOL(NCC3)/DT(NCC3))
	        !D24NC4=0.5D0*(VOLBYDT+VOL(NCC4)/DT(NCC4))
	
	        !D2NC1=E2NC1*D24NC1
	        !D2NC2=E2NC2*D24NC2        
        	!D2NC3=E2NC3*D24NC3
        	!D2NC4=E2NC4*D24NC4
        
			!D4NC1=E4NC1*D24NC1
        	!D4NC2=E4NC2*D24NC2        
        	!D4NC3=E4NC3*D24NC3
        	!D4NC4=E4NC4*D24NC4 
        	 
        	D25=D2NC1*(C5(NCC1)-C5I)+		&
       	   		D2NC2*(C5(NCC2)-C5I)+			&
       	    	D2NC3*(C5(NCC3)-C5I)+			&
       	    	D2NC4*(C5(NCC4)-C5I)
			D45=D4NC1*(D2Q5(NCC1)-D2Q5I)+	&
       	   		D4NC2*(D2Q5(NCC2)-D2Q5I)+		&
       	    	D4NC3*(D2Q5(NCC3)-D2Q5I)+		&
       	    	D4NC4*(D2Q5(NCC4)-D2Q5I)   			 
        	D245(I)=D25-D45    
  		ENDDO  
  		
  		RETURN		
	ENDSUBROUTINE DISSIPATION
	
	SUBROUTINE VISCOUS
		IMPLICIT NONE
		DIV23=2.0D0/3.0D0
		ALPHA=0.00052D0
		BETA=0.2D0
		CONST1=21.479D0
		CONST2=21.394D0
		CONST3=1.5D0
		CMU=0.09D0
		
		UN=0.0D0
		VN=0.0D0
		TN=0.0D0
		RN=0.0D0
		RHON=0.0D0
		
		VISU=0.0D0
		VISV=0.0D0
		VIST=0.0D0
		VISR=0.0D0
		ITEST=0
		
		DO I =1,NELES
			N1=NOD(I,1)
			N2=NOD(I,2)
			N3=NOD(I,3)
			N4=NOD(I,4)
			
			UI=U(I)
			VI=V(I)
			TI=T(I)
			RHOI=C1(I)
			RI=C5(I)/RHOI
			
			UN(N1)=UN(N1)+UI
			VN(N1)=VN(N1)+VI
			TN(N1)=TN(N1)+TI
			RHON(N1)=RHON(N1)+RHOI
			RN(N1)=RN(N1)+RI
			
			UN(N2)=UN(N2)+UI
			VN(N2)=VN(N2)+VI
			TN(N2)=TN(N2)+TI
			RHON(N2)=RHON(N2)+RHOI
			RN(N2)=RN(N2)+RI
			
			UN(N3)=UN(N3)+UI
			VN(N3)=VN(N3)+VI
			TN(N3)=TN(N3)+TI
			RHON(N3)=RHON(N3)+RHOI
			RN(N3)=RN(N3)+RI
			
			UN(N4)=UN(N4)+UI
			VN(N4)=VN(N4)+VI
			TN(N4)=TN(N4)+TI
			RHON(N4)=RHON(N4)+RHOI
			RN(N4)=RN(N4)+RI
			
			ITEST(N1)=ITEST(N1)+1
			ITEST(N2)=ITEST(N2)+1
			ITEST(N3)=ITEST(N3)+1
			ITEST(N4)=ITEST(N4)+1			
		ENDDO
		
		DO I=1,NGHOSTS
			NP=NPARENT(I)
			NG=NGHOST(I)
			NT=NTYPE(I)
			NS=NSIDE(I)
			
			UNG=U(NG)
			VNG=V(NG)
			TNG=T(NG)
			RHONG=C1(NG)
			RNG=C5(NG)/RHONG
			
			IF(NS==1)THEN
				N1=NOD(NP,1)
				N2=NOD(NP,2)
			ELSEIF(NS==2)THEN
				N1=NOD(NP,2)
				N2=NOD(NP,3)
			ELSEIF(NS==3)THEN
				N1=NOD(NP,3)
				N2=NOD(NP,4)
			ELSEIF(NS==4)THEN
				N1=NOD(NP,4)
				N2=NOD(NP,1)
			ENDIF
			
			UN(N1)=UN(N1)+UNG
			VN(N1)=VN(N1)+VNG
			TN(N1)=TN(N1)+TNG
			RHON(N1)=RHON(N1)+RHONG
			RN(N1)=RN(N1)+RNG
			
			UN(N2)=UN(N2)+UNG
			VN(N2)=VN(N2)+VNG
			TN(N2)=TN(N2)+TNG
			RHON(N2)=RHON(N2)+RHONG
			RN(N2)=RN(N2)+RNG
			
			IF(NT==4)THEN
				VN(N1)=0.0D0
				VN(N2)=0.0D0
			ELSEIF(NT==5)THEN
				UN(N1)=0.0D0
				UN(N2)=0.0D0
				VN(N1)=0.0D0
				VN(N2)=0.0D0
				RN(N1)=RTOLERANCE*RTRED
				RN(N2)=RTOLERANCE*RTRED
			ENDIF
			
			ITEST(N1)=ITEST(N1)+1
			ITEST(N2)=ITEST(N2)+1	
		ENDDO
		
		DO I=1,NODES
			DFLOATITEST=DFLOAT(ITEST(I))
			UN(I)=UN(I)/DFLOATITEST
			VN(I)=VN(I)/DFLOATITEST
			TN(I)=TN(I)/DFLOATITEST
			RHON(I)=RHON(I)/DFLOATITEST
			RN(I)=RN(I)/DFLOATITEST
			IF(RN(I) .LT. RTOLERANCE) RN(I)=RTOLERANCE
		ENDDO
		
		DO I=1,NELES
			N1=NOD(I,1)
			N2=NOD(I,2)
			N3=NOD(I,3)
			N4=NOD(I,4)
			
			NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
			
			!FACE1
			SX=SC1X(I)
			SY=SC1Y(I)
			
			XA=XCC(I);		YA=YCC(I);
			XB=X(N1);		YB=Y(N1);
			XC=XCC(NCC1);	YC=YCC(NCC1);
			XD=X(N2);		YD=Y(N2);
			
			UA=U(I);	VA=V(I);	TA=T(I)
			UB=UN(N1);	VB=VN(N1);	TB=TN(N1)
			UC=U(NCC1);	VC=V(NCC1);	TC=T(NCC1)
			UD=UN(N2);	VD=VN(N2);	TD=TN(N2)
			
			RHOA=C1(I);		RA=C5(I)/C1(I)
			RHOB=RHON(N1);	RB=RN(N1)
			RHOC=C1(NCC1);	RC=C5(NCC1)/C1(NCC1)
			RHOD=RHON(N2);	RD=RN(N2) 
			
			DNR = (XA-XC)*(YB-YD)-(XB-XD)*(YA-YC)
                
			ANR = (UA-UC)*(YB-YD)-(UB-UD)*(YA-YC)
			DUX = ANR/DNR
                
            ANR = (VA-VC)*(YB-YD)-(VB-VD)*(YA-YC)
			DVX = ANR/DNR
                
            ANR = (TA-TC)*(YB-YD)-(TB-TD)*(YA-YC)
			DTX = ANR/DNR
			
			ANR = (RA-RC)*(YB-YD)-(RB-RD)*(YA-YC)
			DRX = ANR/DNR
				
            ANR = -(UA-UC)*(XB-XD)+(UB-UD)*(XA-XC)
			DUY = ANR/DNR	
                
            ANR = -(VA-VC)*(XB-XD)+(VB-VD)*(XA-XC)
			DVY = ANR/DNR
                
            ANR = -(TA-TC)*(XB-XD)+(TB-TD)*(XA-XC)
			DTY = ANR/DNR
			
			ANR = -(RA-RC)*(XB-XD)+(RB-RD)*(XA-XC)
			DRY = ANR/DNR
			
			UAVG = 0.25D0*(UA+UB+UC+UD)
			VAVG = 0.25D0*(VA+VB+VC+VD)
			TAVG = 0.25D0*(TA+TB+TC+TD)
			RHOAVG = 0.25D0*(RHOA+RHOB+RHOC+RHOD)
			RAVG = 0.25D0*(RA+RB+RC+RD)
			
			DIVVEL = DUX + DVY
			
			amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(TAVG**0.6D0)
            AK = AMU*CP/PRANDTL
            
            CHI = RHOAVG*RAVG/(AMU*CMU)
            FMU = DTANH(ALPHA*CHI*CHI)/DTANH(BETA*CHI*CHI)
            AMUT = RHOAVG*RAVG*FMU
            AKT = AMUT*CP/PRANDTLT
            
            AMUEFF = AMU + AMUT
            AKEFF = AK + AKT
            
            TXX = AMUEFF*(DUX + DUX - DIV23*DIVVEL)
            TYY = AMUEFF*(DVY + DVY - DIV23*DIVVEL)
            TXY = AMUEFF*(DUY + DVX)
            
            QXX = UAVG*TXX + VAVG*TXY + AKEFF*DTX
            QYY = UAVG*TXY + VAVG*TYY + AKEFF*DTY
            
            VISU(I) = TXX*SX + TXY*SY
            VISV(I) = TXY*SX + TYY*SY
            VIST(I) = QXX*SX + QYY*SY
            VISR(I) = AMUEFF*DRX*SX + AMUEFF*DRY*SY
            
            !FACE2
            SX=SC2X(I)
			SY=SC2Y(I)
			
			XA=XCC(I);		YA=YCC(I);
			XB=X(N2);		YB=Y(N2);
			XC=XCC(NCC2);	YC=YCC(NCC2);
			XD=X(N3);		YD=Y(N3);
			
			UA=U(I);	VA=V(I);	TA=T(I)
			UB=UN(N2);	VB=VN(N2);	TB=TN(N2)
			UC=U(NCC2);	VC=V(NCC2);	TC=T(NCC2)
			UD=UN(N3);	VD=VN(N3);	TD=TN(N3)
			
			RHOA=C1(I);		RA=C5(I)/C1(I)
			RHOB=RHON(N2);	RB=RN(N2)
			RHOC=C1(NCC2);	RC=C5(NCC2)/C1(NCC2)
			RHOD=RHON(N3);	RD=RN(N3) 
			
			DNR = (XA-XC)*(YB-YD)-(XB-XD)*(YA-YC)
                
			ANR = (UA-UC)*(YB-YD)-(UB-UD)*(YA-YC)
			DUX = ANR/DNR
                
            ANR = (VA-VC)*(YB-YD)-(VB-VD)*(YA-YC)
			DVX = ANR/DNR
                
            ANR = (TA-TC)*(YB-YD)-(TB-TD)*(YA-YC)
			DTX = ANR/DNR
			
			ANR = (RA-RC)*(YB-YD)-(RB-RD)*(YA-YC)
			DRX = ANR/DNR
				
            ANR = -(UA-UC)*(XB-XD)+(UB-UD)*(XA-XC)
			DUY = ANR/DNR	
                
            ANR = -(VA-VC)*(XB-XD)+(VB-VD)*(XA-XC)
			DVY = ANR/DNR
                
            ANR = -(TA-TC)*(XB-XD)+(TB-TD)*(XA-XC)
			DTY = ANR/DNR
			
			ANR = -(RA-RC)*(XB-XD)+(RB-RD)*(XA-XC)
			DRY = ANR/DNR
			
			UAVG = 0.25D0*(UA+UB+UC+UD)
			VAVG = 0.25D0*(VA+VB+VC+VD)
			TAVG = 0.25D0*(TA+TB+TC+TD)
			RHOAVG = 0.25D0*(RHOA+RHOB+RHOC+RHOD)
			RAVG = 0.25D0*(RA+RB+RC+RD)
			
			DIVVEL = DUX + DVY
			
			amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(TAVG**0.6D0)
            AK = AMU*CP/PRANDTL
            
            CHI = RHOAVG*RAVG/(AMU*CMU)
            FMU = DTANH(ALPHA*CHI*CHI)/DTANH(BETA*CHI*CHI)
            AMUT = RHOAVG*RAVG*FMU
            AKT = AMUT*CP/PRANDTLT
            
            AMUEFF = AMU + AMUT
            AKEFF = AK + AKT
            
            TXX = AMUEFF*(DUX + DUX - DIV23*DIVVEL)
            TYY = AMUEFF*(DVY + DVY - DIV23*DIVVEL)
            TXY = AMUEFF*(DUY + DVX)
            
            QXX = UAVG*TXX + VAVG*TXY + AKEFF*DTX
            QYY = UAVG*TXY + VAVG*TYY + AKEFF*DTY
            
            VISU(I) = VISU(I) + TXX*SX + TXY*SY
            VISV(I) = VISV(I) + TXY*SX + TYY*SY
            VIST(I) = VIST(I) + QXX*SX + QYY*SY
            VISR(I) = VISR(I) + AMUEFF*DRX*SX + AMUEFF*DRY*SY
            
            !FACE3
            SX=SC3X(I)
			SY=SC3Y(I)
			
			XA=XCC(I);		YA=YCC(I);
			XB=X(N3);		YB=Y(N3);
			XC=XCC(NCC3);	YC=YCC(NCC3);
			XD=X(N4);		YD=Y(N4);
			
			UA=U(I);	VA=V(I);	TA=T(I)
			UB=UN(N3);	VB=VN(N3);	TB=TN(N3)
			UC=U(NCC3);	VC=V(NCC3);	TC=T(NCC3)
			UD=UN(N4);	VD=VN(N4);	TD=TN(N4)
			
			RHOA=C1(I);		RA=C5(I)/C1(I)
			RHOB=RHON(N3);	RB=RN(N3)
			RHOC=C1(NCC3);	RC=C5(NCC3)/C1(NCC3)
			RHOD=RHON(N4);	RD=RN(N4) 
			
			DNR = (XA-XC)*(YB-YD)-(XB-XD)*(YA-YC)
                
			ANR = (UA-UC)*(YB-YD)-(UB-UD)*(YA-YC)
			DUX = ANR/DNR
                
            ANR = (VA-VC)*(YB-YD)-(VB-VD)*(YA-YC)
			DVX = ANR/DNR
                
            ANR = (TA-TC)*(YB-YD)-(TB-TD)*(YA-YC)
			DTX = ANR/DNR
			
			ANR = (RA-RC)*(YB-YD)-(RB-RD)*(YA-YC)
			DRX = ANR/DNR
				
            ANR = -(UA-UC)*(XB-XD)+(UB-UD)*(XA-XC)
			DUY = ANR/DNR	
                
            ANR = -(VA-VC)*(XB-XD)+(VB-VD)*(XA-XC)
			DVY = ANR/DNR
                
            ANR = -(TA-TC)*(XB-XD)+(TB-TD)*(XA-XC)
			DTY = ANR/DNR
			
			ANR = -(RA-RC)*(XB-XD)+(RB-RD)*(XA-XC)
			DRY = ANR/DNR
			
			UAVG = 0.25D0*(UA+UB+UC+UD)
			VAVG = 0.25D0*(VA+VB+VC+VD)
			TAVG = 0.25D0*(TA+TB+TC+TD)
			RHOAVG = 0.25D0*(RHOA+RHOB+RHOC+RHOD)
			RAVG = 0.25D0*(RA+RB+RC+RD)
			
			DIVVEL = DUX + DVY
			
			amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(TAVG**0.6D0)
            AK = AMU*CP/PRANDTL
            
            CHI = RHOAVG*RAVG/(AMU*CMU)
            FMU = DTANH(ALPHA*CHI*CHI)/DTANH(BETA*CHI*CHI)
            AMUT = RHOAVG*RAVG*FMU
            AKT = AMUT*CP/PRANDTLT
            
            AMUEFF = AMU + AMUT
            AKEFF = AK + AKT
            
            TXX = AMUEFF*(DUX + DUX - DIV23*DIVVEL)
            TYY = AMUEFF*(DVY + DVY - DIV23*DIVVEL)
            TXY = AMUEFF*(DUY + DVX)
            
            QXX = UAVG*TXX + VAVG*TXY + AKEFF*DTX
            QYY = UAVG*TXY + VAVG*TYY + AKEFF*DTY
            
            VISU(I) = VISU(I) + TXX*SX + TXY*SY
            VISV(I) = VISV(I) + TXY*SX + TYY*SY
            VIST(I) = VIST(I) + QXX*SX + QYY*SY
            VISR(I) = VISR(I) + AMUEFF*DRX*SX + AMUEFF*DRY*SY
            
            !FACE4
            SX=SC4X(I)
			SY=SC4Y(I)
			
			XA=XCC(I);		YA=YCC(I);
			XB=X(N4);		YB=Y(N4);
			XC=XCC(NCC4);	YC=YCC(NCC4);
			XD=X(N1);		YD=Y(N1);
			
			UA=U(I);	VA=V(I);	TA=T(I)
			UB=UN(N4);	VB=VN(N4);	TB=TN(N4)
			UC=U(NCC4);	VC=V(NCC4);	TC=T(NCC4)
			UD=UN(N1);	VD=VN(N1);	TD=TN(N1)
			
			RHOA=C1(I);		RA=C5(I)/C1(I)
			RHOB=RHON(N4);	RB=RN(N4)
			RHOC=C1(NCC4);	RC=C5(NCC4)/C1(NCC4)
			RHOD=RHON(N1);	RD=RN(N1) 
			
			DNR = (XA-XC)*(YB-YD)-(XB-XD)*(YA-YC)
                
			ANR = (UA-UC)*(YB-YD)-(UB-UD)*(YA-YC)
			DUX = ANR/DNR
                
            ANR = (VA-VC)*(YB-YD)-(VB-VD)*(YA-YC)
			DVX = ANR/DNR
                
            ANR = (TA-TC)*(YB-YD)-(TB-TD)*(YA-YC)
			DTX = ANR/DNR
			
			ANR = (RA-RC)*(YB-YD)-(RB-RD)*(YA-YC)
			DRX = ANR/DNR
				
            ANR = -(UA-UC)*(XB-XD)+(UB-UD)*(XA-XC)
			DUY = ANR/DNR	
                
            ANR = -(VA-VC)*(XB-XD)+(VB-VD)*(XA-XC)
			DVY = ANR/DNR
                
            ANR = -(TA-TC)*(XB-XD)+(TB-TD)*(XA-XC)
			DTY = ANR/DNR
			
			ANR = -(RA-RC)*(XB-XD)+(RB-RD)*(XA-XC)
			DRY = ANR/DNR
			
			UAVG = 0.25D0*(UA+UB+UC+UD)
			VAVG = 0.25D0*(VA+VB+VC+VD)
			TAVG = 0.25D0*(TA+TB+TC+TD)
			RHOAVG = 0.25D0*(RHOA+RHOB+RHOC+RHOD)
			RAVG = 0.25D0*(RA+RB+RC+RD)
			
			DIVVEL = DUX + DVY
			
			amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(TAVG**0.6D0)
            AK = AMU*CP/PRANDTL
            
            CHI = RHOAVG*RAVG/(AMU*CMU)
            FMU = DTANH(ALPHA*CHI*CHI)/DTANH(BETA*CHI*CHI)
            AMUT = RHOAVG*RAVG*FMU
            AKT = AMUT*CP/PRANDTLT
            
            AMUEFF = AMU + AMUT
            AKEFF = AK + AKT
            
            TXX = AMUEFF*(DUX + DUX - DIV23*DIVVEL)
            TYY = AMUEFF*(DVY + DVY - DIV23*DIVVEL)
            TXY = AMUEFF*(DUY + DVX)
            
            QXX = UAVG*TXX + VAVG*TXY + AKEFF*DTX
            QYY = UAVG*TXY + VAVG*TYY + AKEFF*DTY
            
            VISU(I) = VISU(I) + TXX*SX + TXY*SY
            VISV(I) = VISV(I) + TXY*SX + TYY*SY
            VIST(I) = VIST(I) + QXX*SX + QYY*SY
            VISR(I) = VISR(I) + AMUEFF*DRX*SX + AMUEFF*DRY*SY
        	
		ENDDO
		
		RETURN
	ENDSUBROUTINE VISCOUS
	
	SUBROUTINE SOURCETERMS
		IMPLICIT none
		DIV23=2.0D0/3.0D0
		ALPHA=0.00052D0
		BETA=0.2D0
		CONST1=21.479D0
		CONST2=21.394D0
		CONST3=1.5D0
		CMU=0.09D0
		
		SOURCE3=0.0D0
		SOURCE5=0.0D0
		ITEST=0
		
		DO I=1,NELES
			NCC1=NC1(I)		
            NCC2=NC2(I)
            NCC3=NC3(I)
            NCC4=NC4(I)            
            XA=XCC(NCC1);XB=XCC(NCC2);XC=XCC(NCC3);XD=XCC(NCC4)            
            YA=YCC(NCC1);YB=YCC(NCC2);YC=YCC(NCC3);YD=YCC(NCC4)                        
			UA=U(NCC1);UB=U(NCC2);UC=U(NCC3);UD=U(NCC4)			
			VA=V(NCC1);VB=V(NCC2);VC=V(NCC3);VD=V(NCC4)			
			RA=C5(NCC1)/C1(NCC1);RB=C5(NCC2)/C1(NCC2)		
			RC=C5(NCC3)/C1(NCC3);RD=C5(NCC4)/C1(NCC4)
			
			QA=DSQRT(UA*UA+VA*VA)
			QB=DSQRT(UB*UB+VB*VB)
			QC=DSQRT(UC*UC+VC*VC)
			QD=DSQRT(UD*UD+VD*VD)
						
			DNR = (XA-XC)*(YB-YD)-(XB-XD)*(YA-YC)
            
            ANR = (QA-QC)*(YB-YD)-(QB-QD)*(YA-YC)
			DQX = ANR/DNR
			
			ANR = (RA-RC)*(YB-YD)-(RB-RD)*(YA-YC)
			DRX = ANR/DNR
			
			ANR = (UA-UC)*(YB-YD)-(UB-UD)*(YA-YC)
			DUX = ANR/DNR
			
			ANR = (VA-VC)*(YB-YD)-(VB-VD)*(YA-YC)
			DVX = ANR/DNR
			
			ANR = -(QA-QC)*(XB-XD)+(QB-QD)*(XA-XC)
			DQY = ANR/DNR
			
			ANR = -(RA-RC)*(XB-XD)+(RB-RD)*(XA-XC)
			DRY = ANR/DNR
			
			ANR = -(UA-UC)*(XB-XD)+(UB-UD)*(XA-XC)
			DUY = ANR/DNR
			
			ANR = -(VA-VC)*(XB-XD)+(VB-VD)*(XA-XC)
			DVY = ANR/DNR
			
			RHO=C1(I)
			RI=C5(I)/RHO
			AMU=1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(I)**0.6D0)
			CHI=RHO*RI/(CMU*AMU)
			ANR=ALPHA*CHI*CHI
			DNR=BETA*CHI*CHI
            FMU=DTANH(ANR)/DTANH(DNR)
            ANR=DSQRT(CHI)
            DNR=CHI*DSQRT(2.0D0*DSQRT(CMU))
		 	FR2=DTANH(ANR)/DTANH(DNR)	
			AMUT=FMU*RHO*RI
			PHI = DQX*DRX + DQY*DRY
			
			DEE= 0.0D0
			IF(PHI>0.0D0) DEE = DRX*DRX + DRY*DRY
			
			DIVVEL = DUX + DVY
			
			PK =  AMUT/RHO*(2.0*(DUX*DUX+DVY*DVY)+          &
			      (DUY+DVX)*(DUY+DVX)-DIV23*DIVVEL*DIVVEL)
			PK = DMAX1(PK,0.0D0)
			
			SOURCE5(I) = (RHO*(CONST1-CONST2*FR2)*DSQRT(RI*PK) -   	  &
	 		              RHO*CONST3*DEE)*(AMUT/DMAX1(AMUT,0.001D0*AMU))**npseudo
		ENDDO
		RETURN
	ENDSUBROUTINE SOURCETERMS
	
	SUBROUTINE POSTPROCESS
		IMPLICIT NONE
		
		PRINT*, 'Postprocessing.'
		print*,' '

		C1=0.0D0
		C2=0.0D0
		C3=0.0D0
		C4=0.0D0
		C5=0.0D0
		ITEST=0
		
		DO I = 1, NELES
        	DO J = 1, 4
        	  K = NOD(I,J)
        	  ITEST(K) = ITEST(K) + 1
        	  C1(K)=C1(K)+P(I)
        	  C2(K)=C2(K)+U(I)
        	  C3(K)=C3(K)+V(I)
        	  C4(K)=C4(K)+T(I)
        	  c5(k)=c5(k)+r2(i)
        	END DO
      	END DO
      	
      	DO I =1,NGHOSTS
      		NG=NGHOST(I)
      		NP=NPARENT(I)
      		NT=NTYPE(I)
      		NS=NSIDE(I)      		
      		IF(NS==1)THEN
      			N1=NOD(NP,1)
      			N2=NOD(NP,2)
      		ELSEIF(NS==2)THEN
      			N1=NOD(NP,2)
      			N2=NOD(NP,3)
      		ELSEIF(NS==3)THEN
      			N1=NOD(NP,3)
      			N2=NOD(NP,4)      	
      		ELSEIF(NS==4)THEN
      			N1=NOD(NP,4)
      			N2=NOD(NP,1)      	
      		ENDIF      		
      		C1(N1)=C1(N1)+P(NG)
      		C2(N1)=C2(N1)+U(NG)
      		C3(N1)=C3(N1)+V(NG)
      		C4(N1)=C4(N1)+T(NG)
      		c5(n1)=c5(n1)+r2(ng)
      		C1(N2)=C1(N2)+P(NG)
      		C2(N2)=C2(N2)+U(NG)
      		C3(N2)=C3(N2)+V(NG)
      		C4(N2)=C4(N2)+T(NG)
      		c5(n2)=c5(n2)+r2(ng)
      		
      		ITEST(N1)=ITEST(N1)+1
      		ITEST(N1)=ITEST(N1)+1 	
      	ENDDO
      	
      	DO I = 1, NODES
      		DFLOATITEST=DFLOAT(ITEST(I))
          C1(I)=C1(I)/DFLOATITEST
          C2(I)=C2(I)/DFLOATITEST
          C3(I)=C3(I)/DFLOATITEST
          C4(I)=C4(I)/DFLOATITEST
          c5(i)=c5(i)/dfloatitest
      	ENDDO
      	
      	OPEN(UNIT=26,FILE='t_plot.dat')
      	OPEN(28,FILE='pressures_wall.dat')
      	OPEN(35,FILE='contour.dat')
      	open(31,FILE='burnrate.dat')
  		 	
      	WRITE(26,*)'VARIABLES = "X","Y","U","V","P","RHO","MACH","T","RESIDUAL"'
      	WRITE(26,*)'ZONE F=FEPOINT,ET=QUADRILATERAL,N =',		&
     				 NODES,',E =',NELES
      	DO I = 1,NODES
        	PNN=C1(I)
        	UNN=C2(I)
        	VNN=C3(I)
        	TNN=C4(I)
        	resnn=c5(i)
        	AMN=DSQRT(UNN*UNN+VNN*VNN)/DSQRT(GAMA*RR*TNN)
        	RON=PNN/(RR*TNN)   
        	WRITE(26,*)X(I),Y(I),UNN,VNN,PNN*0.00001D0,RON,AMN,TNN,resnn
        	WRITE(35,*)X(I),Y(I),UNN,VNN,PNN*0.00001D0,RON,AMN,TNN,resnn
      	ENDDO
      	DO IE = 1,NELES
        	WRITE(26,*)NOD(IE,1),NOD(IE,2),NOD(IE,3),NOD(IE,4)
      	ENDDO   
      	CLOSE(26)
      	
      	  		
  		CLOSE(28)
  		CLOSE(35)
  		CLOSE(31)
	ENDSUBROUTINE POSTPROCESS

ENDMODULE AXITURB


PROGRAM MAINSOLVER
	USE AXITURB
	IMPLICIT none
		
		OPEN(25,FILE='flow.in')
		READ(25,*) GAMA,AMOL,LOCAL,NITER,ISTART,CFL,PRANDTL,PRANDTLT,P0,T0,PAMB,&
				AMACHINIT,TWALL,PINF,DELTABL,RTRED,NPSEUDO     
		CLOSE(25)!problem specific
		print*,'Data read from:  ',FNAMEGRID,'  and  ',FNAMEBC
		PRINT*,'Refresh convergence data file? (Y/N)'
		READ*,CHOICE
		
		IF(istart==0)THEN
			open(27,file='history.txt')
			OPEN(32,FILE='convergence.dat')
		ELSE 
			open(27,file='history.txt',access='APPEND')
			OPEN(32,FILE='convergence.dat',ACCESS='APPEND')
		ENDIF
		
		CALL GEOMETRY

		gth=35.0d0
		RUNIV=8314.0D0		
		RR=RUNIV/AMOL		
    	CV=RR/(GAMA-1.0D0)
    	CP=CV*GAMA
      	RHO0=P0/(RR*T0)!PROBLEM SPECIFIC
		S0=P0/(RHO0**GAMA)
		TINF=T0/(1.0D0+0.5D0*(GAMA-1.0D0)*AMINF*AMINF)
		VISU=0.0D0
		VISV=0.0D0
		VIST=0.0D0
		VISR=0.0D0
		SOURCE5=0.0D0
		ITEST=0
		CFLTURB=0.5D0
		ALPHA1=0.25D0
		ALPHA2=1.0D0/3.0D0
		ALPHA3=0.5D0
		ALPHA4=1.0D0
		NITEROLD=0
		NOLD=0
		CONVKSC=1.0D0/98066.5D0
		
		IF (ISTART==0)THEN!problem specific
    		DYNRATIO=1.0D0+(GAMA-1.0D0)*0.5D0*AMACHINIT**2.0D0
    		TINIT=T0/DYNRATIO
        	PINIT=P0/(DYNRATIO**(GAMA/(GAMA-1.0D0)))
			UINIT=AMACHINIT*DSQRT(GAMA*RR*TINIT)
        	VINIT=0.0D0
			EINIT=CV*TINIT+0.5D0*(UINIT**2.0D0+VINIT**2.0D0)
        	RHO0=P0/(RR*T0)
			S0=P0/(RHO0**GAMA)
        	RHOINIT=(PINIT/(RR*TINIT))
        	AMUINIT = 1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(TINIT**0.6D0)
        	RINIT=AMUINIT/RHOINIT        	
        	RTOLERANCE=RINIT*RTRED
        	niterold=0
        	nold=0
		
			DO I=1,NELES        	
				C1(I)=RHOINIT
            	C2(I)=RHOINIT*UINIT
            	C3(I)=RHOINIT*VINIT
            	C4(I)=RHOINIT*EINIT
            	C5(I)=RHOINIT*RINIT
    		ENDDO
    	
    	ELSEIF (ISTART==1)THEN
    		nold=1
    		OPEN(55,FILE='restart.in')
        	DO I=1,NELES
				READ(55,*) C1(I),C2(I),C3(I),C4(I),C5(I)
        	ENDDO 
        		read(55,*) RTRED
        		read(55,*) npseudo
				READ(55,*) RINIT
				read(55,*) NITEROLD
				
            CLOSE(55)  
            RTOLERANCE=RINIT*RTRED
		ENDIF
		
		DO I=1,NELES    
				COLD1(I)=C1(I)
				COLD2(I)=C2(I)
				COLD3(I)=C3(I)
				COLD4(I)=C4(I)
				COLD5(I)=C5(I)
		ENDDO
		
		!History file
		call date_and_time(date,time)
		if(istart==0)then
			write(27,*) "New simulation: Viscous-Turbulence solver (Two-D)"
			
		else
			write(27,*) 'Continued simulation from ',niterold
		endif
			write(27,*) 'On    ',date,'  at    ',time
			write(27,*) "Using"
			write(27,*)"Grid file: ",fnamegrid
			write(27,*) "Boundary conditions file: ",fnamebc
			write(27,*) "Parameters:"		
			write(27,*)'GAMA,AMOL,LOCAL,NITER,ISTART,CFL,PRANDTL,PRANDTLT,P0,T0,PAMB,&
				AMACHINIT,TWALL,PINF,DELTABL,RTRED,NPSEUDO '
			write(27,2007)GAMA,AMOL,LOCAL,NITER,ISTART,CFL,PRANDTL,PRANDTLT,P0,T0,PAMB,&
				AMACHINIT,TWALL,PINF,DELTABL,RTRED,NPSEUDO 
		print *,'GAMA,AMOL,LOCAL,NITER,ISTART,CFL,PRANDTL,PRANDTLT,P0,T0,PAMB,&
				AMACHINIT,TWALL,PINF,DELTABL,RTRED,NPSEUDO '
100		print 2007,GAMA,AMOL,LOCAL,NITER,ISTART,CFL,PRANDTL,PRANDTLT,P0,T0,PAMB,&
				AMACHINIT,TWALL,PINF,DELTABL,RTRED,NPSEUDO 
	
2007	format(2f8.2,3i7,10f15.3)
		
		DO ITER=1,NITER
			DTMIN=100000.0D0
			
			
    	 	DO I=1,NELES        	
				RHO=C1(I)
        		U(I)=C2(I)/RHO
        		V(I)=C3(I)/RHO
				VEL=DSQRT(U(I)*U(I)+V(I)*V(I))
				E=C4(I)/RHO
				T(I)=(E-0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
				P(I)=RHO*RR*T(I)
				DT(I)=CFL*DL(I)/(DABS(VEL)+DSQRT(GAMA*RR*T(I)))
				DTMIN=DMIN1(DTMIN,DT(I))
			ENDDO
			
			IF (LOCAL==0)THEN
    			DO I=1,NELES
					DT(I)=DTMIN
    			ENDDO
			ENDIF		
			
			CALL DISSIPATION
        	
        DO I=1,NELES
            RHO=C1(I)
            U(I)=C2(I)/RHO
        	V(I)=C3(I)/RHO
           	E=C4(I)/RHO
			T(I)=(E-0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
            P(I)=RHO*RR*T(I)
            F1(I)=C2(I)
            F2(I)=RHO*U(I)*U(I)+P(I)
            F3(I)=RHO*U(I)*V(I)
            F4(I)=(C4(I)+P(I))*U(I)
            G1(I)=C3(I)
            G2(I)=F3(I)
            G3(I)=RHO*V(I)*V(I)+P(I)
            G4(I)=(C4(I)+P(I))*V(I)
            RT=C5(I)/RHO
            F5(I)=RHO*RT*U(I)
            G5(I)=RHO*RT*V(I)
        ENDDO
        
        DO I=1,NGHOSTS
			NP=NPARENT(I)
            NG=NGHOST(I)
            NT=NTYPE(I)
            NS=NSIDE(I)            
			IF (NT==2)THEN
				U(NG)=DABS(U(NP))
                V(NG)=-V(NP)
                T(NG)=T0-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
                RHO=(S0/(RR*T(NG)))**(1.0D0/(1.0D0-GAMA))                
				P(NG)=RHO*RR*T(NG)
				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=RHO*RT
			   	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)			
			ELSEIF (NT==3)THEN
				IF(NS==1)THEN
					NCOPP=NC3(NP)
				ELSEIF(NS==2)THEN
					NCOPP=NC4(NP)
				ELSEIF(NS==3)THEN
					NCOPP=NC1(NP)
				ELSEIF(NS==4)THEN
					NCOPP=NC2(NP)
				ENDIF	
				U(NG)=DABS(2.0D0*U(NP)-U(NCOPP))
                V(NG)=(2.0D0*V(NP)-V(NCOPP))
                T(NG)=(2.0D0*T(NP)-T(NCOPP))
				P(NG)=(2.0D0*P(NP)-P(NCOPP))
				IF(PAMB>=0.0D0)THEN
                	P(NG)=PAMB
                ENDIF
				RHO=P(NG)/(RR*T(NG))
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))
                C4(NG)=RHO*E		
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=G5(NP)
            ELSEIF (NT==4)THEN		
				U(NG)=U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=-G5(NP)
            ELSEIF (NT==5)THEN
				U(NG)=-U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
                F1(NG)=-F1(NP)
                F2(NG)=2.0D0*P(NP)-F2(NP)
                F3(NG)=-F3(NP)
                F4(NG)=-F4(NP)
                G1(NG)=-G1(NP)
                G2(NG)=-G2(NP)
                G3(NG)=2.0D0*P(NP)-G3(NP)
                G4(NG)=-G4(NP)
                C5(NG)=-C5(NP)
                F5(NG)=-F5(NP)
                G5(NG)=-G5(NP)
			ELSEIF(NT==6)THEN
				DHYD=DABS(YCC(NG)+YCC(NP))
				
				IF(NS==1)THEN
					DNR=DSQRT(SC1X(NP)*SC1X(NP)+SC1Y(NP)*SC1Y(NP))
					ENX=-SC1X(NP)/DNR
					ENY=-SC1Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC3(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC3(N)
					DO WHILE (NC3(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC3(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==2)THEN
					DNR=DSQRT(SC2X(NP)*SC2X(NP)+SC2Y(NP)*SC2Y(NP))
					ENX=-SC2X(NP)/DNR
					ENY=-SC2Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC4(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC4(N)
					DO WHILE (NC4(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC4(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==3)THEN
					DNR=DSQRT(SC3X(NP)*SC3X(NP)+SC3Y(NP)*SC3Y(NP))
					ENX=-SC3X(NP)/DNR
					ENY=-SC3Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC1(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC1(N)
					DO WHILE (NC1(N) .ne. n)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC1(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==4)THEN
					DNR=DSQRT(SC4X(NP)*SC4X(NP)+SC4Y(NP)*SC4Y(NP))
					ENX=-SC4X(NP)/DNR
					ENY=-SC4Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N3)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC2(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC2(N)
					DO WHILE (NC2(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC2(N)
					ENDDO
					GPORT=GPORT/DHYD
				ENDIF
				
				RHO=C1(NP)
				RDOT=ACOMB*((P(NP)*CONVKSC)**(ANCOMB))
				RDOTEB=RDOT
				VEL=RHOPROP/RHO*RDOTEB
				P(NG)=P(NP)
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				AMU=1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				G=(GPORT/(RHOPROP*RDOT))*((RHOPROP*RDOT*DHYD/AMU)**(-0.125D0))
				IF(G .GT. GTH)THEN
					RDOTEB=RDOTEB*(1.0D0+0.023D0*((G**0.8D0)-(GTH**0.8D0)))
				ENDIF
				VEL=RHOPROP/RHO*RDOTEB
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				RHO=P(NG)/(RR*T(NG))
				AMU = 1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C5(NG)=RHO*RT
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)
            ELSEIF(NT==7)THEN!PROBLEMSPECIFIC
	   				P(NG)=PINF
	   				T(NG)=TINF
	   				V(NG)=-V(NP)
	   				RHO=P(NG)/(RR*T(NG))
	   				IF(YCC(NG).LE.DELTABL)THEN
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF)*		  &
	   						  (YCC(NG)/DELTABL)**(1.0D0/7.0D0) - U(NP)
	   				ELSE
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF) - U(NP)
	   				ENDIF
	   				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
					RT= AMU/RHO
					C1(NG)=RHO
					C2(NG)=RHO*U(NG)
					C3(NG)=RHO*V(NG)
					E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
					C4(NG)=RHO*E
            		F1(NG)=C2(NG)
            		F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            		F3(NG)=RHO*U(NG)*V(NG)
            		F4(NG)=(C4(NG)+P(NG))*U(NG)
            		G1(NG)=C3(NG)
            		G2(NG)=F3(NG)
            		G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            		G4(NG)=(C4(NG)+P(NG))*V(NG)
            		C5(NG)=RHO*RT
					F5(NG)=C5(NG)*U(NG)
            		G5(NG)=C5(NG)*V(NG)
            	ENDIF
	   		ENDDO
	   	
			CALL VISCOUS 
			
			CALL SOURCETERMS
			
        DO I=1,NELES
        	NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
			SX1=SC1X(I);SX2=SC2X(I);SX3=SC3X(I);SX4=SC4X(I)
			SY1=SC1Y(I);SY2=SC2Y(I);SY3=SC3Y(I);SY4=SC4Y(I)
			F1I=F1(I);F2I=F2(I);F3I=F3(I);F4I=F4(I);F5I=F5(I)
			G1I=G1(I);G2I=G2(I);G3I=G3(I);G4I=G4(I);G5I=G5(I)
			
    	   R1(I)=0.5D0*((F1I+F1(NCC1))*SX1+&
       					(G1I+G1(NCC1))*SY1+&
       		 			(F1I+F1(NCC2))*SX2+&
     					(G1I+G1(NCC2))*SY2+&
      					(F1I+F1(NCC3))*SX3+& 
     					(G1I+G1(NCC3))*SY3+&
          				(F1I+F1(NCC4))*SX4+& 
     					(G1I+G1(NCC4))*SY4)
    	   R2(I)=0.5D0*((F2I+F2(NCC1))*SX1+&
        				(G2I+G2(NCC1))*SY1+& 
     		 			(F2I+F2(NCC2))*SX2+& 
     					(G2I+G2(NCC2))*SY2+& 
           				(F2I+F2(NCC3))*SX3+& 
     					(G2I+G2(NCC3))*SY3+& 
           				(F2I+F2(NCC4))*SX4+& 
     					(G2I+G2(NCC4))*SY4)
           R3(I)=0.5D0*((F3I+F3(NCC1))*SX1+&
       					(G3I+G3(NCC1))*SY1+& 
         			 	(F3I+F3(NCC2))*SX2+& 
     					(G3I+G3(NCC2))*SY2+& 
         			 	(F3I+F3(NCC3))*SX3+& 
     					(G3I+G3(NCC3))*SY3+& 
         			 	(F3I+F3(NCC4))*SX4+&
     				 	(G3I+G3(NCC4))*SY4)
    	   R4(I)=0.5D0*((F4I+F4(NCC1))*SX1+&
        			 	(G4I+G4(NCC1))*SY1+&           	    
     	 				(F4I+F4(NCC2))*SX2+& 
     					(G4I+G4(NCC2))*SY2+& 
        	 			(F4I+F4(NCC3))*SX3+& 
     					(G4I+G4(NCC3))*SY3+& 
          				(F4I+F4(NCC4))*SX4+& 
     					(G4I+G4(NCC4))*SY4)
     	   R5(I)=0.5D0*((F5I+F5(NCC1))*SX1+&
           		 		(G5I+G5(NCC1))*SY1+&           	    
     	   				(F5I+F5(NCC2))*SX2+& 
     	   				(G5I+G5(NCC2))*SY2+& 
           				(F5I+F5(NCC3))*SX3+& 
     	      			(G5I+G5(NCC3))*SY3+& 
           				(F5I+F5(NCC4))*SX4+& 
    	   				(G5I+G5(NCC4))*SY4)-VOL(I)*SOURCE5(I)
		ENDDO		   	
		
	   	DO I=1,NELES
	   		DTBYVOL=DT(I)/VOL(I)
          C1(I)=COLD1(I)-ALPHA1*DTBYVOL*(R1(I)-D241(I))
          C2(I)=COLD2(I)-ALPHA1*DTBYVOL*(R2(I)-VISU(I)-D242(I))
          C3(I)=COLD3(I)-ALPHA1*DTBYVOL*(R3(I)-VISV(I)-D243(I))
          C4(I)=COLD4(I)-ALPHA1*DTBYVOL*(R4(I)-VIST(I)-D244(I))
  	  	  C5(I)=COLD5(I)-ALPHA1*DTBYVOL*CFLTURB*(R5(I)-VISR(I)-D245(I))
          IF(C5(I)/C1(I).le.RTOLERANCE)   C5(I)=C1(I)*RTOLERANCE
       	ENDDO

		DO I=1,NELES
            RHO=C1(I)
            U(I)=C2(I)/RHO
        	V(I)=C3(I)/RHO
           	E=C4(I)/RHO
			T(I)=(E-0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
            P(I)=RHO*RR*T(I)
            F1(I)=C2(I)
            F2(I)=RHO*U(I)*U(I)+P(I)
            F3(I)=RHO*U(I)*V(I)
            F4(I)=(C4(I)+P(I))*U(I)
            G1(I)=C3(I)
            G2(I)=F3(I)
            G3(I)=RHO*V(I)*V(I)+P(I)
            G4(I)=(C4(I)+P(I))*V(I)
            RT=C5(I)/RHO
            F5(I)=RHO*RT*U(I)
            G5(I)=RHO*RT*V(I)
        ENDDO
        
        
        
        DO I=1,NGHOSTS
			NP=NPARENT(I)
            NG=NGHOST(I)
            NT=NTYPE(I)
            NS=NSIDE(I)            
			IF (NT==2)THEN
				U(NG)=DABS(U(NP))
                V(NG)=-V(NP)
                T(NG)=T0-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
                RHO=(S0/(RR*T(NG)))**(1.0D0/(1.0D0-GAMA))                
				P(NG)=RHO*RR*T(NG)
				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=RHO*RT
			   	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)			
			ELSEIF (NT==3)THEN
				IF(NS==1)THEN
					NCOPP=NC3(NP)
				ELSEIF(NS==2)THEN
					NCOPP=NC4(NP)
				ELSEIF(NS==3)THEN
					NCOPP=NC1(NP)
				ELSEIF(NS==4)THEN
					NCOPP=NC2(NP)
				ENDIF	
				U(NG)=DABS(2.0D0*U(NP)-U(NCOPP))
                V(NG)=(2.0D0*V(NP)-V(NCOPP))
                T(NG)=(2.0D0*T(NP)-T(NCOPP))
				P(NG)=(2.0D0*P(NP)-P(NCOPP))
				IF(PAMB>=0.0D0)THEN
                	P(NG)=PAMB
                ENDIF
				RHO=P(NG)/(RR*T(NG))
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))
                C4(NG)=RHO*E		
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=G5(NP)
            ELSEIF (NT==4)THEN		
				U(NG)=U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=-G5(NP)
            ELSEIF (NT==5)THEN
				U(NG)=-U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
                F1(NG)=-F1(NP)
                F2(NG)=2.0D0*P(NP)-F2(NP)
                F3(NG)=-F3(NP)
                F4(NG)=-F4(NP)
                G1(NG)=-G1(NP)
                G2(NG)=-G2(NP)
                G3(NG)=2.0D0*P(NP)-G3(NP)
                G4(NG)=-G4(NP)
                C5(NG)=-C5(NP)
                F5(NG)=-F5(NP)
                G5(NG)=-G5(NP)
            ELSEIF(NT==6)THEN
				DHYD=DABS(YCC(NG)+YCC(NP))
				
				IF(NS==1)THEN
					DNR=DSQRT(SC1X(NP)*SC1X(NP)+SC1Y(NP)*SC1Y(NP))
					ENX=-SC1X(NP)/DNR
					ENY=-SC1Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC3(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC3(N)
					DO WHILE (NC3(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC3(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==2)THEN
					DNR=DSQRT(SC2X(NP)*SC2X(NP)+SC2Y(NP)*SC2Y(NP))
					ENX=-SC2X(NP)/DNR
					ENY=-SC2Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC4(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC4(N)
					DO WHILE (NC4(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC4(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==3)THEN
					DNR=DSQRT(SC3X(NP)*SC3X(NP)+SC3Y(NP)*SC3Y(NP))
					ENX=-SC3X(NP)/DNR
					ENY=-SC3Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC1(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC1(N)
					DO WHILE (NC1(N) .ne. n)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC1(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==4)THEN
					DNR=DSQRT(SC4X(NP)*SC4X(NP)+SC4Y(NP)*SC4Y(NP))
					ENX=-SC4X(NP)/DNR
					ENY=-SC4Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N3)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC2(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC2(N)
					DO WHILE (NC2(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC2(N)
					ENDDO
					GPORT=GPORT/DHYD
				ENDIF
				
				RHO=C1(NP)
				RDOT=ACOMB*((P(NP)*CONVKSC)**(ANCOMB))
				RDOTEB=RDOT
				VEL=RHOPROP/RHO*RDOTEB
				P(NG)=P(NP)
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				AMU=1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				G=(GPORT/(RHOPROP*RDOT))*((RHOPROP*RDOT*DHYD/AMU)**(-0.125D0))
				IF(G .GT. GTH)THEN
					RDOTEB=RDOTEB*(1.0D0+0.023D0*((G**0.8D0)-(GTH**0.8D0)))
				ENDIF
				VEL=RHOPROP/RHO*RDOTEB
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				RHO=P(NG)/(RR*T(NG))
				AMU = 1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C5(NG)=RHO*RT
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)
            ELSEIF(NT==7)THEN!PROBLEMSPECIFIC
	   				P(NG)=PINF
	   				T(NG)=TINF
	   				V(NG)=-V(NP)
	   				RHO=P(NG)/(RR*T(NG))
	   				IF(YCC(NG).LE.DELTABL)THEN
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF)*		  &
	   						  (YCC(NG)/DELTABL)**(1.0D0/7.0D0) - U(NP)
	   				ELSE
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF) - U(NP)
	   				ENDIF
	   				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
					RT= AMU/RHO
					C1(NG)=RHO
					C2(NG)=RHO*U(NG)
					C3(NG)=RHO*V(NG)
					E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
					C4(NG)=RHO*E
            		F1(NG)=C2(NG)
            		F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            		F3(NG)=RHO*U(NG)*V(NG)
            		F4(NG)=(C4(NG)+P(NG))*U(NG)
            		G1(NG)=C3(NG)
            		G2(NG)=F3(NG)
            		G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            		G4(NG)=(C4(NG)+P(NG))*V(NG)
            		C5(NG)=RHO*RT
					F5(NG)=C5(NG)*U(NG)
            		G5(NG)=C5(NG)*V(NG)
            	ENDIF
	   		ENDDO
	   		
        	CALL SOURCETERMS

        DO I=1,NELES
        	NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
			SX1=SC1X(I);SX2=SC2X(I);SX3=SC3X(I);SX4=SC4X(I)
			SY1=SC1Y(I);SY2=SC2Y(I);SY3=SC3Y(I);SY4=SC4Y(I)
			F1I=F1(I);F2I=F2(I);F3I=F3(I);F4I=F4(I);F5I=F5(I)
			G1I=G1(I);G2I=G2(I);G3I=G3(I);G4I=G4(I);G5I=G5(I)
			
    	   R1(I)=0.5D0*((F1I+F1(NCC1))*SX1+&
       					(G1I+G1(NCC1))*SY1+&
       		 			(F1I+F1(NCC2))*SX2+&
     					(G1I+G1(NCC2))*SY2+&
      					(F1I+F1(NCC3))*SX3+& 
     					(G1I+G1(NCC3))*SY3+&
          				(F1I+F1(NCC4))*SX4+& 
     					(G1I+G1(NCC4))*SY4)
    	   R2(I)=0.5D0*((F2I+F2(NCC1))*SX1+&
        				(G2I+G2(NCC1))*SY1+& 
     		 			(F2I+F2(NCC2))*SX2+& 
     					(G2I+G2(NCC2))*SY2+& 
           				(F2I+F2(NCC3))*SX3+& 
     					(G2I+G2(NCC3))*SY3+& 
           				(F2I+F2(NCC4))*SX4+& 
     					(G2I+G2(NCC4))*SY4)
           R3(I)=0.5D0*((F3I+F3(NCC1))*SX1+&
       					(G3I+G3(NCC1))*SY1+& 
         			 	(F3I+F3(NCC2))*SX2+& 
     					(G3I+G3(NCC2))*SY2+& 
         			 	(F3I+F3(NCC3))*SX3+& 
     					(G3I+G3(NCC3))*SY3+& 
         			 	(F3I+F3(NCC4))*SX4+&
     				 	(G3I+G3(NCC4))*SY4)
    	   R4(I)=0.5D0*((F4I+F4(NCC1))*SX1+&
        			 	(G4I+G4(NCC1))*SY1+&           	    
     	 				(F4I+F4(NCC2))*SX2+& 
     					(G4I+G4(NCC2))*SY2+& 
        	 			(F4I+F4(NCC3))*SX3+& 
     					(G4I+G4(NCC3))*SY3+& 
          				(F4I+F4(NCC4))*SX4+& 
     					(G4I+G4(NCC4))*SY4)
     	   R5(I)=0.5D0*((F5I+F5(NCC1))*SX1+&
           		 		(G5I+G5(NCC1))*SY1+&           	    
     	   				(F5I+F5(NCC2))*SX2+& 
     	   				(G5I+G5(NCC2))*SY2+& 
           				(F5I+F5(NCC3))*SX3+& 
     	      			(G5I+G5(NCC3))*SY3+& 
           				(F5I+F5(NCC4))*SX4+& 
    	   				(G5I+G5(NCC4))*SY4)-VOL(I)*SOURCE5(I)
		ENDDO	
	   	
	   	
	   	
	   	DO I=1,NELES
	   		DTBYVOL=DT(I)/VOL(I)
          C1(I)=COLD1(I)-ALPHA2*DTBYVOL*(R1(I)-D241(I))
          C2(I)=COLD2(I)-ALPHA2*DTBYVOL*(R2(I)-VISU(I)-D242(I))
          C3(I)=COLD3(I)-ALPHA2*DTBYVOL*(R3(I)-VISV(I)-D243(I))
          C4(I)=COLD4(I)-ALPHA2*DTBYVOL*(R4(I)-VIST(I)-D244(I))
  	  	  C5(I)=COLD5(I)-ALPHA2*DTBYVOL*CFLTURB*(R5(I)-VISR(I)-D245(I))
          IF(C5(I)/C1(I).le.RTOLERANCE)	C5(I)=C1(I)*RTOLERANCE
       	ENDDO
        
        DO I=1,NELES
            RHO=C1(I)
            U(I)=C2(I)/RHO
        	V(I)=C3(I)/RHO
           	E=C4(I)/RHO
			T(I)=(E-0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
            P(I)=RHO*RR*T(I)
            F1(I)=C2(I)
            F2(I)=RHO*U(I)*U(I)+P(I)
            F3(I)=RHO*U(I)*V(I)
            F4(I)=(C4(I)+P(I))*U(I)
            G1(I)=C3(I)
            G2(I)=F3(I)
            G3(I)=RHO*V(I)*V(I)+P(I)
            G4(I)=(C4(I)+P(I))*V(I)
            RT=C5(I)/RHO
            F5(I)=RHO*RT*U(I)
            G5(I)=RHO*RT*V(I)
        ENDDO
        
        DO I=1,NGHOSTS
			NP=NPARENT(I)
            NG=NGHOST(I)
            NT=NTYPE(I)
            NS=NSIDE(I)            
			IF (NT==2)THEN
				U(NG)=DABS(U(NP))
                V(NG)=-V(NP)
                T(NG)=T0-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
                RHO=(S0/(RR*T(NG)))**(1.0D0/(1.0D0-GAMA))                
				P(NG)=RHO*RR*T(NG)
				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=RHO*RT
			   	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)			
			ELSEIF (NT==3)THEN
				IF(NS==1)THEN
					NCOPP=NC3(NP)
				ELSEIF(NS==2)THEN
					NCOPP=NC4(NP)
				ELSEIF(NS==3)THEN
					NCOPP=NC1(NP)
				ELSEIF(NS==4)THEN
					NCOPP=NC2(NP)
				ENDIF	
				U(NG)=DABS(2.0D0*U(NP)-U(NCOPP))
                V(NG)=(2.0D0*V(NP)-V(NCOPP))
                T(NG)=(2.0D0*T(NP)-T(NCOPP))
				P(NG)=(2.0D0*P(NP)-P(NCOPP))
				IF(PAMB>=0.0D0)THEN
                	P(NG)=PAMB
                ENDIF
				RHO=P(NG)/(RR*T(NG))
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))
                C4(NG)=RHO*E		
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=G5(NP)
            ELSEIF (NT==4)THEN		
				U(NG)=U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=-G5(NP)
            ELSEIF (NT==5)THEN
				U(NG)=-U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
                F1(NG)=-F1(NP)
                F2(NG)=2.0D0*P(NP)-F2(NP)
                F3(NG)=-F3(NP)
                F4(NG)=-F4(NP)
                G1(NG)=-G1(NP)
                G2(NG)=-G2(NP)
                G3(NG)=2.0D0*P(NP)-G3(NP)
                G4(NG)=-G4(NP)
                C5(NG)=-C5(NP)
                F5(NG)=-F5(NP)
                G5(NG)=-G5(NP)
            ELSEIF(NT==6)THEN
				DHYD=DABS(YCC(NG)+YCC(NP))
				
				IF(NS==1)THEN
					DNR=DSQRT(SC1X(NP)*SC1X(NP)+SC1Y(NP)*SC1Y(NP))
					ENX=-SC1X(NP)/DNR
					ENY=-SC1Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC3(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC3(N)
					DO WHILE (NC3(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC3(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==2)THEN
					DNR=DSQRT(SC2X(NP)*SC2X(NP)+SC2Y(NP)*SC2Y(NP))
					ENX=-SC2X(NP)/DNR
					ENY=-SC2Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC4(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC4(N)
					DO WHILE (NC4(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC4(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==3)THEN
					DNR=DSQRT(SC3X(NP)*SC3X(NP)+SC3Y(NP)*SC3Y(NP))
					ENX=-SC3X(NP)/DNR
					ENY=-SC3Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC1(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC1(N)
					DO WHILE (NC1(N) .ne. n)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC1(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==4)THEN
					DNR=DSQRT(SC4X(NP)*SC4X(NP)+SC4Y(NP)*SC4Y(NP))
					ENX=-SC4X(NP)/DNR
					ENY=-SC4Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N3)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC2(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC2(N)
					DO WHILE (NC2(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC2(N)
					ENDDO
					GPORT=GPORT/DHYD
				ENDIF
				
				RHO=C1(NP)
				RDOT=ACOMB*((P(NP)*CONVKSC)**(ANCOMB))
				RDOTEB=RDOT
				VEL=RHOPROP/RHO*RDOTEB
				P(NG)=P(NP)
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				AMU=1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				G=(GPORT/(RHOPROP*RDOT))*((RHOPROP*RDOT*DHYD/AMU)**(-0.125D0))
				IF(G .GT. GTH)THEN
					RDOTEB=RDOTEB*(1.0D0+0.023D0*((G**0.8D0)-(GTH**0.8D0)))
				ENDIF
				VEL=RHOPROP/RHO*RDOTEB
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				RHO=P(NG)/(RR*T(NG))
				AMU = 1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C5(NG)=RHO*RT
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)
            ELSEIF(NT==7)THEN!PROBLEMSPECIFIC
	   				P(NG)=PINF
	   				T(NG)=TINF
	   				V(NG)=-V(NP)
	   				RHO=P(NG)/(RR*T(NG))
	   				IF(YCC(NG).LE.DELTABL)THEN
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF)*		  &
	   						  (YCC(NG)/DELTABL)**(1.0D0/7.0D0) - U(NP)
	   				ELSE
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF) - U(NP)
	   				ENDIF
	   				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
					RT= AMU/RHO
					C1(NG)=RHO
					C2(NG)=RHO*U(NG)
					C3(NG)=RHO*V(NG)
					E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
					C4(NG)=RHO*E
            		F1(NG)=C2(NG)
            		F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            		F3(NG)=RHO*U(NG)*V(NG)
            		F4(NG)=(C4(NG)+P(NG))*U(NG)
            		G1(NG)=C3(NG)
            		G2(NG)=F3(NG)
            		G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            		G4(NG)=(C4(NG)+P(NG))*V(NG)
            		C5(NG)=RHO*RT
					F5(NG)=C5(NG)*U(NG)
            		G5(NG)=C5(NG)*V(NG)
            	ENDIF
	   		ENDDO
	   		
        	CALL SOURCETERMS
        
        DO I=1,NELES
        	NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
			SX1=SC1X(I);SX2=SC2X(I);SX3=SC3X(I);SX4=SC4X(I)
			SY1=SC1Y(I);SY2=SC2Y(I);SY3=SC3Y(I);SY4=SC4Y(I)
			F1I=F1(I);F2I=F2(I);F3I=F3(I);F4I=F4(I);F5I=F5(I)
			G1I=G1(I);G2I=G2(I);G3I=G3(I);G4I=G4(I);G5I=G5(I)
			
    	   R1(I)=0.5D0*((F1I+F1(NCC1))*SX1+&
       					(G1I+G1(NCC1))*SY1+&
       		 			(F1I+F1(NCC2))*SX2+&
     					(G1I+G1(NCC2))*SY2+&
      					(F1I+F1(NCC3))*SX3+& 
     					(G1I+G1(NCC3))*SY3+&
          				(F1I+F1(NCC4))*SX4+& 
     					(G1I+G1(NCC4))*SY4)
    	   R2(I)=0.5D0*((F2I+F2(NCC1))*SX1+&
        				(G2I+G2(NCC1))*SY1+& 
     		 			(F2I+F2(NCC2))*SX2+& 
     					(G2I+G2(NCC2))*SY2+& 
           				(F2I+F2(NCC3))*SX3+& 
     					(G2I+G2(NCC3))*SY3+& 
           				(F2I+F2(NCC4))*SX4+& 
     					(G2I+G2(NCC4))*SY4)
           R3(I)=0.5D0*((F3I+F3(NCC1))*SX1+&
       					(G3I+G3(NCC1))*SY1+& 
         			 	(F3I+F3(NCC2))*SX2+& 
     					(G3I+G3(NCC2))*SY2+& 
         			 	(F3I+F3(NCC3))*SX3+& 
     					(G3I+G3(NCC3))*SY3+& 
         			 	(F3I+F3(NCC4))*SX4+&
     				 	(G3I+G3(NCC4))*SY4)-VOL(I)*SOURCE3(I)
    	   R4(I)=0.5D0*((F4I+F4(NCC1))*SX1+&
        			 	(G4I+G4(NCC1))*SY1+&           	    
     	 				(F4I+F4(NCC2))*SX2+& 
     					(G4I+G4(NCC2))*SY2+& 
        	 			(F4I+F4(NCC3))*SX3+& 
     					(G4I+G4(NCC3))*SY3+& 
          				(F4I+F4(NCC4))*SX4+& 
     					(G4I+G4(NCC4))*SY4)
     	   R5(I)=0.5D0*((F5I+F5(NCC1))*SX1+&
           		 		(G5I+G5(NCC1))*SY1+&           	    
     	   				(F5I+F5(NCC2))*SX2+& 
     	   				(G5I+G5(NCC2))*SY2+& 
           				(F5I+F5(NCC3))*SX3+& 
     	      			(G5I+G5(NCC3))*SY3+& 
           				(F5I+F5(NCC4))*SX4+& 
    	   				(G5I+G5(NCC4))*SY4)-VOL(I)*SOURCE5(I)
		ENDDO	
	   	
	   	DO I=1,NELES
	   		DTBYVOL=DT(I)/VOL(I)
          C1(I)=COLD1(I)-ALPHA3*DTBYVOL*(R1(I)-D241(I))
          C2(I)=COLD2(I)-ALPHA3*DTBYVOL*(R2(I)-VISU(I)-D242(I))
          C3(I)=COLD3(I)-ALPHA3*DTBYVOL*(R3(I)-VISV(I)-D243(I))
          C4(I)=COLD4(I)-ALPHA3*DTBYVOL*(R4(I)-VIST(I)-D244(I))
  	  	  C5(I)=COLD5(I)-ALPHA3*DTBYVOL*CFLTURB*(R5(I)-VISR(I)-D245(I))
          IF(C5(I)/C1(I).le.RTOLERANCE)	C5(I)=C1(I)*RTOLERANCE
       	ENDDO
	   	
	   DO I=1,NELES
            RHO=C1(I)
            U(I)=C2(I)/RHO
        	V(I)=C3(I)/RHO
           	E=C4(I)/RHO
			T(I)=(E-0.5D0*(U(I)*U(I)+V(I)*V(I)))/CV
            P(I)=RHO*RR*T(I)
            F1(I)=C2(I)
            F2(I)=RHO*U(I)*U(I)+P(I)
            F3(I)=RHO*U(I)*V(I)
            F4(I)=(C4(I)+P(I))*U(I)
            G1(I)=C3(I)
            G2(I)=F3(I)
            G3(I)=RHO*V(I)*V(I)+P(I)
            G4(I)=(C4(I)+P(I))*V(I)
            RT=C5(I)/RHO
            F5(I)=RHO*RT*U(I)
            G5(I)=RHO*RT*V(I)
        ENDDO
        
        DO I=1,NGHOSTS
			NP=NPARENT(I)
            NG=NGHOST(I)
            NT=NTYPE(I)
            NS=NSIDE(I)            
			IF (NT==2)THEN
				U(NG)=DABS(U(NP))
                V(NG)=-V(NP)
                T(NG)=T0-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP
                RHO=(S0/(RR*T(NG)))**(1.0D0/(1.0D0-GAMA))                
				P(NG)=RHO*RR*T(NG)
				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=RHO*RT
			   	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)			
			ELSEIF (NT==3)THEN
				IF(NS==1)THEN
					NCOPP=NC3(NP)
				ELSEIF(NS==2)THEN
					NCOPP=NC4(NP)
				ELSEIF(NS==3)THEN
					NCOPP=NC1(NP)
				ELSEIF(NS==4)THEN
					NCOPP=NC2(NP)
				ENDIF	
				U(NG)=DABS(2.0D0*U(NP)-U(NCOPP))
                V(NG)=(2.0D0*V(NP)-V(NCOPP))
                T(NG)=(2.0D0*T(NP)-T(NCOPP))
				P(NG)=(2.0D0*P(NP)-P(NCOPP))
				IF(PAMB>=0.0D0)THEN
                	P(NG)=PAMB
                ENDIF
				RHO=P(NG)/(RR*T(NG))
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))
                C4(NG)=RHO*E		
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=G5(NP)
            ELSEIF (NT==4)THEN		
				U(NG)=U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	C5(NG)=C5(NP)
            	F5(NG)=F5(NP)
            	G5(NG)=-G5(NP)
            ELSEIF (NT==5)THEN
				U(NG)=-U(NP)
                V(NG)=-V(NP)
                T(NG)=T(NP)
                RHO=C1(NP)
				P(NG)=P(NP)
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
                F1(NG)=-F1(NP)
                F2(NG)=2.0D0*P(NP)-F2(NP)
                F3(NG)=-F3(NP)
                F4(NG)=-F4(NP)
                G1(NG)=-G1(NP)
                G2(NG)=-G2(NP)
                G3(NG)=2.0D0*P(NP)-G3(NP)
                G4(NG)=-G4(NP)
                C5(NG)=-C5(NP)!QUESTION SHOULD THE FLUX AT THE WALLS BE MADE ZERO?
                F5(NG)=-F5(NP)
                G5(NG)=-G5(NP)
            ELSEIF(NT==6)THEN
				DHYD=DABS(YCC(NG)+YCC(NP))
				
				IF(NS==1)THEN
					DNR=DSQRT(SC1X(NP)*SC1X(NP)+SC1Y(NP)*SC1Y(NP))
					ENX=-SC1X(NP)/DNR
					ENY=-SC1Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC3(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)	
					DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC3(N)
					DO WHILE (NC3(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N1)+Y(N2))-(Y(N3)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC3(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==2)THEN
					DNR=DSQRT(SC2X(NP)*SC2X(NP)+SC2Y(NP)*SC2Y(NP))
					ENX=-SC2X(NP)/DNR
					ENY=-SC2Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC4(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC4(N)
					DO WHILE (NC4(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N2)+Y(N3))-(Y(N1)+Y(N4)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC4(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==3)THEN
					DNR=DSQRT(SC3X(NP)*SC3X(NP)+SC3Y(NP)*SC3Y(NP))
					ENX=-SC3X(NP)/DNR
					ENY=-SC3Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC1(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC1(N)
					DO WHILE (NC1(N) .ne. n)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N3)+Y(N4))-(Y(N1)+Y(N2)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC1(N)
					ENDDO
					GPORT=GPORT/DHYD
				ELSEIF(NS==4)THEN
					DNR=DSQRT(SC4X(NP)*SC4X(NP)+SC4Y(NP)*SC4Y(NP))
					ENX=-SC4X(NP)/DNR
					ENY=-SC4Y(NP)/DNR
					n1=nod(np,1);n2=nod(np,2);n3=nod(np,3);n4=nod(np,4)	
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N3)+Y(N2)))
					GPORT=C1(NP)*U(NP)*DELY
					N=NC2(NP)
					n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
					DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
					GPORT=GPORT+C1(N)*U(N)*DELY
					N=NC2(N)
					DO WHILE (NC2(N) .ne. N)
						n1=nod(n,1);n2=nod(n,2);n3=nod(n,3);n4=nod(n,4)
						DELY=0.5D0*((Y(N4)+Y(N1))-(Y(N2)+Y(N3)))
						GPORT=GPORT+C1(N)*U(N)*DELY
						N=NC2(N)
					ENDDO
					GPORT=GPORT/DHYD
				ENDIF
				
				RHO=C1(NP)
				RDOT=ACOMB*((P(NP)*CONVKSC)**(ANCOMB))
				RDOTEB=RDOT
				VEL=RHOPROP/RHO*RDOTEB
				P(NG)=P(NP)
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				AMU=1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				G=(GPORT/(RHOPROP*RDOT))*((RHOPROP*RDOT*DHYD/AMU)**(-0.125D0))
				IF(G .GT. GTH)THEN
					RDOTEB=RDOTEB*(1.0D0+0.023D0*((G**0.8D0)-(GTH**0.8D0)))
				ENDIF
				VEL=RHOPROP/RHO*RDOTEB
				U(NG)=VEL*ENX+VEL*ENX-U(NP)
				V(NG)=VEL*ENY+VEL*ENY-V(NP)
				T(NG)=TFLAME-0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))/CP 
				RHO=P(NG)/(RR*T(NG))
				AMU = 1.0D0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
				RT= AMU/RHO
				C5(NG)=RHO*RT
				C1(NG)=RHO
				C2(NG)=RHO*U(NG)
				C3(NG)=RHO*V(NG)
				E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
				C4(NG)=RHO*E
            	F1(NG)=C2(NG)
            	F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            	F3(NG)=RHO*U(NG)*V(NG)
            	F4(NG)=(C4(NG)+P(NG))*U(NG)
            	G1(NG)=C3(NG)
            	G2(NG)=F3(NG)
            	G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            	G4(NG)=(C4(NG)+P(NG))*V(NG)
            	F5(NG)=C5(NG)*U(NG)
            	G5(NG)=C5(NG)*V(NG)
            ELSEIF(NT==7)THEN!PROBLEMSPECIFIC
	   				P(NG)=PINF
	   				T(NG)=TINF
	   				V(NG)=-V(NP)
	   				RHO=P(NG)/(RR*T(NG))
	   				IF(YCC(NG).LE.DELTABL)THEN
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF)*		  &
	   						  (YCC(NG)/DELTABL)**(1.0D0/7.0D0) - U(NP)
	   				ELSE
	   					U(NG)= 2.0D0*AMINF*DSQRT(GAMA*RR*TINF) - U(NP)
	   				ENDIF
	   				amu = 1.0d0*(0.00000011848D0)*DSQRT(AMOL)*(T(NG)**0.6D0)
					RT= AMU/RHO
					C1(NG)=RHO
					C2(NG)=RHO*U(NG)
					C3(NG)=RHO*V(NG)
					E=CV*T(NG)+0.5D0*(U(NG)*U(NG)+V(NG)*V(NG))		
					C4(NG)=RHO*E
            		F1(NG)=C2(NG)
            		F2(NG)=RHO*U(NG)*U(NG)+P(NG)
            		F3(NG)=RHO*U(NG)*V(NG)
            		F4(NG)=(C4(NG)+P(NG))*U(NG)
            		G1(NG)=C3(NG)
            		G2(NG)=F3(NG)
            		G3(NG)=RHO*V(NG)*V(NG)+P(NG)
            		G4(NG)=(C4(NG)+P(NG))*V(NG)
            		C5(NG)=RHO*RT
					F5(NG)=C5(NG)*U(NG)
            		G5(NG)=C5(NG)*V(NG)
            	ENDIF
	   		ENDDO
	   		
        	CALL SOURCETERMS
        
        DO I=1,NELES
        	NCC1=NC1(I)
			NCC2=NC2(I)
			NCC3=NC3(I)
			NCC4=NC4(I)
			SX1=SC1X(I);SX2=SC2X(I);SX3=SC3X(I);SX4=SC4X(I)
			SY1=SC1Y(I);SY2=SC2Y(I);SY3=SC3Y(I);SY4=SC4Y(I)
			F1I=F1(I);F2I=F2(I);F3I=F3(I);F4I=F4(I);F5I=F5(I)
			G1I=G1(I);G2I=G2(I);G3I=G3(I);G4I=G4(I);G5I=G5(I)
			
    	   R1(I)=0.5D0*((F1I+F1(NCC1))*SX1+&
       					(G1I+G1(NCC1))*SY1+&
       		 			(F1I+F1(NCC2))*SX2+&
     					(G1I+G1(NCC2))*SY2+&
      					(F1I+F1(NCC3))*SX3+& 
     					(G1I+G1(NCC3))*SY3+&
          				(F1I+F1(NCC4))*SX4+& 
     					(G1I+G1(NCC4))*SY4)
    	   R2(I)=0.5D0*((F2I+F2(NCC1))*SX1+&
        				(G2I+G2(NCC1))*SY1+& 
     		 			(F2I+F2(NCC2))*SX2+& 
     					(G2I+G2(NCC2))*SY2+& 
           				(F2I+F2(NCC3))*SX3+& 
     					(G2I+G2(NCC3))*SY3+& 
           				(F2I+F2(NCC4))*SX4+& 
     					(G2I+G2(NCC4))*SY4)
           R3(I)=0.5D0*((F3I+F3(NCC1))*SX1+&
       					(G3I+G3(NCC1))*SY1+& 
         			 	(F3I+F3(NCC2))*SX2+& 
     					(G3I+G3(NCC2))*SY2+& 
         			 	(F3I+F3(NCC3))*SX3+& 
     					(G3I+G3(NCC3))*SY3+& 
         			 	(F3I+F3(NCC4))*SX4+&
     				 	(G3I+G3(NCC4))*SY4)
    	   R4(I)=0.5D0*((F4I+F4(NCC1))*SX1+&
        			 	(G4I+G4(NCC1))*SY1+&           	    
     	 				(F4I+F4(NCC2))*SX2+& 
     					(G4I+G4(NCC2))*SY2+& 
        	 			(F4I+F4(NCC3))*SX3+& 
     					(G4I+G4(NCC3))*SY3+& 
          				(F4I+F4(NCC4))*SX4+& 
     					(G4I+G4(NCC4))*SY4)
     	   R5(I)=0.5D0*((F5I+F5(NCC1))*SX1+&
           		 		(G5I+G5(NCC1))*SY1+&           	    
     	   				(F5I+F5(NCC2))*SX2+& 
     	   				(G5I+G5(NCC2))*SY2+& 
           				(F5I+F5(NCC3))*SX3+& 
     	      			(G5I+G5(NCC3))*SY3+& 
           				(F5I+F5(NCC4))*SX4+& 
    	   				(G5I+G5(NCC4))*SY4)-VOL(I)*SOURCE5(I)
		ENDDO		
	   	
	   	DO I=1,NELES
	   		DTBYVOL=DT(I)/VOL(I)
          C1(I)=COLD1(I)-ALPHA4*DTBYVOL*(R1(I)-D241(I))
          C2(I)=COLD2(I)-ALPHA4*DTBYVOL*(R2(I)-VISU(I)-D242(I))
          C3(I)=COLD3(I)-ALPHA4*DTBYVOL*(R3(I)-VISV(I)-D243(I))
          C4(I)=COLD4(I)-ALPHA4*DTBYVOL*(R4(I)-VIST(I)-D244(I))
  	  	  C5(I)=COLD5(I)-ALPHA4*DTBYVOL*CFLTURB*(R5(I)-VISR(I)-D245(I))
          IF(C5(I)/C1(I).le.RTOLERANCE)	C5(I)=C1(I)*RTOLERANCE
       	ENDDO
		
		CDIFF1=1.0E-10
        CDIFF2=1.0E-10
        CDIFF3=1.0E-10
        CDIFF4=1.0E-10
        CDIFF5=1.0E-20
   		IF(MOD(ITER,100)==0)THEN
          	DO I=1,NELES
             	CDIFF1=DMAX1(DABS(COLD1(I)-C1(I)),CDIFF1)
               	CDIFF2=DMAX1(DABS(COLD2(I)-C2(I)),CDIFF2)
               	CDIFF3=DMAX1(DABS(COLD3(I)-C3(I)),CDIFF3)
               	CDIFF4=DMAX1(DABS(COLD4(I)-C4(I)),CDIFF4)
               	CDIFF5=DMAX1(DABS(COLD5(I)-C5(I)),CDIFF5)
           	ENDDO
                
         	WRITE(32,11) CDIFF1,CDIFF2,CDIFF3,CDIFF4,CDIFF5
11			FORMAT(5E15.5)
            PRINT 10,ITER+NITEROLD,' ',CDIFF1,' ',CDIFF2,' '	&
          				,CDIFF3,' ',CDIFF4,' ',CDIFF5
10     		FORMAT(I10,A5,E12.6,A5,E12.6,A5,E12.6,A5,E12.6,A5,E12.5)
		ENDIF
		
       		DO I=1,NELES
            	COLD1(I)=C1(I)
                COLD2(I)=C2(I)
                COLD3(I)=C3(I)
                COLD4(I)=C4(I)
                COLD5(I)=C5(I)
			ENDDO
			
		IF(MOD(NITEROLD+ITER,20000)==0)THEN
			WRITE(CHOICE,"(A7,I10,A3)") "restart",NITEROLD+ITER,".in"
			PRINT*, CHOICE
			OPEN(31,FILE=CHOICE)
        	DO I=1,NELES
				WRITE(31,*) C1(I),C2(I),C3(I),C4(I),C5(I)
        	ENDDO   
        		write(31,*) RTRED
        		write(31,*) npseudo
        		wRITE(31,*) RINIT
        		WRITE(31,*) ITER+NITEROLD         	
      		CLOSE(31)
		endif	
	
		ENDDO		
		
		NOLD=1
		NITEROLD=NITEROLD+NITER
		PRINT*,' '
        PRINT*,'Continue further? (Y/N)'
        READ*, CHOICE
		PRINT*,' '
        IF(CHOICE=='y'.OR.CHOICE=='Y')THEN
        	write(27,*) "Simulation continued further:"
	  		PRINT*,'Modify, which of the following?'
     		PRINT*,'CFL (C)    Rtolerance switch (T)    Source term pseudo (S)',  &
     			'    None (N)' 
          	READ*,CHOICE
			IF(CHOICE=='c'.OR.CHOICE=='C')THEN
	  			PRINT*,'CFL new:'
	  			READ*,  CFL
	  			write(27,*) "CFL modified to: ",CFL
			ELSEIF(CHOICE=='t'.OR.CHOICE=='T')THEN
	  			PRINT*,'R Tolerance reduction factor:'
	  			print*,'(Current value: ',rtred,')'
				READ*,  RTRED
				write(27,*) "RTRED modified to: ",RTRED
				rtolerance=rinit*rtred
			ELSEIF(CHOICE=='s'.OR.CHOICE=='S')THEN
				PRINT*,'Turbulence source term:'
				print*,'Current pseudo value: ',npseudo
	  			READ*,  nPseudo
	  			write(27,*) "NPSEUDO modified to: ",NPSEUDO
	  		ELSEIF(CHOICE=='ct'.OR.CHOICE=='CT')THEN
	  			PRINT*,'CFL new:'
	  			READ*,  CFL
	  			PRINT*,'R Tolerance reduction factor:'
	  			print*,'(Current value: ',rtred,')'
				READ*,  RTRED
				write(27,*) "CFL modified to: ",CFL
				write(27,*) "RTRED modified to: ",RTRED
				rtolerance=rinit*rtred
			ELSEIF(CHOICE=='cs'.OR.CHOICE=='CS')THEN
				PRINT*,'CFL new:'
	  			READ*,  CFL
				PRINT*,'Turbulence source term:'
				print*,'Current pseudo value: ',npseudo
	  			READ*,  nPseudo
	  			write(27,*) "CFL modified to: ",CFL
	  			write(27,*) "NPSEUDO modified to: ",NPSEUDO
			ELSEIF(CHOICE=='ts'.OR.CHOICE=='TS')THEN
				PRINT*,'R Tolerance reduction factor:'
	  			print*,'(Current value: ',rtred,')'
				READ*,  RTRED
				PRINT*,'Turbulence source term:'
				print*,'Current pseudo value: ',npseudo
	  			READ*,  nPseudo
	  			write(27,*) "RTRED modified to: ",RTRED
	  			write(27,*) "NPSEUDO modified to: ",NPSEUDO
			ELSEIF(CHOICE=='cts'.OR.CHOICE=='CTS')THEN
				PRINT* ,'CFL new:'
	  			READ*, CFL
				PRINT*,'R Tolerance reduction factor:'
	  			print*,'(Current value: ',rtred,')'
				READ*,  RTRED
				rtolerance=rinit*rtred
				PRINT*,'Turbulence source term:'
				print*,'Current pseudo value: ',npseudo
	  			READ*,  nPseudo
	  			write(27,*) "CFL modified to: ",CFL
	  			write(27,*) "RTRED modified to: ",RTRED
	  			write(27,*) "NPSEUDO modified to: ",NPSEUDO
			ENDIF
		
			PRINT*,' '
			PRINT*,'Number of iterations:'
			READ*, NITER
			PRINT*,' '
			write(27,*) "Continued with another ",niter,' iterations' 
    		GOTO 100
		ENDIF
		write(27,*) "Simulation ended."
		write(27,*) " "
		
 		CLOSE(32)
 		
 		OPEN(31,FILE='restart.in')
        DO I=1,NELES
			WRITE(31,*) C1(I),C2(I),C3(I),C4(I),C5(I)
        ENDDO   
        	write(31,*) RTRED
        	write(31,*) npseudo
        	write(31,*) rinit 
        	write(31,*) niterold        	
      	CLOSE(31)
 		
		CALL POSTPROCESS!problem specific
		
		STOP
		
ENDPROGRAM MAINSOLVER


        
        
        !do i = 1, neles
		!	if (c5(i) .ne. c5(i)) then
		!	print*,'after jameson:'
		!		print*, 'iter: ',iter,'stage2', ',index: ',i
		!		print*,'Values:'
		!		print*,'ind:',i,' : ',cold5(i),c5(i)
		!		print*,'ind:',nc1(i),' : ',cold5(nc1(i)),c5(nc1(i))
		!		print*,'ind:',nc2(i),' : ',cold5(nc2(i)),c5(nc2(i))
		!		print*,'ind:',nc3(i),' : ',cold5(nc3(i)),c5(nc3(i))
		!		print*,'ind:',nc4(i),' : ',cold5(nc4(i)),c5(nc4(i))
				
			
		!	pause
		!	endif
		!enddo
		
		
		
	   		!print*,'after bdry'
			!do j=1,nghosts
			!np=nparent(j)
			!ng=nghost(j)
			!nt=ntype(j)
			!ns=nside(j)
			!	print*,np,c5(ng),c5(np),cold5(np),nt
			!enddo
			
			!pause

			!print*,'after viscous'
	   		
	   		!do i=1,neles
	   		!	print*,visu(i),visv(i),vist(i),visr(i)
	   		!enddo
        	!pause
	   		
	   		


        
