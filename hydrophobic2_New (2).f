C C       SOFT PARTICLE WITH HYDROPHILIC INNER CORE
C         OHSHIMA'S THEORY

c       ModifZEd in FEBRUARY 2018

	    implicit double precision(a-h,o-z)
 
	DIMENSION  
     $ AN(4,4),BN(4,4),CN(4,4),DN(4),JN(4,5000),B(4),
     $ A(4,4),DW(4,5000),EN(4,4,5000),FF(4,4),dxl(5000),dxr(5000),
     $ dyt(5000),
     $ dyb(5000),x(5000),y(5000),dx(5000),dy(5000),u(5000),
     $	u1(5000),
     $ v1(5000),v(5000),u2(5000),v2(5000),
     $  vp(5000), staws(5000) ,cdd(5000),dwb(2500,5000),
     $ addc1(5000),Ex(5000)
     $ ,steady(2500,5000),ana(5000),bmf(5000),bmf1(5000)
     @ ,phi(5000),phi1(5000),gmf(5000)
     @ ,gmf1(5000),hmf(5000),hmf1(5000),amf(5000),amf1(5000)
     @ ,fmf1(5000),fmf(5000),theta(5000),As(2500,5000)
     @ ,f(2500,5000),g(2500,5000),h1(2500,5000)
     @ ,S(2500,5000),T(2500,5000),CdmD(2500,5000)
     @ ,cdmE(2500,5000),cdms(2500,5000),P(2500,5000),
     @ hS(2500),HC(2500),dphi(5000),
     @ hS1(2500),HC1(2500),PP(2500),AA(2500),PHIP(2500),
     @ PHIS1(2500), PHIS2(2500), PHIS3(2500), PHIC1(2500), 
     @ PHIC2(2500),BFG(2500)
     @	,AFG(2500) ,TFG(2500)




ccccccccccccccccccccccccccccccccccccccc
c      ESSENTIAL PARAMETERS
ccccccccccccccccccccccccccccccccccccccc
	epse=78.54*8.854181787E-12
	visc=0.001
	ANo=6.023*1.0E23
	eC=1.602*1.0E-19
	Phi_0=0.0258
	CLam=7.35E-3
	ALam=7.146E-3

	Factor1=(aNo*eC/CLam)*(epse*Phi_0)/(visc)
	Factor2=(aNo*eC/ALam)*(epse*Phi_0)/(visc)

       Fara=96500
       b11=0.5*1.0E-5
     

cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     GRID GENERATION ALONG RADIAL (r) DIRECTION (X)
cCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        n1=1000
	n11=n1+1
	n12=n1+2
	 
		
 

 	 nAstar=31
 	 nAstar1=Nastar+1
 	 nBstar=41
 	 
	x(1)=0.0001
	DO I=1,n12
	DX(I)=0.025
	end do



	dx(1)=dx(2)
	dx(n12)=dx(n11)
	 
	do i=2,n12
	x(i)=x(i-1)+DX(i)
	end do
	
	OPEN(1,FILE='Dx.dat')
	DO i=1,n12
	WRITE(1,*)i,dx(i),x(i) 
 
	END DO
	CLOSE(1)
cccccccccccccccccccccccccccccccccccccccccccccc
c      save the needfull data files
cccccccccccccccccccccccccccccccccccccccccccccc

 	OPEN(222,FILE='OLD_kappa_b_MOB.dat')

	OPEN(1110,FILE='C_MOB.dat')
	OPEN(111,FILE='kappa_b_MOB.dat')
	OPEN(1111,FILE='T_MOB_d.dat')
cccccccccccccccccccccccccccccccccccc
	E_ratio=1.0
	

C     See the following (Line no 99 to  135)
        pH=6.49301 ! pH will be same as that of corrsponing PZM

        d=50.0       !kappa b for the corresponding PZM
        
        pka=4.58       !Change accordingly for Sio_2/Tio_2
        pkb=8.6        !Change accordingly for Sio_2/Tio_2
        aNt=0.5*1E-7   !Adjust  accordingly


        epse=80.0*8.854181787E-12
	phi0=0.0258

         

 
  	aNc=10.0
        bNc=10.0


         b1=15E-9

       
       

       SMa=(Fara*aNc*b1*b1)/(epse*phi0)
       SMb=(Fara*bNc*b1*b1)/(epse*phi0)
        
        pKa2=7.8!7 !7.8
	pkb2=-4.95!2 !-4.95



     
	Sig=Fara*aNt*b1/(epse*phi0)

	write(*,*)"b1,sig SMA,SMB=",b1,sig, SMA,SMB
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


11	continue



C	Concentration Calulation...
        Fara=96485
	epse=80.0*8.854181787E-12
	phi0=0.0258
	Fact1=sqrt((2.0*Fara)/(epse*Phi0))
c	write(*,*) Fact1
        aR=10E-9
        con1=(d/(Fact1*aR))*(d/(Fact1*aR))
ccccccccccccccccccccccccccccccccccccccccccccccccccccc


c	pause



	
                                          ! scaled PEL charge

	bt=1               !softness paramter beta
	ab=0.75            ! scaled inner core radius
	al=ab*bt
          

 

 

ccccccccccccccccccccNUMERICAL SOLUTION FOR POTENTIALcccccccccccccccccccc

 55      continue   

		
       itrn=0
       time=dt
       mn=1
       itrm=0
       mn1=0
       mstr=0
       m=1
        
       


 9        continue
5555    continue

555     continue

        itrn=itrn+1
        itrm=itrm+1
        mn1=mn1+1



c	----BOUNDARY CONDITIONS---------------------------
  

c	********** INNER CORE SURFACE *********************
 
	Aka2=10**(-pka2)
	Akb2=10**(-pkb2)
	Akc2=10**(-pH)

	Deno1=Akb2*(Akc2*Akc2*exp(-2*phi(Nastar)))
	Deno2=(Akc2*exp(-phi(Nastar)))+Aka2
	aNum=-Aka2+Akb2*(Akc2*Akc2*exp(-2*phi(Nastar)))

 
        fph=(aNum/(Deno1+Deno2))
     
c	write(*,*)"fFactor=",fpH


c	pause

        Sig1=Sig*fPh

      	  phi1(nAstar)=phi1(nAstar+1) 
     $  +Sig1*dx(nastar)

      	  phi(nAstar)=phi(nAstar+1) 
     $  +Sig1*dx(nastar)


c	********** r TENDS TO INFINITY ******************

 
          phi1(n12)=0.0
          phi(n12)=0.0

	 



 
c **************************************************
c       	   PHI CAL		                   *
c   	   	FOR WHOLE  REGION (r>GAMMA)                *
c ***************************************************
 
        DO 2915 i=nAstar1,n11

 
 


        AFG(i)=SMA*(1.0/(1.0+(10.0**(-pH+pKa))*exp(-phi(i))))
        BFG(i)=SMB*(1.0/(1.0+(10.0**(pH-pKB))*exp(phi(i))))
        TFG(i)= -(AFG(i)-BFG(i))


	i1=i-1

	dxr(i)=0.5*(dx(i)+dx(i+1))
	dxl(i)=0.5*(dx(i)+dx(i-1))
	 
 
cccccccccccccccnhenccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
	 if( (i.ge.nAstar1).and.(i.le.nBstar)) then

	AN(3,3)= 1.0/(dxl(i)*dxl(i))  
     $           -1.0/(dxl(i)*x(i))

 
 	BN(3,3)= -2.0/(dxl(i)*dxl(i)) -d*d
 
	CN(3,3)= 1.0/(dxl(i)*dxl(i))  
     $           +1.0/(dxl(i)*x(i))


        DN(3)=-TFG(i)!+d*d*Phi(i)
	else
cccccccccccccccnhenccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 	
   
	AN(3,3)= 1.0/(dxl(i)*dxl(i))  
     $           -1.0/(dxl(i)*x(i))

 
 	BN(3,3)= -2.0/(dxl(i)*dxl(i)) -d*d 
 
	CN(3,3)= 1.0/(dxl(i)*dxl(i))  
     $           +1.0/(dxl(i)*x(i))


        DN(3)=0.0!+d*d*Phi(i)
	end if
c 	write(*,*) AN(3,3),BN(3,3),CN(3,3),DN(3)
ccccccccccccccccccccccccccd ifcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 

 

        if(i.eq.nAstar1)dn(3)=dn(3)-an(3,3)*Phi1(nAstar1-1)
        if(i.eq.n11)dn(3)=dn(3)-cn(3,3)*Phi1(n12)
        if(i.ne.nAstar1)go to 2889
        cdd(nAstar1)=cn(3,3)/bn(3,3)
        addc1(nAstar1)=dn(3)/bn(3,3)
        go to 2915
2889      continue
        cdd(i)=cn(3,3)/(bn(3,3)-an(3,3)*cdd(i1))
        addc1(i)=(dn(3)-an(3,3)*addc1(i1))/(bn(3,3)-an(3,3)*cdd(i1))
2915      continue
        Phi1(n11)= addc1(n11)
        n45=n11-nAstar1
        do 2829 i=1,n45
        ii=n11-i
        i21=ii+1
        Phi1(ii)=addc1(ii)-cdd(ii)*Phi1(i21)


ccc	write(*,*) Phi1(i)

2829     continue
 





	bigphi=0.0
 
	 do 7654 i=Nastar1,n11
 
	  phiat=abs( (phi1(i)-phi(i))/phi(i))
  
	  if( phiat.gt.bigphi) bigphi=phiat
	  phi(i)=phi1(i)


7654	continue
 
	if(itrn.le.3001) then

  
c 	 write(*,*) 'bigphi itrn=',bigphi,itrn


	 goto 555
	endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c        CALCULATION COMPLETE FOR POTENTIAL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        do i=Nastar1,n11
        dphi(i)=(Phi(i+1)-Phi(i))/(1.0*dx(i))
        end do

	 
c        do i=2,n1
c	dphi(i)=(1.0 /(2.0*dx(i)))*
c     $  (-3.0*phi(i)+4.0*phi(i+1)-phi(i+2) )
c	end do

        dphi(n12)=(Phi(n12)-Phi(n11))/(dx(1))
c        dphi(n11)=dphi(n12) 
        dphi(Nastar)=-SIG!(Phi(2)-Phi(1))/(dx(1))


	open(133,file='dphi_dr.dat')
	do i=Nastar,n12
 	write(133,*)x(i), dphi(i)  
 	end do
	close(133)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       DATA FILES FOR NUMERICALLY COMPUTED POTENTIAL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


	If(d.eq.3.5) then
	open(13333,file='2phi.dat')
	do i=Nastar,n12
 	write(13333,*)x(i), phi(i)  
 	end do
	close(13333)
	end if


	If(d.eq.1) then
	open(1333,file='1phi.dat')	  
	do i=Nastar,n12
 	write(1333,*)x(i), phi(i)  
 	end do
	close(1333)
	end if

	If(d.eq.5) then
	open(1334,file='5phi.dat')	  
	do i=Nastar,n12
 	write(1334,*)x(i), phi(i)  
 	end do
	close(1334)
	end if
 

	If(d.eq.10) then
	open(1336,file='10phi.dat')	  
	do i=Nastar,n12
 	write(1336,*)x(i), phi(i)  
 	end do
	close(1336)
	end if

	If(d.eq.50) then
	open(1337,file='50phi.dat')	  
	do i=Nastar,n12
 	write(1337,*)x(i), phi(i)  
 	end do
	close(1337)
	end if

 	If(d.eq.100) then
	open(13375,file='100phi.dat')	  
	do i=Nastar,n12
 	write(13375,*)x(i), phi(i)  
 	end do
	close(13375)
	end if
	
	


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Mobility calculation based on Ohshima's semi-analytical formula
C                        Ohshima, 1994,1995 paper
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       List of parameters
        bts1=0.0    !(bts1:Nastar,n12)
        
1234   continue
        
        bts2=(1.0/ab)*bts1
        alp1=1.001!(alp1:Nbstar,n100)
        alp2=ab*alp1
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                     Mi's values
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	aM11=((3.0*ab)/(2.0*alp1*alp1))-((3.0*ab*ab)/(2.0*alp1*alp1))+
     $	((ab*ab*ab)/(2.0+(4.0*bts2)))+(1.0/(1.0+(2.0*bts2)))-
     $  ((ab*ab*ab*ab*bts2)/(2.0+(4.0*bts2)))-((ab*bts2)/(1.0+(2.0*bts2)
     $ )) +((3.0*ab*ab*ab*bts2)/(2.0+(4.0*bts2))) + ((3.0*bts2)/
     $  (1.0+(2.0*bts2)))

	aM11=aM11*cosh(alp1-alp2)

	aM12=((3.0*ab)/(2.0*alp1*alp1*alp1))-((3.0*ab*ab)/(2.0*alp1))+
     $ ((ab*ab*ab)/((2.0*alp1)+(4.0*bts2*alp1)))+
     $ (1.0/(alp1+(2.0*bts2*alp1)))-((ab*ab*bts2*alp2*alp2)/((2.0*alp1)
     $ +(4.0*bts2*alp1)))-((alp2*bts2)/(1.0+(2.0*bts2)))+
     $ ((3.0*ab*ab*ab*bts2)/((2.0*alp1)+(4.0*bts2*alp1)))+((3.0*bts2)/
     $ (alp1+(2.0*bts2*alp1)))


	aM12=aM12*sinh(alp1-alp2)

	aM1=aM11-aM12

	aM21=((3.0*ab*ab)/((2.0*alp1*alp2)+(4.0*alp1*alp2*bts2)))+
     $ ((ab*ab*ab)/(2.0+(4.0*bts2)))+(1.0/(1.0+(2.0*bts2)))+
     $ ((3.0*bts2*ab)/((alp1*alp1)+(2.0*bts2*alp1*alp1))) +
     $ ((3.0*ab*ab*ab*bts2)/(2.0+(4.0*bts2)))+
     $  ((3.0*bts2)/(1.0+(2.0*bts2)))

	aM21=aM21*cosh(alp1-alp2)

	aM22=((3.0*ab*ab)/((2.0*alp1)+(4.0*bts2*alp1)))+
     $ ((3.0*ab*ab*bts2)/(alp1+(2.0*bts2*alp1)))+
     $ ((alp2*bts2)/(1.0+(2.0*bts2)))+
     $ ( (alp2*alp2*ab*ab*bts2)/((2.0*alp1)+(4.0*alp1*bts2)) )


	aM22=aM22*sinh(alp1-alp2)

       aM2=aM21+aM22-((3.0*ab*ab)/((2.0*alp1*alp2)+(4.0*alp1*alp2*bts2))
     $ )-((9.0*ab*bts2)/((2.0*alp1*alp1)+(4.0*alp1*alp1*bts2)))
     $ +((3.0*ab*bts2)/((2.0*alp1*alp1)+(4.0*alp1*alp1*bts2)))


	aM31=(1.0/(1.0+(2.0*bts2)))+((3.0*bts2)/(1.0+(2.0*bts2)))-
     $ ((ab*bts2)/(1.0+(2.0*bts2)))

         aM31=aM31*cosh(alp1-alp2)

         aM32= (1.0/(alp1+(2.0*bts2*alp1)))+((3.0*bts2)/
     $ (alp1+(2.0*alp1*bts2)))-((alp2*bts2)/(1.0+(2.0*bts2)))

       aM32=aM32*sinh(alp1-alp2)


	aM3=aM31-aM32-ab


         aM41=(1.0/(1.0+(2.0*bts2)))-(bts2/(1.0+(2.0*bts2)))+((3.0*bts2)
     $   /((alp1*alp2)+(2.0*alp1*alp2*bts2)))

         aM41=aM41*sinh(alp1-alp2)

         aM42=(1.0/(alp1+(2.0*alp1*bts2)))+((3.0*bts2)/
     $ (alp2+(2.0*alp2*bts2)))-(bts2/(alp1+(2.0*alp1*bts2)))

         aM42=aM42*cosh(alp1-alp2)

         aM4=aM41-aM42+(1.0/alp1)+((2.0*alp1)/((3.0*ab)+(6.0*ab*bts2)))
     $   +((alp2*ab)/(3.0+(6.0*bts2)))


         aM51=((3.0*ab*ab*bts2)/(2.0+(4.0*bts2)))-((3.0*ab*bts2)/
     $ (2.0+(4.0*bts2)))

         aM51=aM51*cosh(alp1-alp2)

         aM52=((3.0*ab*bts2)/((2.0*alp1)+(4.0*bts2*alp1)))-((3.0*ab*
     $ bts2*alp2)/(2.0+(4.0*bts2)))

         aM52=aM52*sinh(alp1-alp2)

         aM5=aM51+aM52-((alp1*alp2*bts2)/(1.0+(2.0*bts2)))-((ab*ab*alp2
     $  *alp2*bts2)/(2.0+(4.0*bts2)))


	

         
         

         
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Gr(i)=PP(i)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	  DO i=1,n12
	  PP(i)=-(1.0+ab*ab*ab*(0.5/(x(i)*x(i)*x(i))) )*dPhi(i)
 	  end do



	aMu1=0.0
	DO i=nBstar,n12
	aMu1=aMu1+(d*d/9.0)*(  3.0*(1.0-x(i)*x(i))
     $	-2*(aM2/aM1)*(1.0-x(i)*x(i)*x(i))    )*PP(i)*dx(i)

     	end do
 


        
	aMu2=0.0
	  DO i=nAstar,n12
	aMu2=aMu2+(2.0/3.0)*(aM3/aM1)*(d*d/(alp1*alp1))*
     $	( (1.0+0.5*x(i)*x(i)*x(i))    )*PP(i)*dx(i)
	end do
 
	  DO i=1,n12
	AA(i)=alp1*(x(i)-ab)
       	end do
	
	
	aMu3=0.0
	  DO i=nastar,nBstar

        

	aMu3=aMu3-(2.0*d*d/(3.0*alp1*alp1))
     $	*(  1-1.5*(ab/(aM1*alp1*alp1))*
     $   ( (aM3+(aM4*alp1*x(i)))*cosh(AA(i))-
     $     (aM4+(aM3*alp1*x(i)))*sinh(AA(i)) )   )
     $	*PP(i)*dx(i)
        end do
     
        aMu4=0.0
            DO i=nastar,nBstar
            
            
	  
	  aMu4=aMu4-  ( ((2.0*d*d)/(3.0*alp1*alp1))*(aM5/aM1)*(  (x(i)/
     $	 alp1)-(3.0/(alp1*alp1*alp2)) )
     $   *sinh(AA(i))*PP(i)
     $   *dx(i)  ) - ( (2.0*d*d/(3.0*alp1*alp1))*(aM5/aM1)*(((3.0*x(i))/
     $	 ( alp1*alp2))-(1.0/(alp1*alp1)) )*cosh(AA(i)
     $ )*PP(i)  *dx(i)  )
       	end do
     
       
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Mobility calculation based on Ohshima's semi-analytical formula
c       for hydrophilic particle
C

       bt=alp1
       al=alp2

c                Ohshima, 1994,1995 paper
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                     Li's values
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	aL11=1.0+(0.5*ab*ab*ab)+(1.5/(bt*bt))*ab-((1.5/(bt*bt))*ab*ab)
	aL11=aL11*cosh(bt-al)

	aL12=1.0-(1.5*ab*ab)+(0.5*ab*ab*ab)+((1.5/(bt*bt))*ab)
	aL12=aL12*(1.0/bt)*sinh(bt-al)

	aL1=aL11-aL12

	aL2=1.0+(0.5*ab*ab*ab)+((1.5/(bt*bt))*ab)
	aL2=aL2*cosh(bt-al)
	aL2=aL2+1.5*ab*ab*(1.0/bt)*sinh(bt-al)-(1.5/(bt*bt))*ab


	aL3=cosh(bt-al)-((1.0/bt)*sinh(bt-al))-ab

	aL4=sinh(bt-al)-((1.0/bt)*cosh(bt-al))
	aL4=aL4+((al/3.0)*ab)+((2.0/3.0)*(bt/ab))+1.0/bt
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      list of parameters
         bt=alp1     !(alp1:Nbstar,n100)
         al=ab*bt
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       Gr(i)=PP(i)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



	OMu1=0.0
	DO i=nBstar,n12
	OMu1=OMu1+(d*d/9.0)*(  3.0*(1-x(i)*x(i))
     $	-2*(aL2/aL1)*(1.0-x(i)*x(i)*x(i))    )*PP(i)*dx(i)
	end do
	
	

	OMu2=0.0
	  DO i=nAstar,n12
	OMu2=OMu2+(2.0/3.0)*(aL3/aL1)*(d*d/(bt*bt))*
     $	( (1.0+0.5*x(i)*x(i)*x(i))    )*PP(i)*dx(i)
	end do

	  DO i=1,n12
	AA(i)=bt*(x(i)-ab)
	end do
	OMu3=0.0
	  DO i=nastar,nBstar
	OMu3=OMu3-(2.0*d*d/(3.0*bt*bt))
     $	*(  1-1.5*(ab/(aL1*bt*bt))*
     $   ( (aL3+aL4*bt*x(i))*cosh(AA(i))-
     $     (aL4+aL3*bt*x(i))*sinh(AA(i)) )   )
     $	*PP(i)*dx(i)
	end do
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	write(1110,*)con1, aMu1+aMu2+aMu3+aMu4

 	write(111,*)pH, aMu1+aMu2+aMu3+aMu4

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     See the following (Line no 666 to  706)
C      Important part_involves the calculation of effective and total charges
 	
 	
C     calculation of Qeff (CQ) and Q_fix (AQ) star
        DO  i=nAstar1,n11
        AFG(i)=SMA*(1.0/(1.0+(10.0**(-pH+pKa))*exp(-phi(i))))
        BFG(i)=SMB*(1.0/(1.0+(10.0**(pH-pKB))*exp(phi(i))))
        TFG(i)= -(AFG(i)-BFG(i))
       ENd do



        Q_fixS=0.0
        DO i=nastar1,nbstar
        Q_fixS=Q_fixS+(3.0/(1.0-ab*ab*ab))*TFG(i)*x(i)*x(i)*dx(i)
        end do
        AQ=Q_fixS  !AQ is the Q_fix star, charge due to PEL immobile charges
        
        
        Q_eff1=0.0
       DO i=nastar1,nbstar
       Q_eff1=Q_eff1-(3.0/(1.0-ab*ab*ab))*(0.5*(phi(i)+phi(i+1)))
     $  *d*d*x(i)*x(i)*dx(i)
        end do
        BQ=Q_eff1      !BQ is the neutralisation of charge due
                       !to penetration of mobile ions within PEL
                       
        CQ=AQ+BQ       ! is the effetive charge within PEL, i.e.,Q_eff


C
C        OPEN(1222,FILE='pH_Qfix_Q_eff_STAR.dat')
C        write(1222,*)pH, Q_fixS,Q_eff1,CQ, Q_eff1+Q_fixS
        
        
        OPEN(1223,FILE='akb_Qfix_STAR.dat')
        write(1223,*)d, AQ
        OPEN(1224,FILE='bkb_Q_eff_STAR.dat')
        write(1224,*)d, CQ
 	write(*,*)  d,AQ,BQ,CQ

C         if(pH.le.12) then
C	 pH= pH+1
C	itrn=1
C	go to 11
C	end if

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	stop
	end


	
	
	

     







