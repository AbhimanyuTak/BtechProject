/*
*--------------------------------------------------------------------*
*	MPI related declerations:

	Include 'mpif.h'
*/
        #include<math.h>
        int nproc,ierr;//status[mpi_status_size] no include file
	int tag=1,llim[101],ulim[101],rank=0;
/*	Include 'nxlim.h'
*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*
*   common definitions for XYfdtd.f                                  *
*--------------------------------------------------------------------*
        implicit none


*  define some fundamental constants
*  c=speed of light in vacuum, pi, eps0=permittivity of free space
*  xmu0=permeability of free space
*
*--------------------------------------------------------------------*/
        const int llimp=0,ulimp=420,nyy=2012;// array 110=ulimp-llimp+1 nyy=110
        int i,j,n;
	//const double c=2.99795,pi=3.14159265,eps0=8.854,xmu0=1.25663706,qe=1.602176487,cmasse=9.10938215,akb=1.3806503;
//	const double c=2.99795*pow(10.0,8.0),pi=3.14159265,eps0=8.854*pow(10.0,-12.0),xmu0=1.25663706*pow(10.0,-6.0),qe=1.602176487*pow(10.0,-19.0),cmasse=9.10938215*pow(10.0,-31.0),akb=1.3806503*pow(10.0,-23.0);
        const double c=2.99795e8,pi=3.14159265,eps0=8.854e-12,xmu0=1.25663706e-6,
                     qe=1.602176487e-19,cmasse=9.10938215e-31,akb=1.3806503e-23;
        double xd1,xd2,yd1,yd2,xxpr,timp,xxpr1;
        double ERMS[110][110],erms2[110][110],exrms[110][110],exrms2[110][110],eyrms[110][110],eyrms2[110][110];
        double den[110][110],dent[110][110];
        double exi[110][110],eyi[110][110],exi1[110][110],eyi1[110][110],exs[110][110],eys[110][110];
        double hzi[110][110],hzs[110][110];
        double vx[110][110],vy[110][110],ext[110][110],eyt[110][110];
        double exs1[110][110],exs2[110][110],eys1[110][110],eys2[110][110];
        double t,dt,ds,dte,dtm,dteds,dtmds,dtsi,DTAC;
        double c1,c2,c3,c4,xmid[500],ymid[500],sgdx0[500],sgdy0[500];
        double E0,FREQ,OMEG;
        double om2snu,QSM,crec,dnma,xfront;
        double FNUM,RECOMB,EMOB,EDIF,ETEM,nu2omeg2;
        double PARC,ACCEL,CORNUI,ndifmax;
        double DINI[100],TIMD,PRESSURE,TEMP0;
        double DENG0,fnucor,amu,nu_omeg,e2_epsm;
        double DIFFUSION[110][110],SIO[110][110],SIOT[110][110];
        double pow1[110][110],powt,frqio[110][110],puis[110][110],fnuimax,cene;
	double denp[110][110],difa,dife,frq[110][110];
        int IORDER,IABSOR,naccel;
        int nmaxwell,nx,ny,nperdt,nx1,ny1,nec;
        int isig,nlamb;
        int nstep,nini,imid[100],inul,jmid[100];



/*	parameter(c=2.99795d8,pi=3.14159265d0,
     >          eps0=8.854d-12,xmu0=1.25663706d-6,
     >          qe=1.602176487d-19,cmasse=9.10938215d-31,
     >          AKB=1.3806503d-23 )

        parameter(tag=1)
	dimension llim(0:100),ulim(0:100)
      	parameter(110=110)
*
*  array definitions for a simple wave propagation study
*


      	COMMON/DIFID/ DIFFUSION(llimp:ulimp,NYY),FRQIO(llimp:ulimp,NYY)

      	COMMON/FRMST/ ERMS(llimp:ulimp,110),ERMS2(llimp:ulimp,110),
     &	frq(llimp:ulimp,110)
      	COMMON/FRMSX/ EXRMS(llimp:ulimp,110),EXRMS2(llimp:ulimp,110)
      	COMMON/FRMSY/ EYRMS(llimp:ulimp,110),EYRMS2(llimp:ulimp,110)
      	COMMON/PLAD/den(llimp:ulimp,110),dent(llimp:ulimp,110),
     &	DENP(llimp:ulimp,NYY),DNMA
      	COMMON/EINC/ EXI(llimp:ulimp,NYY),EYI(llimp:ulimp,NYY)
      	COMMON/EIN1/ EXI1(llimp:ulimp,NYY),EYI1(llimp:ulimp,NYY)
      	COMMON/ESCA/ EXS(llimp:ulimp,NYY),EYS(llimp:ulimp,NYY)
      	COMMON/HINC/ HZI(llimp:ulimp,NYY)
      	COMMON/HSCA/ HZS(llimp:ulimp,NYY)
      	COMMON/VELO/ VX(llimp:ulimp,NYY),VY(llimp:ulimp,NYY)
      	common/ESCAT/EXT(llimp:ulimp,NYY),EYT(llimp:ulimp,NYY)
      	common/ESAVX/EXS1(llimp:ulimp,NYY),EXS2(llimp:ulimp,NYY)
      	common/ESAVY/EYS1(llimp:ulimp,NYY),EYS2(llimp:ulimp,NYY)
c      	COMMON/POWER/POW(llimp:ulimp,NYY),POWT,PUIS(llimp:ulimp,NYY)


      	common/STUF1/t,dt,ds,dte,dtm,dteds,dtmds,dtsi,DTAC
      	COMMON/STUF2/c1,c2,c3,c4,id1,id2,jd1,jd2,idom,n
      	COMMON/STUF3/xmid(500),ymid(500),sgdx0(500),sgdy0(500),DINI(100)
      	COMMON/STUF4/NINI,IMID(100),JMID(100)
      	Common/stuf5/xd1,xd2,yd1,yd2,xxpr,timp,xxpr1
      	COMMON/wavdf/E0,FREQ,OMEG,e2_epsm
      	COMMON/CONST/om2snu,QSM,fnucor,amu,nu_omeg
      	COMMON/FRQU/FNUM,RECOMB,CREC,CENE,EMOB,EDIF,ETEM,nu2omeg2
      	COMMON/DIFU/ DIFE,DIFA
      	COMMON/ACCE/PARC,ACCEL,CORNUI,ndifmax,naccel
      	COMMON/NMAXW/NMAXWELL
      	COMMON/METH/ISIG,IABSOR

      	COMMON/DENSI/TIMD
      	COMMON/DIM/nx,ny,nperdt,nlamb,nx1,ny1,nec,nstep
      	COMMON/GAS1/PRESSURE,TEMP0
      	COMMON/GAS2/ DENG0
      	COMMON/SION/SIO(llimp:ulimp,NYY),SIOT(llimp:ulimp,NYY)
      	COMMON/POSF/XFRONT,fnuimax
	common/commach/ rank,llim,ulim,nproc
        SAVE

*--------------------------------------------------------------------*/
