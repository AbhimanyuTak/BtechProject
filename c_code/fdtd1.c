/*!=============================================================
      program FDTDPLASMA
*
*  Two dimensional parallelised (MPI-based) Yee algorithm for the TM case
*  Absorbing (Mur 2nd order) boundary conditions are applied.
*  Plasma described by ambilor diffusion + ioniZation
!=============================================================
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "xyfdtd.h"
#include <time.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>

int KELEC,KRMS,nmod,NECRIR,icpling,K;
double output[421][2012];	//ulimp -> llimp:ulimp ??
double slambx,slamby,tstop, radius = 0;
int length,istart,iend,icc;
time_t itime;	//long
      
FILE *fptr,*fptr1,*fptr2,*fptr3, *file_a, *file_xelec, *file_yelec, *file_eyi;

void SETUP();
void INC_EFIELD();
void ADV_EFIELD_vel();
void ADV_HFIELD();
void RMS(int k);
void MR_MUR();
void ECRIR();
void zero();

int main()
{
//!=============================================================
//char *date= (char *) malloc(sizeof(char) * 30);
//fptr=fopen("output1.txt","w");
//time(&itime);
//date=ctime(&itime);
//fprintf(fptr,"%s %c","Program was started at",date);
//fclose(fptr);
struct timeval begin, end, total_start, total_end;
double t1, t2, t_inc_efield, t_adv_efield_vel, t_mr_mur, t_hfield;
gettimeofday(&begin, NULL);
t1 = begin.tv_usec;
gettimeofday(&total_start, NULL);
      nini=1;
      fptr=fopen("xstart.dat","r");

      char *line = NULL;
      size_t len = 0;
      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      tstop=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      nlamb=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      slambx=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      slamby=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      E0=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      PRESSURE=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      FREQ=atof(line);
      len=0;


      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      IABSOR=atof(line);
      len=0;

//      getline(&line, &len, fptr);
//      getline(&line, &len, fptr);
          
      if(nini<=0)
	  nini=1;
      K=1;
      do{
      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      xmid[K]=atof(line);
      len=0;


      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      ymid[K]=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      sgdx0[K]=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      sgdy0[K]=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      DINI[K]=atof(line);
      len=0; 
      K++;
      }while(K<=nini);
      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      icpling=atof(line);
      len=0;


      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      nmaxwell=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      ndifmax=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      naccel=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      NECRIR=atof(line);
      len=0; 
    
      icpling=0;
      
      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      crec=atof(line);
      len=0;

      getline(&line, &len, fptr);
      getline(&line, &len, fptr);
      len=10;
      getline(&line, &len, fptr);
      cene=atof(line);
      len=0; 

      fclose(fptr);

      gettimeofday(&end, NULL);
	  printf("Reading Input time: %f s\n", ((end.tv_sec - begin.tv_sec) + 
	  	((end.tv_usec - begin.tv_usec)/1000000.0)));

//!=============================================================
       gettimeofday(&begin, NULL);
	   t1 = begin.tv_usec;

	   nx=nlamb*slambx;
       ny=nlamb*slamby;
       nperdt=2.0*nlamb;
       nx1=nx-1;
       ny1=ny-1;
       t=0;
       nec=0;
       nstep=0;
//!=============================================================

//!=============================================================
      
	  SETUP();
	  gettimeofday(&end, NULL);
	  //t2 = end.tv_usec;
	  printf("Initialization time: %f s\n", ((end.tv_sec - begin.tv_sec) + 
              ((end.tv_usec - begin.tv_usec)/1000000.0)));
      //ELEC_DENS();
//*----



      printf("computing up to time  %f \n",tstop);
        
//!=============================================================
      fptr=fopen("output.txt","w");


       if (rank==0)
	{
//C--------------- Initial parameters ouput-----
		fprintf(fptr,"nx=%d   ny=%d",nx,ny);
		fprintf(fptr,"DS=%f   DT=%f",ds,dt);
		fprintf(fptr,"Freq=%f   Omega=%f",FREQ,OMEG);
		fprintf(fptr,"Time period=%f",1.0/FREQ);
		fprintf(fptr,"Lambda=%f",3.0e8/FREQ);
  	   
//!put 0	
		fprintf(fptr,"Collision Freq=%f",FNUM);
		fprintf(fptr,"Recombination Coef=%f",RECOMB);
		fprintf(fptr,"Mobility=%f",EMOB);
		fprintf(fptr,"Electron Temp= (eV)%f",ETEM);
		fprintf(fptr,"DIFFUSION Coef=%f",EDIF);
		fprintf(fptr,"Initial Gas/Neutral density=%f",DENG0);
		fprintf(fptr,"Cutoff-density=%f",(eps0*cmasse/pow(qe,2))*(pow(OMEG,2)+pow(FNUM,2)));// check formula
                fprintf(fptr, " Tstop %f",tstop); 
        	fclose(fptr);
	}
//!=============================================================

      KELEC=0;
      n=0;
int f=0;

fptr1=fopen("eyi1.dat","w");
fptr2=fopen("eyt.dat","w");
fptr3=fopen("hzs.dat","w");
file_eyi=fopen("file_eyi.dat","w");
//!=============================================================
      do{

         n=n+1;
         printf("%d \n", n);
         nstep += 1;
         //TIMD = 10;
         if(TIMD>tstop){
         	printf("%e %d %d\n", dt,nx,ny); 
         	printf("OUTPUT 1\n");
      		gettimeofday(&total_end, NULL);
      		//total_end = end.tv_usec;
      		printf("INC_EFIELD time: %f s\n", t_inc_efield);
      		printf("ADV_EFIELD_vel time: %f s\n", t_adv_efield_vel);
     	    printf("MR_MUR time: %f s\n", t_mr_mur);
      		printf("ADV_HFIELD time: %f s\n", t_hfield);
	  		printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + 
	  			((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
		 	exit(0);
		 }
         
         gettimeofday(&begin, NULL);
		 INC_EFIELD();
        fprintf(file_eyi,"%e \n", eyi1[nx/2][ny/2]);
         gettimeofday(&end, NULL);
	  	 t_inc_efield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

		 fprintf(fptr1,"%d %.15f %.15f\n",n, eyi[nx/2][ny/2], eyi1[nx/2][ny/2]);

		 gettimeofday(&begin, NULL);
		 ADV_EFIELD_vel();
		 gettimeofday(&end, NULL);
         t_adv_efield_vel += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

		 fprintf(fptr2,"%d %e %e \n",n, eyi1[nx/2][ny/2], eys[nx/2][ny/2]);
         //fprintf(fptr2,"%d %f\n",n, eyt[nx/2][ny/2]);
         gettimeofday(&begin, NULL);
	  	 MR_MUR();
	  	 gettimeofday(&end, NULL);
		 t_mr_mur += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

         t = t + dt/2.0;

         gettimeofday(&begin, NULL);
		 ADV_HFIELD();
		 gettimeofday(&end, NULL);
		 t_hfield += ((end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0));

		 fprintf(fptr3,"%d %.15f \n",n,hzs[nx/2][ny/2]);
         t = t + dt/2.0;

         KRMS=1;
         if(n==1000){
		   	//ECRIR();
		   	file_a=fopen("file_a.dat","w");
		   	//fprintf(fptr,"%s",tiret);	
        	//fprintf(fptr,"%s","  PLASMA DENSITY ");
        	//fprintf(file_a,"%d %d \n",nx,ny);
        	for(j=1;j<=ny;j++){
	   			for(i=1;i<=nx;i++){
   	       			fprintf(file_a,"%e ",den[i][j]);
   	       		}
   	       		fprintf(file_a,"\n");
	   		}
   	       	fclose(file_a);
   	       	printf("Hello\n");
		   	//exit(0);
		  }

		  double e_total;                                                                                                                                 
		  if(n==1000){
		   	//ECRIR();
		   	file_xelec=fopen("file_xelec.dat","w");
		   	//file_yelec=fopen("file_yelec.dat","w");
		   	//fprintf(fptr,"%s",tiret);	
        	//fprintf(fptr,"%s","  PLASMA DENSITY ");
        	//fprintf(file_a,"%d %d \n",nx,ny);
        	for(j=1;j<=ny;j++){
	   			for(i=1;i<=nx;i++){
	   				e_total = sqrt(pow(eyt[i][j],2) + pow(ext[i][j],2));
   	       			fprintf(file_xelec,"%e ",e_total);
   	       			//fprintf(file_yelec,"%e ",eyt[i][j]);
   	       		}
   	       		fprintf(file_xelec,"\n");
   	       		//fprintf(file_yelec,"\n");
	   		}
   	       	fclose(file_xelec);
   	       	//fclose(file_yelec);
   	       	printf("Hello\n");
		   	exit(0);
		  }
//!=============================================================
         if(fmod(n,nperdt)==0) 
         {
				KRMS=2;
                //RMS(KRMS);
	 }
//!=============================================================

         if(fmod(n,nperdt*nmaxwell)==0) 
	     {

		  /*if(icpling!=0)
		         ELEC_DENS();
		   else
		   {
		    */

		     DTAC=1.0/FREQ;
             TIMD=(float)(n)/(float)(nperdt)*nmaxwell*DTAC;

		   //}
		 
		   //if(nx/2>=llim[rank] && nx/2<=ulim[rank])
	             //front();// ! check
		 
//!=============================================================XXX
		   KELEC=KELEC+1;
		   
		   if(fmod(KELEC,NECRIR)==0)
		   {
		     nec=nec+1;
		     ECRIR();
		   }

//!=============================================================XXXXX
		   if(rank==0)
		   {
			   //printf("%d %f %f %f %f",n,TIMD,DTAC,ACCEL,dnma);
	
			   if(KELEC==1)
			      fptr=fopen("RES.out","w");
			   if(KELEC>1)
			   { 
				fptr=fopen("RES.out","a");
				if(fptr!=NULL)
				   fprintf(fptr,"%d %f %f %f %f",n,TIMD,DTAC,ACCEL,dnma);
				else
				   goto l100;
				 fclose(fptr);
			   }
                    }
	}  
l100:
//!=============================================================XXXXX

         nmod=2*nperdt;


      }while(1);

	  fclose(fptr1);
	  fclose(fptr2);
	  fclose(fptr3);
	  fclose(file_eyi);
      //fptr=fopen("OUTPUT2.txt","w");
      //time(&itime);
      //date=ctime(&itime);
      //fprintf(fptr,"%s %c","Program was ends at",date);
      //fclose(fptr);

      //exit(0);
      printf("OUTPUT 2\n");
      gettimeofday(&total_end, NULL);
      		//total_end = end.tv_usec;
      		printf("INC_EFIELD time: %f s\n", t_inc_efield);
      		printf("ADV_EFIELD_vel time: %f s\n", t_adv_efield_vel);
     	    printf("MR_MUR time: %f s\n", t_mr_mur);
      		printf("ADV_HFIELD time: %f s\n", t_hfield);
	  		printf("Total time: %f s\n", ((total_end.tv_sec - total_start.tv_sec) + 
	  			((total_end.tv_usec - total_start.tv_usec)/1000000.0)));
      
      return 0;
}
/* 
!=============================================================

!=============================================================
c     this SUBROUTINE initialiZes the computations
*/
      void SETUP()
      {

        static int k;
        static double ardiy,ardix,dinig;
        static double xd0,yd0,xxi,yyj;
        static double aaa[2012],bbb[2012],ccc[2012],ddd[2012];


		istart = 1;
		iend = nx;
		/// check with sir///
       dt=1.0/(float)(nperdt)/FREQ;
       ds=c/(float)(nlamb)/FREQ;
//!=============================================================
       OMEG=2.0*pi*FREQ;
       QSM=qe/cmasse;
       e2_epsm=qe*qe/eps0/cmasse;
/*
!!
!  AIR DAT
!  FNUM=electron-neutral coll frequency
!  RECOMB=electron-ion recombination coefficient
!  EMOB=electron mobility
!  ETEM=electron temperature
!  DIFE=electron diffusion coeff
*/
      FNUM=5.3*pow(10,9)*PRESSURE;   //d9 ??
      RECOMB=crec*1.0*pow(10,-13);    // d-13 ??
      EMOB=QSM/FNUM;
      ETEM=2.0*abs(cene);
      dife=EMOB*ETEM;
      difa=dife/100.0;
      nu2omeg2=pow(OMEG,2)+pow(FNUM,2);
      nu_omeg=FNUM/OMEG;
//!!

//!================Initial density location =======================

        for(K=1;K<=nini;K++)
        {
          if(xmid[K]>0) 
           imid[K]=xmid[K]*nx;
          if(ymid[K]>0) 
           jmid[K]=ymid[K]*ny;
          if(xmid[K]==0) 
           imid[K]=nx*0.5;
          if(ymid[K]==0) 
           jmid[K]=ny/2;
        }

//!=============================================================

      TEMP0=300.0;
      DENG0=PRESSURE/760.0*101300.0/akb/TEMP0;


      	  radius = nx/5*ds;
      	  printf("radius: %f\n", radius);
      	  printf("ds: %f\n", ds);
         for(j=1;j<=ny;j++)
         {
           for(i=1;i<=nx;i++)      
           {
             denp[i][j]=0.0;
             den[i][j]=0.0;
             if(pow((i-(nx/2))*ds,2)+pow((j-(ny/2))*ds,2)<=pow(radius,2))
	              den[i][j]=1e20;
	          	  //printf("hello\n");
           }
         }

         for(j=1;j<=ny;j++)
         {
          for(i=1;i<=nx;i++)      
           {
            ERMS[i][j]=E0/sqrt(2.0);
           }
          }      

      for(K=1;K<=nini;K++)
      {
        xd0=ds*imid[K];
        yd0=ds*jmid[K];
 
        if(xmid[K]<0) 
         xd0=-xmid[K];
        if(ymid[K]<0) 
         yd0=-ymid[K];
      
//!================Initial density, Gaussian, defined =======================
//!make DEN and DENP =0
        for(i=1;i<=nx;i++)
        {
          xxi=ds*i;
          ardix=0.0;
          if(sgdx0[K]>0)
            ardix=(-pow((xxi-xd0),2))/2.0/sgdx0[K]/sgdx0[K];
          for(j=1;j<=ny;j++)
          {    
		yyj=ds*j;
	        ardiy=0.0;
		 if(sgdy0[K]>0) 
		   ardiy=-pow((yyj-yd0),2)/2.0/sgdy0[K]/sgdy0[K];
                dinig=DINI[K]*exp(ardix+ardiy);
                if(dinig<=DINI[K]*1.0*exp(-2)) 
    		   dinig=0;
                //den[i][j]=den[i][j]+dinig;
                //denp[i][j]=den[i][j];

          }
         }
      }

//c     find relative delay for x, y, Z cell displacement

      dte=dt/eps0;
      dtm=dt/xmu0;
      dteds=dte/ds;
      dtmds=dtm/ds;
      //printf("dtm: %f \n", dtm);
      //printf("xmu0: %f \n", xmu0);
      //printf("dtmds: %f \n", dtmds);
/*c
c   mur constants
c
*/
      c1=(c*dt-ds)/(c*dt+ds);
      c2=2.0*ds/(c*dt+ds);
      c3=(c*dt*c*dt)/(2.0*ds*(c*dt+ds));


      }
//!=============================================================
      void INC_EFIELD()
      {
/*c3\
	fptr4=fopen("file_eyi.dat","w");
c   Subroutine to calculate INCIDENT field plane wave
c*/
   //   INCLUDE 'xyfdtd.inc'
        float sine,sine1,x,x0;
        float sineb,sine1b,x0b;
		int i0;
        i0=2;
        x0=(i0-1)*ds;
        x0b=(nx-1)*ds;
        for(i=1;i<=nx;i++)
        {
          x=(i-1)*ds;
          sine=0.0;
          sine1=0.0;
          sineb=0.0;
          sine1b=0.0;
          if(x<=(x0+c*t))  sine = sin(OMEG*(t-(x-x0)/c));
          if(x<=(x0+c*(t+dt)))  sine1 = sin(OMEG*(t+dt-(x-x0)/c));
          if(x>=x0b-c*t)  sineb = sin(OMEG*(t+(x-x0b)/c));
          if(x>=(x0b-c*(t+dt)))  sine1b = sin(OMEG*(t+dt+(x-x0b)/c));
          for(j=1;j<=ny;j++) 
          {
//	      	eyi[i][j] =   E0*(sine+sineb);
            //eyi1[i][j] =  E0*(sine1+sine1b);
            eyi1[i][j] =  E0*(sine1b);
            //file_xelec=fopen("file_xelec.dat","w");

           //eyi[i][j] =   E0*(sine);
           //eyi1[i][j] =  E0*(sine1);
          }
        }

      }
//!=============================================================
      void ADV_EFIELD_vel()
      {
      
        double aa,beta,alpha,gamma1,extk,eytk,const7,const8;
        double omp2,qmdt,const3,const4,const5,const6;
        double aaa[2012],bbb[2012],ccc[2012],ddd[2012];
        double fff[2012],ggg[2012],hhh[2012],jjj[2012];
        int ierror,icountm,icountm1,icounts,icounts1;

//!!keep the variables as it is below

          //printf("Beta: %f \n", beta);      	
          qmdt=qe/cmasse*dt;  
          aa=FNUM*dt/2.0; 
          alpha=(1.0-aa)/(1.0+aa);
          gamma1=1+aa;
          const3=.50*qe*qe/eps0/cmasse;
          const4=dt*dt/4.0/gamma1;
          const5=1.0/(1.0+beta);
          const6=1.0-beta;
          const7=.25*dte*(1.0+alpha);
          const8=1.0/2.0/gamma1;
          //printf("const5: %f \n", const5);

//!---
      for(j=2;j<=ny-1;j++)
      {
        for(i=1;i<=nx-1;i++)
         {
          omp2=(den[i][j]+den[i+1][j])*const3;
          beta=omp2*const4;

          extk=ext[i][j];

            exs[i][j]=const5*( exs[i][j]*(const6)+qe*(den[i][j]+den[i+1][j])*vx[i][j]*const7
            	      -(exi[i][j]+exi1[i][j])*beta+(hzs[i][j]-hzs[i][j-1])*dteds); 

            ext[i][j]=exs[i][j]+exi1[i][j];
            vx[i][j]=vx[i][j]*alpha - qmdt*(ext[i][j]+extk)*const8;
 
         }
      }

        for(j=1;j<=ny-1;j++)
        {
         for(i=2;i<=nx-1;i++)
          {
           omp2=(den[i][j]+den[i][j+1])*const3;
           beta=omp2*const4;

           eytk=eyt[i][j];

            eys[i][j]=const5*(eys[i][j]*(const6)+qe*(den[i][j]+den[i][j+1])*vy[i][j]*const7
            	      -(eyi[i][j]+eyi1[i][j])*beta-(hzs[i][j]-hzs[i-1][j])*dteds);

            eyt[i][j]=eys[i][j]+eyi1[i][j];
            vy[i][j]=vy[i][j]*alpha - qmdt*(eyt[i][j]+eytk)*const8;

          }
        }


//!! just one statement from the three loops (check I)

        if (rank==0)
        {
         for( j=1;j<=ny;j++)
          {
           for (i=1;i<=nx;i++)   
            {
             eyt[i][j]=eys[i][j]+eyi1[i][j];
             ext[i][j]=exs[i][j]+exi1[i][j];
            }
          }
        }

        if (rank>0 && rank< (nproc-1))
        {
         for(j=1;j<=ny;j++)
          {
           for(i=0;i<=nx+1;i++)
            {
             eyt[i][j]=eys[i][j]+eyi1[i][j];
             ext[i][j]=exs[i][j]+exi1[i][j];
            }
        
          } 
        }

        if (rank==(nproc-1))
        {
         for(j=1;j<=ny;j++)
         {
          for(i=0;i<=nx;i++)
           {
            eyt[i][j]=eys[i][j]+eyi1[i][j];
            ext[i][j]=exs[i][j]+exi1[i][j];
           }
         }      
        }

      }
//!=============================================================
      void ADV_HFIELD()
      {
/*
c
c  Subroutine to step 
c  currently for vacuum
c
      INCLUDE 'xyfdtd.inc'
*/

        for(j=1;j<=ny;j++)
          for(i=1;i<=nx;i++)
            hzs[i][j]=hzs[i][j]+(-(eys[i+1][j]-eys[i][j])+(exs[i][j+1]-exs[i][j]))*dtmds;
 
	}
/*    
        END SUBROUTINE ADV_HFIELD
!=============================================================
*/
      void RMS(int k)
      {
/*
c
c   Subroutine to calculate 
c
      INCLUDE 'xyfdtd.inc'
*/
       double z1,z2,p1,p2,p3,p4;


        for(j=2;j<=ny;j++)
	{        
	  for(i=istart;i<=ulim[rank];i++)
	  {
            z1=(ext[i][j]*ext[i][j]+ext[i-1][j]*ext[i-1][j])*.5;
            z2=(eyt[i][j]*eyt[i][j]+eyt[i][j-1]*eyt[i][j-1])*.5;
            erms2[i][j]=erms2[i][j]+(z1+z2)/(float)nperdt;
            exrms2[i][j]=exrms2[i][j]+(z1)/(float)(nperdt);
            eyrms2[i][j]=eyrms2[i][j]+(z2)/(float)(nperdt);
	  }
	}
        if(k==2)
	{
          for(i=1;i<=nx;i++)
          {
	    for(j=1;j<=ny;j++)
	    {
              ERMS[i][j] = sqrt(erms2[i][j]);  
              erms2[i][j]=0.0;
              exrms[i][j] = sqrt(exrms2[i][j]);
              exrms2[i][j]=0.0;
              eyrms[i][j] = sqrt(eyrms2[i][j]);  
              eyrms2[i][j]=0.0;
            } 
	  }
       }
    

      }

//!=============================================================

      void MR_MUR()
      {
/*
c
c   Subroutine to apply Mr. Murs absorbing boundary conditions
c

      INCLUDE 'xyfdtd.inc'
*/
        double csym;
        double aaa[2012],bbb[2012],ccc[2012],ddd[2012];
        double fff[2012],ggg[2012],hhh[2012],jjj[2012];


        csym=1.0;
        if(IABSOR==2) csym=0.0;

        for(j=1;j<=ny;j++)
		{
          eys[1][j]=eys1[2][j]+c1*(eys[2][j]-eys1[1][j]);          
          eys[nx][j]=eys1[nx-1][j]+c1*(eys[nx-1][j]-eys1[nx][j]);
        }

        for(i=1;i<=nx;i++)
		{
          exs[i][1]=exs1[i][2] +csym*c1*(exs[i][2]-exs1[i][1]);
          exs[i][ny]=exs1[i][ny-1] +csym*c1*(exs[i][ny-1]-exs1[i][ny]);
        }
        


          for(j=1;j<=ny;j++)
            for(i=1;i<=nx;i++)
	    	{
               exs1[i][j]=exs[i][j];
               eys1[i][j]=eys[i][j];
            }
	}
/*
!---
      END SUBROUTINE MR_MUR
!=============================================================


!-------------------------------------------------------
! Storing Intermediate Results in Resxx.res Files
!-------------------------------------------------------
*/
      void ECRIR()
      {
//      INCLUDE 'xyfdtd.inc'
      int necd,necu,ni,nj,necc,ni1,nj1;
      char fil[50];
      char chr1[2];
      char ch[50][3]={"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49"};
      char chr2[2];
      char chr3[3];
      char chr4[2];
      char tiret[40]="----------------------------------------";
//!-----------------------------
	strcpy(chr4,ch[rank+1]);

	if(nec>=1000) return;

	if(nec<10) 
	{	
	  strcpy(chr1,ch[nec+1]);
	  strcpy(fil,"ascii/re");
       	  strcat(fil,chr4);
	  strcat(fil,"00");
	  strcat(fil,chr1);
	  strcat(fil,".res");
	}
	if(nec>=10 && nec<100) 
	{

	  necd=nec/10;
	  necu=nec-necd*10;
	  strcpy(chr2,ch[necd+1]);
	  strcat(chr2,ch[necu+1]);
/*
chr2[0]=ch[necd+1][0];
chr2[1]=ch[necu+1][0];
*/
          strcpy(fil,"ascii/re");
       	  strcat(fil,chr4);
	  strcat(fil,"0");
	  strcat(fil,chr2);
	  strcat(fil,".res");

	}
	if(nec>=100) 
	{
	  necc=nec/100;
	  necd=(nec-necc*100)/10;
	  necu=nec-necc*100-necd*10;
/*chr3[2]=ch[necc+1][0];
chr3[1]=ch[necd+1][0];
chr3[0]=ch[necu+1][0];
*/
	  strcpy(chr3,ch[necc+1]);
	  strcat(chr3,ch[necd+1]);
	  strcat(chr3,ch[necu+1]);

printf("%s ",chr3);
          strcpy(fil,"ascii/re");
       	  strcat(fil,chr4);
	  strcat(fil,chr3);
	  strcat(fil,".res");
	}

       ni=nx;
       nj=ny-1;

//!-----------------------------
	fptr=fopen(fil,"w");
	fprintf(fptr,"%s",tiret);
	fprintf(fptr,"\n");

	fprintf(fptr,"%f %s",TIMD, "  s");
//!-----------------------------
	fprintf(fptr,"%s",tiret);	
        fprintf(fptr,"%s","  PLASMA DENSITY ");
        fprintf(fptr,"%d %d",ni,nj-1);
        for(j=1;j<=nj;j++)
	   for(i=1;i<=nx;i++)
   	       fprintf(fptr,"%f ",den[i][j]);
      
//!-----------------------------
	fprintf(fptr,"%s",tiret);
        fprintf(fptr,"%s","  RMS E FIELD");
        for(j=1;j<=nj;j++)
	   		for(i=1;i<=nx;i++)
	        	fprintf(fptr,"%f",ERMS[i][j]);
        
	fprintf(fptr,"%s",tiret);
//!-------------------------------
	fclose(fptr);
	}

void zero()
	{
//        include 'xyfdtd.inc'

               for(i=1;i<=nx;i++)
                {
		 for(j=1;j<=nyy;j++) 
		 {
			ERMS[i][j]=0.0;
			erms2[i][j]=0.0;
			exrms[i][j]=0.0;
			exrms2[i][j]=0.0;
			eyrms[i][j]=0.0;
			eyrms2[i][j]=0.0;
			den[i][j]=0.0;
			dent[i][j]=0.0;
			denp[i][j]=0.0;
			exi[i][j]=0.0;
			eyi[i][j]=0.0;
			exi1[i][j]=0.0;
			eyi1[i][j]=0.0;
			exs[i][j]=0.0;
			eys[i][j]=0.0;
			hzi[i][j]=0.0;
			hzs[i][j]=0.0;
			vx[i][j]=0.0;
			vy[i][j]=0.0;
			ext[i][j]=0.0;
			eyt[i][j]=0.0;
			exs1[i][j]=0.0;
			eys1[i][j]=0.0;
			exs2[i][j]=0.0;
			eys2[i][j]=0.0;
			frqio[i][j]=0.0;
			DIFFUSION[i][j]=0.0;
			SIO[i][j]=0.0;
			SIOT[i][j]=0.0;

		}
               }
		PARC=0.0;
       }

//c-----------------------------------------------------------


