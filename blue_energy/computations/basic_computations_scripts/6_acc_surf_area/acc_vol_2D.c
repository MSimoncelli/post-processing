/* This program calcultates the accessible volume, the accessible surface area 
and a matrice giving accessible positions for a 3D structure. */


/* To execute the program, you just have to launch it with an input file giving
the following information (program < inputfile)
first line : name of the file giving the positions of the atoms (> nom)
second line : total number of atoms (> natoms)
third line : diameter of the obstacles in bohrs (> Dobs)
fourth line : length of the box in the three directions in bohrs (> Lx,Ly)
fifth line : zmin and zmax (> zmin,zmax)
sixth line : number of bins in any direction (> nbins)
seventh line : diameter of the probe in bohrs (> Dprobe) */


/* Principle of the calculation : the positions given in input are considered as centers of obstacles,
the box is divided in little areas, each area is accessible (1) or not (0). */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define L 1000


char nom[50];
char nom2[50];
int natoms,nbins,***Acc,***Grad;
double zmin,zmax,dbins,vol,*X,*Y,*Z,Dobs,Dprobe,Lx,Ly,***curve,***curvesol,***curveliq,***normx,***normy,***normz;
double ***MAT1,***MAT2,***MAT3,***dAcc; /* To use gradient3D only */


void read();
void positions();
void filling_Acc();
void surface();
void xyzfile();
void gradient3D(double ***MAT_IN);
void curvature();


char **cmatrix(int nl, int nc);
int **imatrix(int nl, int nc);
double *dvector(int n);
int *ivector(int n);
double **dmatrix(int nl, int nc);
int ***tdimatrix(int X_SIZE, int Y_SIZE, int Z_SIZE);
double ***tddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE);


int main ( void )
{
 int i;

 printf("\nHello!\n");
 read();
 positions();
 filling_Acc();
 surface();
 //curvature();

 printf("\nCalculation done !\n\n");
 printf("\nThe accessible volume is %.0lf bohr^3. The total volume considering Lx, Ly and Lz is %.0lf bohr^3.\n\n",vol,Lx*Ly*(zmax-zmin));

 xyzfile();

 return 0;
}


/*Reads the input file and allocates the vectors and matrices*/
void read()
{
 scanf("%s",nom);      	 		 /* name of the positions' file */
 scanf("%d",&natoms);      		 /* number of atoms in this file */
 scanf("%lf",&Dobs);       		 /* diameter of the obstacles */
 scanf("%lf %lf",&Lx,&Ly);  	  	 /* length of the box in the three directions */
 scanf("%lf %lf",&zmin,&zmax);    	 
 scanf("%d",&nbins);  		         /* number of bins */
 dbins=nbins*1.0;
 scanf("%lf",&Dprobe);       		 /* diameter of the probe */

 X=dvector(natoms+1); Y=dvector(natoms+1); Z=dvector(natoms+1);
 Acc=tdimatrix(nbins+1,nbins+1,nbins+1);
 dAcc=tddmatrix(nbins+1,nbins+1,nbins+1);
 Grad=tdimatrix(nbins+1,nbins+1,nbins+1);
 normx=tddmatrix(nbins+1,nbins+1,nbins+1);
 normy=tddmatrix(nbins+1,nbins+1,nbins+1);
 normz=tddmatrix(nbins+1,nbins+1,nbins+1);
 curvesol=tddmatrix(nbins+1,nbins+1,nbins+1);
 curveliq=tddmatrix(nbins+1,nbins+1,nbins+1);
 curve=tddmatrix(nbins+1,nbins+1,nbins+1);
 MAT1=tddmatrix(nbins+1,nbins+1,nbins+1);
 MAT2=tddmatrix(nbins+1,nbins+1,nbins+1);
 MAT3=tddmatrix(nbins+1,nbins+1,nbins+1);

 printf("\nReading OK!\n\n");

}


/* Filling of X, Y, Z for all the atoms */
void positions()
{ 
 int i;
 FILE *in;
 char ligne[L];

 in=fopen(nom,"r");
	
 for(i=1;i<=natoms;i++)
 {
    fgets(ligne,L,in);
    sscanf(ligne,"%lf %lf %lf",&X[i],&Y[i],&Z[i]);
  }
  
 fclose(in);
 
 printf("\nPositions OK!\n\n");

}


/* Filling of the Acc matrix with 1 or 0 */
void filling_Acc()
{
 int i,j,k,n,access;
 double d,posx,posy,posz,dx,dy,dz;

 for(i=1;i<=nbins;i++)
  {
   printf("Bin %d in the x direction.\n",i);
   for(j=1;j<=nbins;j++)
    {
     for(k=1;k<=nbins;k++)
	{
	 n=1;
	 access=1;
	 posx=(i-0.5)*Lx/(nbins*1.0);
	 posy=(j-0.5)*Ly/(nbins*1.0);
	 posz=((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin;
	 while((access==1)&&(n<natoms))
		{
		  dx=posx-X[n];
		  dy=posy-Y[n];
		  dz=posz-Z[n];
		  /* Periodic boundary conditions */
	          dx/=Lx;
	          dy/=Ly;
	          /*dz/=Lz;*/
		  dx-=round(dx);
		  dy-=round(dy);
		  /*dz-=round(dz);*/
		  dx*=Lx;
		  dy*=Ly;
		  /*dz*=Lz;*/
		  d=dx*dx+dy*dy+dz*dz;
		  d=sqrt(d);
		  if(d<((Dprobe+Dobs)/2.0)) access=0;
		  n++;
		}
	 Acc[i][j][k]=access;
         dAcc[i][j][k]=access*1.0;
	 if(access==1) vol+=(Lx/(nbins*1.0))*(Ly/(nbins*1.0))*((zmax-zmin)/(nbins*1.0));
	}
    }
  }

}


/* Writes a xyz file with element F if accessible and O if not accessible */
void xyzfile(void)
{
 int i,j,k;
 
 FILE *out;

 out=fopen("access_vol.xyz","w");

 fprintf(out,"%d\n",nbins*nbins*nbins);
 fprintf(out,"pas 1\n");

 for(i=1;i<=nbins;i++)
  {
   for(j=1;j<=nbins;j++)
    {
     for(k=1;k<=nbins;k++)
      {
       if(Acc[i][j][k]==1) fprintf(out,"F %lf %lf %lf\n",(i-0.5)*Lx/(nbins*1.0),(j-0.5)*Ly/(nbins*1.0),((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin);
       if(Acc[i][j][k]==0) fprintf(out,"O %lf %lf %lf\n",(i-0.5)*Lx/(nbins*1.0),(j-0.5)*Ly/(nbins*1.0),((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin);
      }
    }
  } 
 
 fclose(out);
 
 out=fopen("access_surf.xyz","w");

 fprintf(out,"%d\n",nbins*nbins*nbins);
 fprintf(out,"pas 1\n");

 for(i=1;i<=nbins;i++)
  {
   for(j=1;j<=nbins;j++)
    {
     for(k=1;k<=nbins;k++)
      {
       if(Grad[i][j][k]==1) fprintf(out,"F %lf %lf %lf\n",(i-0.5)*Lx/(nbins*1.0),(j-0.5)*Ly/(nbins*1.0),((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin);
       if(Grad[i][j][k]==0) fprintf(out,"O %lf %lf %lf\n",(i-0.5)*Lx/(nbins*1.0),(j-0.5)*Ly/(nbins*1.0),((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin);
      }
    }
  } 
 
 fclose(out);
  
}


/* Calculates the accessible surface area from the Acc matrix */
/* For the method see Faraday Discussions, 144, 223-243 (2010) - equation 21.a. */
/* And "The Lattice Boltzmann Equation for Fluid Dynamics and Beyond", S. Succi - pp 67-69.*/
void surface()
{
 int i,j,k,nsurf;
 double ***grad_x,***grad_y,***grad_z;
 double Surf,Surfsol,Surfliq,norm_grad,cs2;
 FILE *out_surfpoints;

 Surf=0.0;
 Surfsol=0.0;
 Surfliq=0.0;
 nsurf=0;	/* Number of grid points which are on a surface */
 cs2=1.0/3.0; 

 grad_x=tddmatrix(nbins+1,nbins+1,nbins+1);
 grad_y=tddmatrix(nbins+1,nbins+1,nbins+1);
 grad_z=tddmatrix(nbins+1,nbins+1,nbins+1);
 
 gradient3D(dAcc);

 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             grad_x[i][j][k]=MAT1[i][j][k];
             grad_y[i][j][k]=MAT2[i][j][k];
             grad_z[i][j][k]=MAT3[i][j][k];
 	     norm_grad=sqrt(grad_x[i][j][k]*grad_x[i][j][k]+grad_y[i][j][k]*grad_y[i][j][k]+grad_z[i][j][k]*grad_z[i][j][k]);
 	     if(norm_grad!=0) {Grad[i][j][k]=1; nsurf++; normx[i][j][k]=grad_x[i][j][k]/norm_grad; normy[i][j][k]=grad_y[i][j][k]/norm_grad; normz[i][j][k]=grad_z[i][j][k]/norm_grad;} 
	     Surf+=(1.0/cs2)*norm_grad*Lx*Ly*(zmax-zmin)/(nbins*nbins*nbins*1.0);
	     if(Acc[i][j][k]==0) Surfsol+=(1.0/cs2)*norm_grad*Lx*Ly*(zmax-zmin)/(nbins*nbins*nbins*1.0);
	     if(Acc[i][j][k]==1) Surfliq+=(1.0/cs2)*norm_grad*Lx*Ly*(zmax-zmin)/(nbins*nbins*nbins*1.0);
            }
        }
    }
 
/* Write the surface points' positions in the surf_points.out file */
 out_surfpoints=fopen("surf_points.out","w");
 fprintf(out_surfpoints,"%d %lf\n",nsurf,Surf);
 for(i=1;i<=nbins;i++)
  {
   for(j=1;j<=nbins;j++)
    {
     for(k=1;k<=nbins;k++)
	{
  	 if(Grad[i][j][k]==1) fprintf(out_surfpoints,"%lf %lf %lf %lf %lf %lf\n",(i-0.5)*Lx/(nbins*1.0),(j-0.5)*Ly/(nbins*1.0),((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin,grad_x[i][j][k],grad_y[i][j][k],grad_z[i][j][k]); 
	}
    }
  }
 fclose(out_surfpoints);

 printf("\nThe accessible surface area is: %lf bohrs2.\n",Surf); 
 printf("\nThe accessible surface area (solid side) is: %lf bohrs2.\n",Surfsol); 
 printf("\nThe accessible surface area (liquid side) is: %lf bohrs2.\n",Surfliq); 

 /* Write the accessible points positions in the vol_points.out file */
 out_surfpoints=fopen("vol_points.out","w");
 fprintf(out_surfpoints,"%d %lf %lf %lf\n",nbins*nbins*nbins,Lx,Ly,zmax-zmin);
 for(i=1;i<=nbins;i++)
  {
   for(j=1;j<=nbins;j++)
    {
     for(k=1;k<=nbins;k++)
	{
  	 fprintf(out_surfpoints,"%lf %lf %lf %d\n",(i-0.5)*Lx/(nbins*1.0),(j-0.5)*Ly/(nbins*1.0),((k-0.5)*(zmax-zmin)/(nbins*1.0))+zmin,Acc[i][j][k]); 
	}
    }
  }
 fclose(out_surfpoints);

}


/* Calculate the curvature in each position of the grid */
void curvature()
{
 int i,j,k,cn,cm;
 double mean_curve,mean_curve_sol,mean_curve_liq;
 double cs2,w_center,w_1st,w_2nd,***temp_curve; 
 FILE *out; 

 temp_curve=tddmatrix(nbins+1,nbins+1,nbins+1);
 
 w_center=12.0/36.0;
 w_1st=2.0/36.0;
 w_2nd=1.0/36.0;
 cs2=1.0/3.0;

 mean_curve=0.0;
 mean_curve_sol=0.0;
 mean_curve_liq=0.0;

 gradient3D(normx);
 
 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             curve[i][j][k]+=MAT1[i][j][k];
             if(Acc[i][j][k]==0) curvesol[i][j][k]+=MAT1[i][j][k];
             if(Acc[i][j][k]==1) curveliq[i][j][k]+=MAT1[i][j][k];
            }
        }
    }             

 gradient3D(normy);
 
 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             curve[i][j][k]+=MAT2[i][j][k];
             if(Acc[i][j][k]==0) curvesol[i][j][k]+=MAT2[i][j][k];
             if(Acc[i][j][k]==1) curveliq[i][j][k]+=MAT2[i][j][k];
            }
        }
    }             
  
 gradient3D(normz);
 
 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             curve[i][j][k]+=MAT3[i][j][k];
             if(Acc[i][j][k]==0) curvesol[i][j][k]+=MAT3[i][j][k];
             if(Acc[i][j][k]==1) curveliq[i][j][k]+=MAT3[i][j][k];
            }
        }
    }             

 /* Average over first neighbours */
 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             temp_curve[i][j][k]=0.0;
             /* Central position */
             temp_curve[i][j][k]=w_center*curve[i][j][k];
             /* First neighbours */
             cn=i+1; {if(cn>nbins) cn=1;} temp_curve[i][j][k]+=w_1st*curve[cn][j][k]; 
             cn=i-1; {if(cn<1) cn=nbins;} temp_curve[i][j][k]+=w_1st*curve[cn][j][k]; 
             cn=j+1; {if(cn>nbins) cn=1;} temp_curve[i][j][k]+=w_1st*curve[i][cn][k]; 
             cn=j-1; {if(cn<1) cn=nbins;} temp_curve[i][j][k]+=w_1st*curve[i][cn][k];
             /* 2D conditions are different */
             cn=k+1; {if(cn>nbins) cn=nbins;} temp_curve[i][j][k]+=w_1st*curve[i][j][cn];
             cn=k-1; {if(cn<1) cn=1;} temp_curve[i][j][k]+=w_1st*curve[i][j][cn];
             /* 2nd neighbours */
             cn=i+1; cm=j+1; {if(cn>nbins) cn=1;} {if(cm>nbins) cm=1;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][cm][k];
             cn=i+1; cm=j-1; {if(cn>nbins) cn=1;} {if(cm<1) cm=nbins;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][cm][k];
             cn=i-1; cm=j+1; {if(cn<1) cn=nbins;} {if(cm>nbins) cm=1;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][cm][k];
             cn=i-1; cm=j-1; {if(cn<1) cn=nbins;} {if(cm<1) cm=nbins;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][cm][k];
             /* 2D conditions are different */
             cn=i+1; cm=k+1; {if(cn>nbins) cn=1;} {if(cm>nbins) cm=nbins;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][j][cm];
             cn=i+1; cm=k-1; {if(cn>nbins) cn=1;} {if(cm<1) cm=1;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][j][cm];
             cn=i-1; cm=k+1; {if(cn<1) cn=nbins;} {if(cm>nbins) cm=nbins;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][j][cm];
             cn=i-1; cm=k-1; {if(cn<1) cn=nbins;} {if(cm<1) cm=1;}
                             temp_curve[i][j][k]+=w_2nd*curve[cn][j][cm];
             cn=j+1; cm=k+1; {if(cn>nbins) cn=1;} {if(cm>nbins) cm=nbins;}
                             temp_curve[i][j][k]+=w_2nd*curve[i][cn][cm];
             cn=j+1; cm=k-1; {if(cn>nbins) cn=1;} {if(cm<1) cm=1;}
                             temp_curve[i][j][k]+=w_2nd*curve[i][cn][cm];
             cn=j-1; cm=k+1; {if(cn<1) cn=nbins;} {if(cm>nbins) cm=nbins;}
                             temp_curve[i][j][k]+=w_2nd*curve[i][cn][cm];
             cn=j-1; cm=k-1; {if(cn<1) cn=nbins;} {if(cm<1) cm=1;}
                             temp_curve[i][j][k]+=w_2nd*curve[i][cn][cm];
            }
        }
    }             

 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             curve[i][j][k]=temp_curve[i][j][k]; 
	    }
        }
    }

 out=fopen("curvature_points.out","w");
 
 for(i=1;i<=nbins;i++)
    {
     for(j=1;j<=nbins;j++)
        {
         for(k=1;k<=nbins;k++)
            {
             mean_curve+=curve[i][j][k]/2.0;
             mean_curve_sol+=curvesol[i][j][k]/2.0;
             mean_curve_liq+=curveliq[i][j][k]/2.0;
             fprintf(out,"%d %d %d %lf\n",i-1,j-1,k-1,curve[i][j][k]);
            }
        }
    }
 
 fclose(out);

 printf("The mean curvature is: %lf bohr-1.\n",mean_curve);
 printf("The mean curvature (solid side) is: %lf bohr-1.\n",mean_curve_sol);
 printf("The mean curvature (liquid side) is: %lf bohr-1.\n",mean_curve_liq);

}

void gradient3D(double ***MAT_IN)
{
 int i,j,k,cn,cm;
 double cs2,w_center,w_1st,w_2nd;
 
 w_center=12.0/36.0;
 w_1st=2.0/36.0;
 w_2nd=1.0/36.0;
 cs2=1.0/3.0;

 /* Calculation of the gradient in each position of the grid */
 for(i=1;i<=nbins;i++)
  {
   for(j=1;j<=nbins;j++)
    {
     for(k=1;k<=nbins;k++)
	{
	 MAT1[i][j][k]=0.0;
	 MAT2[i][j][k]=0.0;
	 MAT3[i][j][k]=0.0;
	/* Central position */
	/* Nothing to calculate */
	/* 1st neighbours*/
	 cn=i+1; {if(cn>nbins) cn=1;} MAT1[i][j][k]+=w_1st*(MAT_IN[cn][j][k]-MAT_IN[i][j][k])*dbins/Lx;
	 cn=i-1; {if(cn<1) cn=nbins;} MAT1[i][j][k]-=w_1st*(MAT_IN[cn][j][k]-MAT_IN[i][j][k])*dbins/Lx;
	 cn=j+1; {if(cn>nbins) cn=1;} MAT2[i][j][k]+=w_1st*(MAT_IN[i][cn][k]-MAT_IN[i][j][k])*dbins/Ly;
	 cn=j-1; {if(cn<1) cn=nbins;} MAT2[i][j][k]-=w_1st*(MAT_IN[i][cn][k]-MAT_IN[i][j][k])*dbins/Ly;
	/* 2D conditions are different */
	 cn=k+1; {if(cn>nbins) cn=nbins;} MAT3[i][j][k]+=w_1st*(MAT_IN[i][j][cn]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	 cn=k-1; {if(cn<1) cn=1;} MAT3[i][j][k]-=w_1st*(MAT_IN[i][j][cn]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	/* 2nd neighbours */
	 cn=i+1; {if(cn>nbins) cn=1;} cm=j+1; {if(cm>nbins) cm=1;} 
		MAT1[i][j][k]+=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT2[i][j][k]+=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Ly;	
	 cn=i+1; {if(cn>nbins) cn=1;} cm=j-1; {if(cm<1) cm=nbins;} 
		MAT1[i][j][k]+=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT2[i][j][k]-=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Ly;	
	/* 2D conditions are different */
	 cn=i+1; {if(cn>nbins) cn=1;} cm=k+1; {if(cm>nbins) cm=nbins;} 
		MAT1[i][j][k]+=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT3[i][j][k]+=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);	
	 cn=i+1; {if(cn>nbins) cn=1;} cm=k-1; {if(cm<1) cm=1;} 
		MAT1[i][j][k]+=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT3[i][j][k]-=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	 cn=i-1; {if(cn<1) cn=nbins;} cm=j+1; {if(cm>nbins) cm=1;} 
		MAT1[i][j][k]-=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT2[i][j][k]+=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Ly;
	 cn=i-1; {if(cn<1) cn=nbins;} cm=j-1; {if(cm<1) cm=nbins;} 
		MAT1[i][j][k]-=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT2[i][j][k]-=w_2nd*(MAT_IN[cn][cm][k]-MAT_IN[i][j][k])*dbins/Ly;
	/* 2D conditions are different */
	 cn=i-1; {if(cn<1) cn=nbins;} cm=k+1; {if(cm>nbins) cm=nbins;} 
		MAT1[i][j][k]-=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT3[i][j][k]+=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	 cn=i-1; {if(cn<1) cn=nbins;} cm=k-1; {if(cm<1) cm=1;} 
		MAT1[i][j][k]-=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/Lx; 
		MAT3[i][j][k]-=w_2nd*(MAT_IN[cn][j][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);	
	 cn=j+1; {if(cn>nbins) cn=1;} cm=k+1; {if(cm>nbins) cm=nbins;} 
		MAT2[i][j][k]+=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/Ly; 
		MAT3[i][j][k]+=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	 cn=j+1; {if(cn>nbins) cn=1;} cm=k-1; {if(cm<1) cm=1;} 
		MAT2[i][j][k]+=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/Ly; 
		MAT3[i][j][k]-=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	 cn=j-1; {if(cn<1) cn=nbins;} cm=k+1; {if(cm>nbins) cm=nbins;} 
		MAT2[i][j][k]-=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/Ly; 
		MAT3[i][j][k]+=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	 cn=j-1; {if(cn<1) cn=nbins;} cm=k-1; {if(cm<1) cm=1;} 
		MAT2[i][j][k]-=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/Ly; 
		MAT3[i][j][k]-=w_2nd*(MAT_IN[i][cn][cm]-MAT_IN[i][j][k])*dbins/(zmax-zmin);
	}
     }
  }

}



/********************************************************************************************/
/*Allocation dynamique*/

char **cmatrix(int nl, int nc)
{
  int i;
  char **m;
  m=(char **) malloc(nl*sizeof(char*));
  if (m) { m[0]=(char *) malloc(nl*nc*sizeof(char));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}


int **imatrix(int nl, int nc)
{
  int i;
  int **m;
  m=(int **) malloc(nl*sizeof(int*));
  if (m) { m[0]=(int *) malloc(nl*nc*sizeof(int));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

double *dvector(int n)
{
  return (double *) malloc(n*sizeof(double));
}

int *ivector(int n)
{
  return (int *) malloc(n*sizeof(int));
}

double **dmatrix(int nl, int nc)
{
  int i;
  double **m;
  m=(double **) malloc(nl*sizeof(double*));
  if (m) { m[0]=(double *) malloc(nl*nc*sizeof(double));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

int ***tdimatrix (int X_SIZE, int Y_SIZE, int Z_SIZE)
{
 int ***m;
 int i,j;

 m = (int ***)malloc(sizeof(int **) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (int **)malloc(sizeof(int *) * Y_SIZE);
  
    for (j = 0; j < Y_SIZE; j++)
       m[i][j] = (int *)malloc(sizeof(int) * Z_SIZE);
 }

 return m;
}

double ***tddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE)
{
 double ***m;
 int i,j;

 m = (double ***)malloc(sizeof(double **) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (double **)malloc(sizeof(double *) * Y_SIZE);
  
    for (j = 0; j < Y_SIZE; j++)
       m[i][j] = (double *)malloc(sizeof(double) * Z_SIZE);
 }

 return m;
}





