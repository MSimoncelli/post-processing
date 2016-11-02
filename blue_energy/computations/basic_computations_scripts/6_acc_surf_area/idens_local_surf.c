/* This program calcultates ionic densities as a function of the distance to a surface
from the positions' file from Paul Madden's simulation code (cartesian coordinates). */

/* To execute the program, you just have to launch it with an input file giving
the following information (program < inputfile)
first line : name of the positions's file 
second line : number of configurations in the positions' file
third line : name of the surface positions' file
fourth line : length of the box in the 3 cartesian directions
fifth line : number of different species and number of species considered in the calculation
following lines : number of ions for each species 
following line : maximum length explored and number of boxes into which this length will be divided 
following line : limits in the z direction zmin,zmax. */

/* The code is written without considering the periodic boundary conditions, so be careful or modify the code ! */ 


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define L 1000
#define pi 3.14159265


char nom[200],nom2[200];
int nions,nions_zap,nconfigs,nbox,*nspecies,diffspecies,diffmove,*count;
double Lx,Ly,Lz,**pos,**boxes,maxlength,zmin,zmax;


void read();
void write();
int find(int j);
void normalize();


char **cmatrix(int nl, int nc);
int **imatrix(int nl, int nc);
double *dvector(int n);
int *ivector(int n);
double **dmatrix(int nl, int nc);


int main ( void )
{

 printf("Hello\n");
 read();
 write();
 printf("The output file is surface_density.out. \n");
 
 return 0;
}


/*Reads the input file and allocates the vectors and matrices*/
void read(){
  int i;
  scanf("%s",nom);					/* name of the positions' file */
  scanf("%d",&nconfigs);	      				/* number of configurations in the positions's file */
  scanf("%s",nom2);					/* name of the surface positions' file */
  scanf("%lf %lf %lf",&Lx,&Ly,&Lz);	   		/* lengths of the simulation box */
  scanf("%d %d",&diffspecies,&diffmove);  		/* number of different species */	
  nspecies=ivector(diffspecies+1);  //allocate a vector of size diffspecies+1
  count=ivector(diffspecies+1); //allocate a vector of size diffspecies+1
  for(i=1;i<=diffspecies;i++) scanf("%d",&nspecies[i]); 	/* number of ions for each species */
  scanf("%lf %d",&maxlength,&nbox);			/* number of boxes in the chosen direction */
  scanf("%lf %lf",&zmin,&zmax);	
  nions=0; nions_zap=0;
   //my modification: considers (O,H1,H2)  as a single specie to discard and the ions...the same for (C1,C2,P)
  nions=nspecies[2]+nspecies[3];
  nions_zap=nspecies[1]+nspecies[4]; 
  //for(i=1;i<=diffmove;i++) nions+=nspecies[i];			/* Calculation of the total number of ions */
  //for(i=diffmove+1;i<=diffspecies;i++) nions_zap+=nspecies[i];	/* Calculation of the number of ions not considered */
  nspecies[0]=nions; //whats the meaning of this??->seems NOT useful at all!!!
  printf("%d %d\n",nions,nions_zap);
  boxes=dmatrix(diffmove+1,nbox+1);	/* Boxes to fill in with the number of ions in each box */
  pos=dmatrix(4,nions+1);
  printf("Reading : OK !\n");
}

/* Calculates the number densities and output them */
void write(){
  FILE *in,*in2,*out;
  char ligne[L],ligne2[L];
  int i,j,k,nsurf,ibox,ipoint_surf,zug;
  double surf_area,dx,dy,dz,*surfposx,*surfposy,*surfposz,posx,posy,posz,dist,mindist;
  double norm_grad,*surfgradx,*surfgrady,*surfgradz,scal;

  /* Positions of the surface points */
  in2=fopen(nom2,"r");
  fgets(ligne2,L,in2);
  sscanf(ligne2,"%d %lf",&nsurf,&surf_area); //ok, the first line of the file surf_points.out contains the number of surface poinst and the surface area
  surfposx=dvector(nsurf+1); surfposy=dvector(nsurf+1); surfposz=dvector(nsurf+1);
  surfgradx=dvector(nsurf+1); surfgrady=dvector(nsurf+1); surfgradz=dvector(nsurf+1);
  //read the coordinate of each surface point and the normal versor to the surface in that point 
  for(k=1;k<=nsurf;k++){
    fgets(ligne2,L,in2);	 
    sscanf(ligne2,"%lf %lf %lf %lf %lf %lf",&surfposx[k],&surfposy[k],&surfposz[k],&surfgradx[k],&surfgrady[k],&surfgradz[k]);
  }
  fclose(in2);
 
  /* Filling of boxes= make an histrogram?!? */
  in=fopen(nom,"r");
  //skip the initial O,H1,H2 atoms.
  for(zug=0;zug<nspecies[1];zug++){
  	fgets(ligne,L,in);
  }
  for(i=1;i<=nconfigs;i++){ //loop over configs
    if(i%10==0) printf("Config %d over %d\n",i,nconfigs); //just a matter of UI
    for(j=1;j<=nions;j++){
      fgets(ligne,L,in); //Reads characters from stream and stores them as a C string into str until (L-1) characters have been read
      sscanf(ligne,"%lf %lf %lf",&posx,&posy,&posz);//extract the 3 coordinated of the ION!
      mindist=maxlength*1.1; 
      ipoint_surf=1;
      if((posz>=zmin)&&(posz<=zmax)){
        //look for the point of the surface which is closest to that ion!!!!!
        for(k=1;k<=nsurf;k++){
          dx=posx-surfposx[k]; dy=posy-surfposy[k]; dz=posz-surfposz[k];
          /* Periodic boundary conditions */
          dx/=Lx; dy/=Ly; dz/=Lz;
          //here you do not multiply by 2 and use the function round instead of the truncation to int!
          dx-=round(dx); dy-=round(dy); dz-=round(dz);
          dx*=Lx; dy*=Ly; dz*=Lz;
          dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
          if(dist<=mindist) {mindist=dist; ipoint_surf=k;}
        }
        /* Projection of the vector from the surface to the ion on the surface normal */
        //shortest distance between the point and the surface!
        dx=posx-surfposx[ipoint_surf]; dy=posy-surfposy[ipoint_surf]; dz=posz-surfposz[ipoint_surf];
        //scalar product!
        scal=dx*surfgradx[ipoint_surf]+dy*surfgrady[ipoint_surf]+dz*surfgradz[ipoint_surf];
        //strange, this should be the norm of the orthogonal versor i.e. should be 1
        norm_grad=surfgradx[ipoint_surf]*surfgradx[ipoint_surf]+surfgrady[ipoint_surf]*surfgrady[ipoint_surf]+surfgradz[ipoint_surf]*surfgradz[ipoint_surf];
        norm_grad=sqrt(norm_grad);
        //just normalize the gradient-> in the end this is eqwuibvalent to do the scalar product with the normalized gradient.
        mindist=fabs(scal/norm_grad);
        ibox=0;
        //find the box index in the hystogrma which correspond to the mindist!
        while((mindist>(ibox*maxlength/nbox))&&(ibox<=nbox+1)) ibox++; 

        if((ibox<=nbox)&&(ibox>=1)){
          //boxes[0][] does not distinguish between ionic species.
          //boxes[1][] is the histogram relative to the 1st type off ion
          //boxes[2][] is the hystogrm relative to the 2nd type of ion
          boxes[0][ibox]++; boxes[find(j)][ibox]++; count[find(j)]++; count[0]++;
        }
      }
    }
    for(j=1;j<=nions_zap;j++){
      fgets(ligne,L,in);
    }
  } 
  fclose(in);

 /* Average over the number of configurations, and obtention of the ionic density dividing by the volume of each box */
 out=fopen("surface_density.out","w");
  for(k=1;k<=nbox;k++){
    for(j=0;j<=diffmove;j++){
      boxes[j][k]/=(surf_area*nconfigs*maxlength/nbox); 
      if(j==0) fprintf(out,"%e	%e",k*maxlength/nbox,boxes[j][k]);
      if((j!=0)&&(j!=diffmove)) fprintf(out,"	%e",boxes[j][k]);
      if(j==diffmove) fprintf(out,"	%e\n",boxes[j][k]);
    }
  }
 fclose(out);

 /* Average over the counted ions, and obtention of the ionic probability */
 out=fopen("surface_probability.out","w");
 for(k=1;k<=nbox;k++)
    {
     for(j=0;j<=diffmove;j++)
	{
  //since in the above function you divided, now you multiply to recover the original case!!!
	 boxes[j][k]*=(surf_area*nconfigs*maxlength/nbox);
	 boxes[j][k]/=(count[j]*0.01);
	 if(j==0) fprintf(out,"%e	%e",k*maxlength/nbox,boxes[j][k]);
	 if((j!=0)&&(j!=diffmove)) fprintf(out,"	%e",boxes[j][k]);
	 if(j==diffmove) fprintf(out,"	%e\n",boxes[j][k]);
	}
    }
 fclose(out);

}


/* Find the species ion j belongs to */
int find(int j)
{
 int i;
 //int i, sum;
 //i=1; sum=nspecies[i];
 //while(j>sum) {i++; sum+=nspecies[i];}
 
  if (j<=(nions/2)){
    i=1;
  }
  else{
    i=2;
  }
 return i;
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
