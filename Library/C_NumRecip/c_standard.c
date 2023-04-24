#include <stdio.h>
#include <string.h>

#define N 355  /* Number of files */ 

main()
{
  FILE *fp;
  int j, n;
  char *s, *fitxer[N];
  char *dummy[N];

  fp = fopen("dades.txt","r");
  
  for(j=0; j<N; j++){
    s = fgets(s, N, fp);
    fitxer[j] = s;
    printf("%d\t%s\n",j,fitxer[j]);
  }
  
  for(j=0; j<N; j++){
    printf("%d\t%s\n",j, fitxer[j]);
    getchar();
  }

  for(j=0; j<N; j++){
    s = strcat("duconv -u ", fitxer[j]); 
    dummy[j] = strcat(fitxer[j],".lx");
    s = strcat(s," ");
    s = strcat(s,dummy[j]);
    printf("%s\n",s);
    n = system(s);
    getchar();
  }

  for(j=0; j<N; j++){
    s = strcat("gcc -c ", fitxer[j]); 
    printf("%s\n",s);
    n = system(s);
    getchar();
  }
}  

/*
void fitxer (char file[], char a[], int b, char c[])
{
  char *pFile;
  char num[12];

    pFile = strcat(file,a);

    inttochar(b,num);

    pFile = strcat(file,num);

    pFile = strcat(file,c);
} 
*/

void evol_03(int S, float v[], int n[], float x)
{
  /* This function writes the temporal evolution of the whole set of species */
  int i=0;
  FILE *fp;
  char file[12];

  for (i=1;i<=S;i++) {

    file[0]='\0';
    /* fitxer(file,"evol_",n[i],".dat"); */
    fp = fopen(file,"a+");
    fprintf(fp,"%f\t",x);
    fprintf(fp,"%5.4f\n",v[i]);
    fclose(fp);

  }
} 






