#include <MODEL.h>

/* Output Variables: 
  
    Local Variables, Genuine Variables, The Rest of Model Variables 
  */
 
void AssignLabel_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  int k, i, n, m;
  char * p;

  char ** L = (char **)calloc(Table->LOCAL_STATE_VARIABLES, sizeof(char *));
  for(i = 0; i<Table->LOCAL_STATE_VARIABLES; i++) L[i] = (char *)calloc(5, sizeof(char));

  Defining_Output_Variables_Labels (Table, L);

  Label[0] = '\0';
  char * n_Label = (char *)calloc(10, sizeof(char) );

  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, Table);
  }
  else if (j < Table->LOCAL_STATE_VARIABLES ) {

    if(Table->TYPE_of_MODEL == 3) {
      k = j%Table->LOCAL_STATE_VARIABLES;
      n_Label[0] = '\0';
      if (k>0 && k<=Table->N_E) {
	      sprintf(n_Label, "%d", k-1);
	      i = 1;
      }
      else if (k>Table->N_E &&   k<=2*Table->N_E) {
	      sprintf(n_Label, "%d", k-1-Table->N_E);
	      i = 2;
      }
      else if (k>2*Table->N_E && k<=2*Table->N_E+Table->N_E*Table->N_E) {
	      i = 3;
	      n = (k-1-2*Table->N_E)%Table->N_E;
	      m = (k-1-2*Table->N_E)/Table->N_E;
	      p = strcat(Label, "[");
	      sprintf(n_Label, "%d", n);
	      p = strcat(Label, "]");
	      p = strcat(Label, "[");
	      sprintf(n_Label, "%d", m);
	      p = strcat(Label, "]");
      }
      else {
	      n_Label[0] = '\0';
	      i = 0;
      }
      p = strcat(Label, "n[");
      p = strcat(Label, L[i]);
      if (i>0) {
	      p = strcat(Label, "_");
	      p = strcat(Label, n_Label);
      }
      p = strcat(Label, "]");
    }
    else {
      p = strcat(Label, "n[");
      p = strcat(Label, L[j]);
      p = strcat(Label, "]");
    }
  }
  else if (j < Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->LOCAL_STATE_VARIABLES;

    switch(j)
      {
      case  0:
	      p = strcat(Label , "<n>");       /*  0: Density: No of Individuals per Cell */
        break;
      case  1:
	      p = strcat(Label , "s");         /*  1: SDV: No of Individuals per Cell */
        break;
      case  2:
	      p = strcat(Label , "N");         /*  2: Total No of Individuals */
        break;
      case  3:
	      p = strcat(Label , "n[R]");         /*  3: Resources */
        break;
      case  4:
	      p = strcat(Label , "n[A]");         /*  4: Total Free Consumers */
        break;
      case  5:
	      p = strcat(Label , "n[A_R]");       /*  5: Total Reproductive Consumers */
        break;
      case  6:
	      p = strcat(Label , "n[RA]");        /*  6: Total Handling Consumers */
        break;
      case  7:
	      p = strcat(Label , "n[ARA]");       /*  7: Total Triplets */
        break;
      case  8:
      	p = strcat(Label , "n[C]");       /*  8: Total Population Consumers */
        break;


      default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
      }
  }
  free(n_Label);
  for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++) free(L[i]);
  free(L);
}

void AssignLongLabel_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
{
  int k, i, n, m;
  char * p;

  char ** L = (char **)calloc(Table->LOCAL_STATE_VARIABLES, sizeof(char *));

  for(i = 0; i<Table->LOCAL_STATE_VARIABLES; i++) L[i] = (char *)calloc(5, sizeof(char));

  Defining_Output_Variables_Labels (Table, L);

  char * n_Label = (char *)calloc(10, sizeof(char) );
  Label[0] = '\0';

  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    AssignLabel_to_Model_Variables(j, Label, Table);
  }
  else if (j < Table->LOCAL_STATE_VARIABLES ) {
     if(Table->TYPE_of_MODEL == 3) {
      k = j%Table->LOCAL_STATE_VARIABLES;
      n_Label[0] = '\0';
      if (k>0 && k<=Table->N_E) {
	      sprintf(n_Label, "%d", k-1);
	      i = 1;
      }
      else if (k>Table->N_E &&   k<=2*Table->N_E) {
	      sprintf(n_Label, "%d", k-1-Table->N_E);
	      i = 2;
      }
      else if (k>2*Table->N_E && k<=2*Table->N_E+Table->N_E*Table->N_E) {
	      i = 3;
	      n = (k-1-2*Table->N_E)%Table->N_E;
	      m = (k-1-2*Table->N_E)/Table->N_E;
	      p = strcat(Label, "[");
	      sprintf(n_Label, "%d", n);
	      p = strcat(Label, "]");
	      p = strcat(Label, "[");
	      sprintf(n_Label, "%d", m);
	      p = strcat(Label, "]");
      }
      else {
	      n_Label[0] = '\0';
	      i = 0;
      }
      p = strcat(Label, "n[");
      p = strcat(Label, L[i]);
      if (i>0) {
	      p = strcat(Label, "_");
	      p = strcat(Label, n_Label);
      }
      p = strcat(Label, "]");
    }
    else {
      p = strcat(Label, "n[");
      p = strcat(Label, L[j]);
      p = strcat(Label, "]");
    }
  }
  else if (j < Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->LOCAL_STATE_VARIABLES;

    switch(j) {

    case  0:
      p = strcat(Label , "Average Density");                   /*  0: S */
      break;
    case  1:
      p = strcat(Label , "Standard Deviation of the Density"); /*  0: S */
      break;
    case  2:
      p = strcat(Label , "Total No of Individuals");           /*  0: S */
      break;
    case  3:
	p = strcat(Label , "Total Resources");         /*  3: Resources */
        break;
    case  4:
	p = strcat(Label , "Total Free Consumers");         /*  4: Total Free Consumers */
        break;
    case  5:
      p = strcat(Label , "Total Reproductive Consumers");  /*  5: Total Reproductive Consumers */
        break;
    case  6:
      p = strcat(Label , "Total Handling Consumers");      /*  6: Total Handling Consumers */
      break;
    case  7:
      p = strcat(Label , "Total Triplets");       /*  7: Total Triplets */
      break;
    case  8:
      p = strcat(Label , "Total Population Consumers");       /*  8: Total Population Consumers */
      break;


    default:
      printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
      printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 4\n");
        exit(0);
    }
  }
  free(n_Label);
  for(i=0; i<Table->LOCAL_STATE_VARIABLES; i++) free(L[i]);
  free(L);
}

void Defining_Output_Variables_Labels (Parameter_Table * Table, char ** L)
{
  int i; 
  char * p;
  char * q; 

  switch(Table->TYPE_of_MODEL)
    {
    case  0:
      p = strcat(L[0], "R");
      break;
    case  1:
      p = strcat(L[0], "R");

    break;
    case  2:
      p = strcat(L[0], "R");
      p = strcat(L[1], "A");
      p = strcat(L[2], "RA");
      p = strcat(L[3], "ARA");

      break;
    case  3:
      p = strcat(L[0], "R");
      p = strcat(L[1], "A");
      p = strcat(L[2], "RA");
      p = strcat(L[3], "ARA");

      break;
    case  4:
      p = strcat(L[0], "R");
      p = strcat(L[1], "A");

      break;
    case  5:
      p = strcat(L[0], "V");
      p = strcat(L[1], "R");
      p = strcat(L[2], "G");

      break;
    case  6:
      p = strcat(L[0], "V");
      p = strcat(L[1], "R");
      p = strcat(L[2], "G");

      break;
    case  7:
      p = strcat(L[0], "R");
      p = strcat(L[1], "A");
      p = strcat(L[2], "RA");

      break;
    case  8:
      p = strcat(L[0], "R");
      p = strcat(L[1], "A");
      p = strcat(L[2], "RA");
      p = strcat(L[3], "ARA");

      break;

    case  9:
      p = strcat(L[0], "A");
      p = strcat(L[1], "RA");

      break;

    case  10:
      p = strcat(L[0], "R");
      p = strcat(L[1], "A");
      p = strcat(L[2], "RA");

      break;

    case  11:
      p = strcat(L[0], "AC");
      p = strcat(L[1], "A");

        break;

    case  12:
      p = strcat(L[0], "A");

      break;

    case  13:
      p = strcat(L[0], "A");
      p = strcat(L[1], "RA");

      break;

    case  14:
      p = strcat(L[0], "A");
      p = strcat(L[1], "RA");
      p = strcat(L[1], "ARA");

      break;

    case  15:
      p = strcat(L[0], "RP");
      p = strcat(L[1], "R");
      p = strcat(L[2], "A");
      p = strcat(L[3], "RA");

      break;

    case  16:    /* DIFFUSION_HII_nD: Local Variables: Handling Consumers for each Resource Type */
      q = (char *)calloc(10, sizeof(char));
      for(i=0; i<Table->No_of_RESOURCES; i++){
        sprintf(q, "%d", i); 
        p = strcat(L[i], "C_");
        p = strcat(L[i], q);
      }
      free(q);
      break;

    case  17:
      p = strcat(L[0], "W");
      p = strcat(L[1], "Q");
      p = strcat(L[2], "F");
      p = strcat(L[3], "WF");

      break;

    case  18:
      p = strcat(L[0], "W");
      p = strcat(L[1], "Q");
      p = strcat(L[2], "F");
      p = strcat(L[3], "WF");

      break;
    
    case  19:
      p = strcat(L[0], "W");
      p = strcat(L[1], "Q");
      p = strcat(L[2], "F");
      p = strcat(L[3], "WF");

      break;
      
    default:
      printf(".... INVALID PARAMETER KEY (key = %d)\n", Table->TYPE_of_MODEL);
      printf(".... The permited correspondences are: 0 to 19\n");
      exit(0);
    }
}
