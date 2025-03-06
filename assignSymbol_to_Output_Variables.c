#include <MODEL.h>

void AssignSymbol_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
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
      p = strcat(Label , "x");         /*  0 */
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
	    p = strcat(Label , "n[A]");             /*  4: Total Free Consumers */
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
      p = strcat(Label , "n[C]");       /*  8: Total Populaiton Consumers  */
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

void AssignCPGPLOT_Symbol_to_Output_Variables(int j, char * Label, Parameter_Table * Table)
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
    AssignSymbol_to_Model_Variables(j, Label, Table);
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
	p = strcat(Label, "\\d");
	p = strcat(Label, n_Label);
	p = strcat(Label, "\\u");
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
      p = strcat(Label , "x");       /*  0 */
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
      p = strcat(Label , "n[C]");       /*  8: Total Populaiton Consumers  */
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
