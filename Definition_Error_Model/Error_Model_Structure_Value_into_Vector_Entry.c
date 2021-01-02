#include <MODEL.h>

void Parameter_Table_into_Vector_Entries_Error_Model ( Parameter_Table * P, gsl_vector * X,
						       int * Index,
						       int No_of_PARAMETERS,
						       int No_of_IC,
						       int No_of_Error_Parameters )
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_Error_Parameters; i++) {
    key = Index[i];
    value = Error_Model_into_Vector_Entry_Table( key, P );
    gsl_vector_set(X, i+No_of_PARAMETERS+No_of_IC, value);
  }
}

double Error_Model_into_Vector_Entry_Table(int j, Parameter_Table * Table)
{
  double value;
  
  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    switch(j)
      {
        case  0:  value = Table->Error_Parameter_0 ; /*  0 */
        break;
        case  1:  value = Table->Error_Parameter_1 ; /*  1 */
        break;
        default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 1\n");
        exit(0);
      }
  }
  else {
    switch(j)
      {
        case  0:  value = Table->Err_0 ; /*  0 */
        break;
        case  1:  value = Table->Err_1 ; /*  1 */
        break;
        case  2:  value = Table->Err_2 ; /*  2 */
        break;
        case  3:  value = Table->Err_3 ; /*  3 */
        break;
        case  4:  value = Table->Err_4 ; /*  4 */
        break;
        case  5:  value = Table->Err_5 ; /*  5 */
        break;
        case  6:  value = Table->Err_6 ; /*  6 */
        break;
        case  7:  value = Table->Err_7 ; /*  7 */
        break;
        case  8:  value = Table->Err_8 ; /*  8 */
        break;
        case  9:  value = Table->Err_9 ; /*  9 */
        break;
        case 10:  value = Table->Err_10 ; /*  10 */
        break;
        case 11:  value = Table->Err_11 ; /*  11 */
        break;
        case 12:  value = Table->Err_12 ; /*  12 */
        break;
        case 13:  value = Table->Err_13 ; /*  13 */
        break;
	case 14:  value = Table->Err_14 ; /*  11 */
        break;
        case 15:  value = Table->Err_15 ; /*  12 */
        break;
        case 16:  value = Table->Err_16 ; /*  13 */
        break;
        default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 13\n");
        exit(0);
      }
  }

  return(value); 
}

double Error_Model_into_Vector_Entry_Model(int j, Parameter_Model * Table)
{
  double value;
  
  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    switch(j)
      {
        case  0:  value = Table->Error_Parameter_0 ; /*  0 */
        break;
        case  1:  value = Table->Error_Parameter_1 ; /*  1 */
        break;
        default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 1\n");
        exit(0);
      }
  }
  else {
    switch(j)
      {
        case  0:  value = Table->Err_0 ; /*  0 */
        break;
        case  1:  value = Table->Err_1 ; /*  1 */
        break;
        case  2:  value = Table->Err_2 ; /*  2 */
        break;
        case  3:  value = Table->Err_3 ; /*  3 */
        break;
        case  4:  value = Table->Err_4 ; /*  4 */
        break;
        case  5:  value = Table->Err_5 ; /*  5 */
        break;
        case  6:  value = Table->Err_6 ; /*  6 */
        break;
        case  7:  value = Table->Err_7 ; /*  7 */
        break;
        case  8:  value = Table->Err_8 ; /*  8 */
        break;
        case  9:  value = Table->Err_9 ; /*  9 */
        break;
        case 10:  value = Table->Err_10 ; /*  10 */
        break;
        case 11:  value = Table->Err_11 ; /*  11 */
        break;
        case 12:  value = Table->Err_12 ; /*  12 */
        break;
        case 13:  value = Table->Err_13 ; /*  13 */
        break;
        case 14:  value = Table->Err_14 ; /*  11 */
        break;
        case 15:  value = Table->Err_15 ; /*  12 */
        break;
        case 16:  value = Table->Err_16 ; /*  13 */
        break;
        default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 13\n");
        exit(0);
      }
  }

  return(value); 
}
