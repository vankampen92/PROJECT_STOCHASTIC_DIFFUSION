#include <MODEL.h>

void Vector_Entries_into_Parameter_Table_Error_Model ( const gsl_vector * X,
						       Parameter_Table * P,
						       int * Index,
						       int No_of_PARAMETERS,
						       int No_of_IC,
						       int No_of_Error_Parameters)
{
  int i;
  int key;
  double value;

  for( i=0; i<No_of_Error_Parameters; i++) {
    key = Index[i];
    value = gsl_vector_get(X, i+No_of_PARAMETERS+No_of_IC);
    Vector_Entry_into_Error_Model_Table( value, key, P );
  }
}

void Vector_Entry_into_Error_Model_Table(double value, int j, Parameter_Table * Table)
{
  
  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    switch(j)
      {
        case  0:  Table->Error_Parameter_0 = value; /*  0 */
        break;
        case  1:  Table->Error_Parameter_1 = value; /*  1 */
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
        case  0:  Table->Err_0 = value; /*  0 */
        break;
        case  1:  Table->Err_1 = value; /*  1 */
        break;
        case  2:  Table->Err_2 = value; /*  2 */
        break;
        case  3:  Table->Err_3 = value; /*  3 */
        break;
        case  4:  Table->Err_4 = value; /*  4 */
        break;
        case  5:  Table->Err_5 = value; /*  5 */
        break;
        case  6:  Table->Err_6 = value; /*  6 */
        break;
        case  7:  Table->Err_7 = value; /*  7 */
        break;
        case  8:  Table->Err_8 = value; /*  8 */
        break;
        case  9:  Table->Err_9 = value; /*  9 */
        break;
        case 10:  Table->Err_10 = value; /*  10 */
        break;
        case 11:  Table->Err_11 = value; /*  11 */
        break;
        case 12:  Table->Err_12 = value; /*  12 */
        break;
        case 13:  Table->Err_13 = value; /*  13 */
        break;
	case 14:  Table->Err_14 = value; /*  11 */
        break;
        case 15:  Table->Err_15 = value; /*  12 */
        break;
        case 16:  Table->Err_16 = value; /*  13 */
        break;
        default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 13\n");
        exit(0);
      }
  }
}

void Vector_Entry_into_Error_Model_Model(double value, int j, Parameter_Model * Table)
{
  
  if (j >= Table->OUTPUT_VARIABLES_GENUINE) {
    j -= Table->OUTPUT_VARIABLES_GENUINE;
    /* The first output variables are the model variables */
    switch(j)
      {
        case  0:  Table->Error_Parameter_0 = value; /*  0 */
        break;
        case  1:  Table->Error_Parameter_1 = value; /*  1 */
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
        case  0:  Table->Err_0 = value; /*  0 */
        break;
        case  1:  Table->Err_1 = value; /*  1 */
        break;
        case  2:  Table->Err_2 = value; /*  2 */
        break;
        case  3:  Table->Err_3 = value; /*  3 */
        break;
        case  4:  Table->Err_4 = value; /*  4 */
        break;
        case  5:  Table->Err_5 = value; /*  5 */
        break;
        case  6:  Table->Err_6 = value; /*  6 */
        break;
        case  7:  Table->Err_7 = value; /*  7 */
        break;
        case  8:  Table->Err_8 = value; /*  8 */
        break;
        case  9:  Table->Err_9 = value; /*  9 */
        break;
        case 10:  Table->Err_10 = value; /*  10 */
        break;
        case 11:  Table->Err_11 = value; /*  11 */
        break;
        case 12:  Table->Err_12 = value; /*  12 */
        break;
        case 13:  Table->Err_13 = value; /*  13 */
        break;
	case 14:  Table->Err_14 = value; /*  14 */
        break;
        case 15:  Table->Err_15 = value; /*  15 */
        break;
        case 16:  Table->Err_16 = value; /*  16 */
        break;
        default:
        printf(".... INVALID OUTPUT VARIABLE KEY [key = %d]\n", j);
        printf(".... The permited correspondences are:\n");
        printf(".... from 0 to 13\n");
        exit(0);
      }
  }
}
