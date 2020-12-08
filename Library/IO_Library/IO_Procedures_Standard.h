void Reading_Standard_Data_Matrix_from_File( char * File_Name,
					     double ** Matrix_Data,
					     int * No_of_Rows, int No_of_Columns,
					     int NAMES, char ** Name_of_Rows,
					     int TIMES, double * Time_Vector );

void Writing_Standard_Data_Matrix(double ** Matrix_Data,
				  int N_row, int N_column,
				  int NAMES, char ** Name_of_Rows,
				  int TIMES, double * Time_Vector);

void Writing_Standard_Data_Matrix_to_File(char * File_Name,
					  double ** Matrix_Data,
					  int N_row, int N_column,
					  int NAMES, char ** Name_of_Rows,
					  int TIMES, double * Time_Vector);
