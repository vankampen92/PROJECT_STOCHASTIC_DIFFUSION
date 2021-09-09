#include <MODEL.h>

void Write_Parameter_Table ( Parameter_Table * P, int no_Par )
{
  double Par_Value;
  int i, k;
  
  printf(" Parameter table:\n" );
	 
  for ( i = 0; i < no_Par; i++){

    k = P->Index[i];

    Par_Value = AssignStructValue_to_VectorEntry (k, P);

    printf("%s \t %4.2f \n", P->Symbol_Parameters[k], Par_Value);
  }

  printf(" End writing parameter table!\n");
}

void write_Parameter_Table___RANGES___VALUES___LATEX ( char * File_Model_Parameters,
						       char * Type_of_Boundary_Parameter_Space,
						       Parameter_Table * P,
						       double * Par_MAX, double * Par_min,
						       int no_Par )
{
  double Par_Value;
  int i, k;
  FILE *fp = fopen( File_Model_Parameters, "w" );
  char Parameter_Name[100];
  char Parameter_Symbol[50];

  printf(" Writing latex parameter table in %s\n",
	 File_Model_Parameters);
  printf(" Parameter Space Definition Code (Boundary File): %s\n",
	 Type_of_Boundary_Parameter_Space);
  
  fprintf(fp, "\\begin{table}\n");
  // fprintf(fp, "\\noindent\n");
  fprintf(fp, "\\rowcolors{1}{green}{pink}\n"); 
  fprintf(fp, "\\centering\n"); 
  fprintf(fp, "\\begin{tabular}{p{5cm}lcccc}\n");
  fprintf(fp, "{\\bf Model Parameter} & {\\bf Symbol} & {\\bf Min. Value} & {\\bf Max. Value} & {\\bf Value}""\\""""\\""\n");
  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");

  for ( i = 0; i < no_Par; i++){

    k = P->Index[i];

    Par_Value = AssignStructValue_to_VectorEntry (k, P);

    AssignLabel_to_Model_Parameters__LATEX(k , Parameter_Name, P);

    AssignLabel_to_Model_Parameters__LATEX__SYMBOL(k , Parameter_Symbol, P);

    fprintf(fp, "%s & %s & %4.2f & %4.2f & %4.2f""\\""""\\""\n", 
	    Parameter_Name, Parameter_Symbol, Par_min[i], Par_MAX[i], Par_Value);
  }

  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");
  

  fprintf(fp, "\\end{tabular}\n");

  //fprintf(fp, "\\vspace{1cm}\n");
  //fprintf(fp, "\\noindent");
  fprintf(fp, "\\caption{Model Parameters for the AIDS-HIV transmission model. Model parameters boundary values defined in a boundary file.}\n");

  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  printf(" End writing latex parameter table!\n");
}

void write_Param_Table_PARAMETER_RANGES___LATEX ( char * File_Model_Parameters,
						  char * Type_of_Boundary_Parameter_Space,
						  Parameter_Table * P, 
						  double * Par_MAX, double * Par_min,
						  int no_Par )
{
  int i, k;
  FILE *fp = fopen( File_Model_Parameters, "w" );
  char Parameter_Name[100];
  char Parameter_Symbol[50];

  printf(" Writing latex parameter table in %s\n",
	 File_Model_Parameters);
  printf(" Parameter Space Definition Code (Boundary File): %s\n",
	 Type_of_Boundary_Parameter_Space);

  fprintf(fp, "\\begin{table}\n");
  // fprintf(fp, "\\noindent\n");
  fprintf(fp, "\\rowcolors{1}{green}{pink}\n"); 
  fprintf(fp, "\\centering\n"); 
  fprintf(fp, "\\begin{tabular}{lccc}\n");  
  fprintf(fp, "{\\bf Model Parameter} & {\\bf Symbol} & {\\bf Min. Value} & {\\bf Max. Value}""\\""""\\""\n"); 
  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");
  
  for ( i = 0; i < no_Par; i++){

    k = P->Index[i];
    
    AssignLabel_to_Model_Parameters__LATEX(k , Parameter_Name, P);

    AssignLabel_to_Model_Parameters__LATEX__SYMBOL(k , Parameter_Symbol, P);
    
    fprintf(fp, "%s & %s & %5.3f & %5.3f""\\""""\\""\n",
	    Parameter_Name, Parameter_Symbol, Par_min[i], Par_MAX[i]);
  }

  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");
  

  fprintf(fp, "\\end{tabular}\n");
  // fprintf(fp, "\\vspace{1cm}\n");
  // fprintf(fp, "\\noindent");

  fprintf(fp, "{\\caption Model Parameter Ranges for the AIDS-HIV transmission model. Model parameters boundary values defined in a boundary file.}\n");

  fprintf(fp, "\\end{table}\n");
  fclose(fp);
}

void write_Param_Table_AVE_STD___LATEX( char * File_Model_Parameters,
					Parameter_Table * P, 
					double * Par_AVE, double * Par_STD, int no_Par )
{
  int i, k;
  FILE *fp = fopen( File_Model_Parameters, "w" );
  char Parameter_Name[100];
  char Parameter_Symbol[50];

  fprintf(fp, "\\begin{table}\n");
  // fprintf(fp, "\\noindent\n");
  fprintf(fp, "\\rowcolors{1}{green}{pink}\n"); 
  fprintf(fp, "\\centering\n"); 
  fprintf(fp, "\\begin{tabular}{lccc}\n");  
  fprintf(fp, "{\\bf Model Parameter} & {\\bf Symbol} & {\\bf Average} & {\\bf Stardard Deviation}""\\""""\\""\n"); 
  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");
  
  for ( i = 0; i < no_Par; i++){

    k = P->Index[i];
    
    AssignLabel_to_Model_Parameters__LATEX(k , Parameter_Name, P);

    AssignLabel_to_Model_Parameters__LATEX__SYMBOL(k , Parameter_Symbol, P);
    
    fprintf(fp, "%s & %s & %5.3f & %5.3f""\\""""\\""\n", Parameter_Name, Parameter_Symbol, Par_AVE[i], Par_STD[i]);
  }

  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");
  

  fprintf(fp, "\\end{tabular}\n");
  //fprintf(fp, "\\vspace{1cm}\n");
  //fprintf(fp, "\\noindent");

  fprintf(fp, "\\caption{Model Parameter Average Values and Deviations for the AIDS-HIV transmission model.}\n");

  fprintf(fp, "\\end{table}\n");
  fclose(fp);
}

void write_Param_Table_AVE_CI___LATEX( char * File_Model_Parameters,
				       Parameter_Table * P, 
				       double * Par_AVE,
				       double ** Confidence_Interval,
				       int no_Par )
{
  int i, k;
  FILE *fp = fopen( File_Model_Parameters, "w" );
  char Parameter_Name[100];
  char Parameter_Symbol[50];

  fprintf(fp, "\\begin{table}\n");
  // fprintf(fp, "\\noindent\n");
  fprintf(fp, "\\rowcolors{1}{green}{pink}\n"); 
  fprintf(fp, "\\centering\n"); 
  fprintf(fp, "\\begin{tabular}{lccc}\n");  
  fprintf(fp, "{\\bf Model Parameter} & {\\bf Symbol} & {\\bf Average} & {\\bf Confidence Interval}""\\""""\\""\n"); 
  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");
  
  for ( i = 0; i < no_Par; i++){

    k = P->Index[i];
    
    AssignLabel_to_Model_Parameters__LATEX(k , Parameter_Name, P);

    AssignLabel_to_Model_Parameters__LATEX__SYMBOL(k , Parameter_Symbol, P);
    
    fprintf(fp, "%s & %s & %5.3f & (%5.3f, %5.3f)""\\""""\\""\n", Parameter_Name, Parameter_Symbol, 
	    Par_AVE[i], Confidence_Interval[i][0], Confidence_Interval[i][1]);
  }
  fprintf(fp, "\\hline");
  fprintf(fp, "\\hline\n");

  fprintf(fp, "\\end{tabular}\n");
  //fprintf(fp, "\\vspace{1cm}\n");
  //fprintf(fp, "\\noindent");

  fprintf(fp, "\\caption{Model Parameter Average Values and Confidence Intervals for the AIDS-HIV transmission model.}\n");

  fprintf(fp, "\\end{table}\n");
  fclose(fp);
}

