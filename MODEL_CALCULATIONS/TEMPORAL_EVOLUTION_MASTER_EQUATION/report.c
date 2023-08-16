#include <MODEL.h>

#include "extern.h"
extern int HELP_INPUT_ARGUMENTS;

void Parameters_from_Command_Line(FILE *fp, Parameter_Table * Table)
{
  fprintf(fp, " M O D E L   P A R A M E T E R S  ( S E A R C H E A B L E ):\n");
  fprintf_Model_Parameters(fp, Table);

  fprintf(fp, " M O D E L   ( C O R E )   P A R A M E T E R S:\n");
#include <include.Parameter_Model.report.c>
  Print_Press_Key(1,0,"."); 

#if defined CPGPLOT_REPRESENTATION
#include <include.CPG.fprintPar.c>
#endif

#include <include.Type_of_Model.report.c>

#include <include.Trend_Control.report.c>

#include <include.Time_Control.report.c>

#include <include.Time_Dependence_Control.report.c>

#include <include.Initial_Conditions.report.c>

#include <include.Parameter_Space.report.c>

  fprintf_Output_Variables(fp, Table);

}

void Parameters_ModelReport(char *File, Parameter_Table * Table)
{
  FILE *fp;
  Parameters_from_Command_Line(stdout, Table);
  fp = fopen(File, "w");
  Parameters_from_Command_Line(fp, Table);
  fclose(fp);
}
