void ArgumentControl(int argc, char **argv);

void ArgumentControl_fprintf(FILE *, Parameter_Table *, int argc, char **argv);

void Parameters_from_Command_Line(FILE *fp, Parameter_Table *P);

void Parameters_ModelReport(char *File, Parameter_Table *P);

int M_O_D_E_L( Parameter_Table * Table );

int M_O_D_E_L___S_T_O( Parameter_Table * Table);
