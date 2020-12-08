typedef struct Generic_Root_Datainfo
{
  double x;
  double x_0;
  double x_1; 
  Parameter_Table * P;
  double (*Function)(const double, void *);
}Generic_Root_Data;
