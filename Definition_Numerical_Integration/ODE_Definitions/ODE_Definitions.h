int jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params);

int function (double t, const double y[], double dydt[], void *params);

void JACOBIAN_Matrix(gsl_matrix * m, const double *y, double t, int W, Parameter_Table * P);

void evaluating_RESTRICTED_JACOBIAN_Matrix(gsl_matrix * m, const double *y,
					   double t, Parameter_Table * P);
