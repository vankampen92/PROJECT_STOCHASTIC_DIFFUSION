 fprintf(fp, "\n");
 fprintf(fp, " Time Dependence Control Default Values:\n");

 fprintf(fp, " -DP %d\t No of time dependent parameters: %d\n", TIME_DEPENDENT_PARAMETERS, TIME_DEPENDENT_PARAMETERS);

 fprintf(fp, " -DC %d\t No of Covariates: %d\n", No_of_COVARIATES, No_of_COVARIATES);

 fprintf(fp, " -D0 %d\t Time dependent parameters of Type 0: %d\n", TYPE_0_PARAMETERS, TYPE_0_PARAMETERS);

 fprintf(fp, " -D1 %d\t Time dependent parameters of Type 1: %d\n", TYPE_1_PARAMETERS, TYPE_1_PARAMETERS);

 fprintf(fp, " -D2 %d\t Time dependent parameters of Type 2: %d\n", TYPE_2_PARAMETERS, TYPE_2_PARAMETERS);

 fprintf(fp, "-P0 0: Parameter Key (%s) 0: \n", Table->Symbol_Parameters[0]);
 fprintf(fp, " [Parameter: 0]: -a0 %d [time dependence pattern]\n", pat_0);

 fprintf(fp, "-P1 1: Parameter Key (%s) 1: \n", Table->Symbol_Parameters[1]);
 fprintf(fp, " [Parameter: 1]: -a1 %d [time dependence pattern]\n", pat_1);

 fprintf(fp, "-P2 2: Parameter Key (%s) 2: \n", Table->Symbol_Parameters[2]);
 fprintf(fp, " [Parameter: 2]: -a2 %d [time dependence pattern]\n", pat_2);

 fprintf(fp, "-P3 3: Parameter Key (%s) 3: \n", Table->Symbol_Parameters[3]);
 fprintf(fp, " [Parameter: 3]: -a3 %d [time dependence pattern]\n", pat_3);

 fprintf(fp, "-P4 4: Parameter Key (%s) 4: \n", Table->Symbol_Parameters[4]);
 fprintf(fp, " [Parameter: 4]: -a4 %d [time dependence pattern]\n", pat_4);

 fprintf(fp, "-P5 5: Parameter Key (%s) 5: \n", Table->Symbol_Parameters[5]);
 fprintf(fp, " [Parameter: 5]: -a5 %d [time dependence pattern]\n", pat_5);

 fprintf(fp, "-P6 6: Parameter Key (%s) 6: \n", Table->Symbol_Parameters[6]);
 fprintf(fp, " [Parameter: 6]: -a6 %d [time dependence pattern]\n", pat_6);

 fprintf(fp, " ... Time Dependence Control structure admits potentially more parameters\n"); 

