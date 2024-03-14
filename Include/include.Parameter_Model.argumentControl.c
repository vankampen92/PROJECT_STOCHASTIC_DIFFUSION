/* Human Parameters */

case 'H': /* Maximum and Minimum Transmission Rate */
	if(argv[argcount][2]=='M')        sscanf(argv[argcount+1],"%d",
					&No_of_CELLS);                   /* M */

	else if(argv[argcount][2]=='X')   sscanf(argv[argcount+1],"%d",
					&No_of_CELLS_X);                 /* X */

	else if(argv[argcount][2]=='Y')   sscanf(argv[argcount+1],"%d",
					&No_of_CELLS_Y);                 /* Y */

	else if(argv[argcount][2]=='u') {
		if(argv[argcount][3]=='\0')   sscanf(argv[argcount+1],"%lf",
						&Mu);                        /* u */
		else if(argv[argcount][3]=='R')  sscanf(argv[argcount+1],"%lf",  
						&Mu);                       /*  u */
		else if(argv[argcount][3]=='C')  sscanf(argv[argcount+1],"%lf",   
						&Mu_C);                      /* uC */
		else if(argv[argcount][3]=='Q')  sscanf(argv[argcount+1],"%lf",
						&Lambda_C_1);                /* uQ */
		else {
			printf(" Error in include.Parameter_Model.argumentControl.c\n");
			printf(" Error at reading input arguments: %s  \n", argv[argcount]);
			printf(" The program will exit\n");
			exit(0);
		}
	}

	else if(argv[argcount][2]=='N')   sscanf(argv[argcount+1],"%d",
					&No_of_INDIVIDUALS);             /* N */

	else if(argv[argcount][2]=='S')   sscanf(argv[argcount+1],"%d",
					&No_of_RESOURCES);               /* S */

	else if(argv[argcount][2]=='K')   sscanf(argv[argcount+1],"%d", 
					&K_R);                      

	else if(argv[argcount][2]=='0')  sscanf(argv[argcount+1],"%lf",
					&Lambda_R_0);                    /* 0 */

	else if(argv[argcount][2]=='1')  {
		if(argv[argcount][3]=='\0')   sscanf(argv[argcount+1],"%lf",
						&Delta_R_0);                /* 1 */
		else if(argv[argcount][3]=='0')   sscanf(argv[argcount+1],"%lf",
						&Nu_C_0);                   /* 10 */
		else if(argv[argcount][3]=='1')   sscanf(argv[argcount+1],"%lf",
						&Chi_C_0);                  /* 11 */
		else if(argv[argcount][3]=='2')   sscanf(argv[argcount+1],"%lf",
						&Eta_C_0);                 /* 12 */
		else if(argv[argcount][3]=='3')   sscanf(argv[argcount+1],"%lf",
						&Mu_C);                    /* 13 */
		else if(argv[argcount][3]=='4')   sscanf(argv[argcount+1],"%d",
						&N_E);                    /* 14 */
		else if(argv[argcount][3]=='5')   sscanf(argv[argcount+1],"%d",
						&f);                    /* 15 */
		else if(argv[argcount][3]=='6')   sscanf(argv[argcount+1],"%d",
						&i_0);                    /* 16 */
		else if(argv[argcount][3]=='7')   sscanf(argv[argcount+1],"%lf",
						&Beta_C);                 /* 17 */
		else if(argv[argcount][3]=='8')   sscanf(argv[argcount+1],"%d",
						&k_E);                    /* 18 */
		else if(argv[argcount][3]=='9')   sscanf(argv[argcount+1],"%lf",
						&Theta_C);                /* 19 */
		else {
			printf(" Error in include.Parameter_Model.argumentControl.c\n");
			printf(" Error at reading input arguments: %s  \n", argv[argcount]);
			printf(" The program will exit\n");
			exit(0);
		}
	}

	else if(argv[argcount][2]=='2') {  
		if(argv[argcount][3]=='\0')   sscanf(argv[argcount+1],"%lf", 
					&Lambda_R_1);              /* 2 */
		else if(argv[argcount][3]=='0') sscanf(argv[argcount+1],"%lf", 
					&Eta_R);                  /* 20 */
		else {
			printf(" Error in include.Parameter_Model.argumentControl.c\n");
			printf(" Error at reading input arguments: %s  \n", argv[argcount]);
			printf(" The program will exit\n");
		exit(0);
		}
	}

	else if(argv[argcount][2]=='p')   {
		if(argv[argcount][3]=='1')  sscanf(argv[argcount+1],"%lf",        
					&p_1);                           /* p1 */
		else if(argv[argcount][3]=='2')  sscanf(argv[argcount+1],"%lf",
					&p_2);                          /* p2 */
		else {
			printf(" Error in include.Parameter_Model.argumentControl.c\n");
			printf(" Error at reading input arguments: %s  \n", argv[argcount]);
			printf(" The program will exit\n");
		exit(0);
		}
	}

	else if(argv[argcount][2]=='3')   sscanf(argv[argcount+1],"%lf",      
					&Delta_R_1);                

	else if(argv[argcount][2]=='4')   sscanf(argv[argcount+1],"%lf",
					&Beta_R);                   /* 4 */

	else if(argv[argcount][2]=='5')   sscanf(argv[argcount+1],"%lf",
					&Lambda_C_0);               /* 5 */

	else if(argv[argcount][2]=='6')   sscanf(argv[argcount+1],"%lf",
					&Delta_C_0);                /* 6 */

	else if(argv[argcount][2]=='7')   sscanf(argv[argcount+1],"%lf",
					&Lambda_C_1);                /* 7 */

	else if(argv[argcount][2]=='8')   sscanf(argv[argcount+1],"%lf",
					&Delta_C_1);                 /* 8 */

	else if(argv[argcount][2]=='9')   sscanf(argv[argcount+1],"%lf",
					&Alpha_C_0);                 /* 9 */

	else {
		printf(" Error in include.Parameter_Model.argumentControl.c\n");
		printf(" Error at reading input arguments: %s  \n", argv[argcount]);
		exit(0);
	}

	skip++;
break;

#include <include.Type_of_Model.argumentControl.c>
