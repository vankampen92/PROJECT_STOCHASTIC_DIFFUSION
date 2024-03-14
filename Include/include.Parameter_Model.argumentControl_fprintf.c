/* Human Parameters */

	    case 'H': /* Maximum and Minimum Transmission Rate */
              if(argv[argcount][2]=='M')         fprintf(fp,"-HM  %d  ",
						  Table->No_of_CELLS);       /* M */

	      else if(argv[argcount][2]=='u') {
		if(argv[argcount][3]=='\0')     fprintf(fp,"-Hu  %g  ", Table->Mu);       /* u */
		
		else if(argv[argcount][3]=='R') fprintf(fp, "-HuR %g  ", Table->Mu);    /*  uR */
						       
		else if(argv[argcount][3]=='C') fprintf(fp, "-HuC %g  ", Table->Mu_C);  /*  uC */
		
		else if(argv[argcount][3]=='Q') fprintf(fp, "-HuQ %g  ", Table->Lambda_C_1);  /*  uQ */
		
		else {
		  printf(" Error in include.Parameter_Model.argumentControl.c\n");
		  printf(" Error at reading input arguments: %s  \n", argv[argcount]);
		  printf(" The program will exit\n");
		  exit(0);
		}
	      }
	      else if(argv[argcount][2]=='N')    fprintf(fp,"-HN  %d  ",
							 Table->No_of_INDIVIDUALS);    /* N */

              else if(argv[argcount][2]=='X')    fprintf(fp,"-HX  %d  ",
							 Table->No_of_CELLS_X);       /* X */
		
	      else if(argv[argcount][2]=='Y')    fprintf(fp,"-HY  %d  ",
							 Table->No_of_CELLS_Y);       /* Y */

	      else if(argv[argcount][2]=='S')    fprintf(fp,"-HS  %d  ",
							 Table->No_of_RESOURCES);     /* S */

	      else if(argv[argcount][2]=='K')   fprintf(fp, "-HK %d  ",
							Table->K_R);                      /*  K */

              else if(argv[argcount][2]=='0')   fprintf(fp, "-H0 %g  ",
							Table->Lambda_R_0);                 /* 6 */

              else if(argv[argcount][2]=='1') {

		if(argv[argcount][3]=='\0')     fprintf(fp, "-H1 %g  ",
							Table->Delta_R_0);                 /* 1 */

		else if(argv[argcount][3]=='0')   fprintf(fp, "-H10  %g",
						        Table->Nu_C_0);                   /* 10 */
		
		else if(argv[argcount][3]=='1')   fprintf(fp, "-H11  %g",
							Table->Chi_C_0);                  /* 11 */

		else if(argv[argcount][3]=='2')   fprintf(fp, "-H12  %g",
							 Table->Eta_C_0);                 /* 12 */

		else if(argv[argcount][3]=='2')   fprintf(fp, "-H13  %g",
							 Table->Mu_C);                 /* 13 */

		else if(argv[argcount][3]=='4')   fprintf(fp, "-H14  %d",
							 Table->N_E);                    /* 14 */

		else if(argv[argcount][3]=='5')   fprintf(fp, "-H15  %d",
							 Table->f);                    /* 15 */

		else if(argv[argcount][3]=='6')   fprintf(fp, "-H16  %d",
							 Table->i_0);                    /* 16 */

		else if(argv[argcount][3]=='7')   fprintf(fp, "-H17  %g",
							 Table->Beta_C);                 /* 17 */

		else if(argv[argcount][3]=='8')   fprintf(fp, "-H18  %d",
							 Table->k_E);                    /* 18 */

		else if(argv[argcount][3]=='9')   fprintf(fp, "-H19  %g",
							 Table->Theta_C);                /* 19 */
		
		else {
		  printf(" Error in include.Parameter_Model.argumentControl_fprintf.c\n");
		  printf(" The program will exit\n");
		  exit(0); 
		}
	      }

              else if(argv[argcount][2]=='2') { 
		if(argv[argcount][3]=='\0')   
		  fprintf(fp, "-H2 %g  ", Table->Lambda_R_1);                /* 2 */

		else if(argv[argcount][3]=='0')
		  fprintf(fp, "-H20 %g ", Table->Eta_R);                     /* 20 */

		else {
		  printf(" Error in include.Parameter_Model.argumentControl.c\n");
		  printf(" Error at reading input arguments: %s  \n", argv[argcount]);
		  printf(" The program will exit\n");
		  exit(0);
		}
		
              } 
	      else if(argv[argcount][2]=='p')   {
		if(argv[argcount][3]=='1')  fprintf(fp, "-Hp1  %lf",                     /* p1 */
						   Table->p_1); 
						        
		else if(argv[argcount][3]=='2')  fprintf(fp, "-Hp2  %lf",
						   Table->p_2);                          /* p2 */
		else {
		  printf(" Error in include.Parameter_Model.argumentControl.c\n");
		  printf(" Error at reading input arguments: %s  \n", argv[argcount]);
		  printf(" The program will exit\n");
		  exit(0);
		}

	      }
		
	      else if(argv[argcount][2]=='3')   fprintf(fp, "-H3 %g  ",
							Table->Delta_R_1);                 /* 3 */

              else if(argv[argcount][2]=='4')   fprintf(fp, "-H4  %g",
							Table->Beta_R);                   /* 4 */

              else if(argv[argcount][2]=='5')   fprintf(fp, "-H5  %g",
							Table->Lambda_C_0);               /* 5 */

	      else if(argv[argcount][2]=='6')   fprintf(fp, "-H6  %g",
						        Table->Delta_C_0);                /* 6 */
		
	      else if(argv[argcount][2]=='7')   fprintf(fp, "-H7  %g",
							Table->Lambda_C_1);                /* 7 */

              else if(argv[argcount][2]=='8')   fprintf(fp, "-H8  %g",
							Table->Delta_C_1);                 /* 8 */

              else if(argv[argcount][2]=='9')   fprintf(fp, "-H9  %g",
							Table->Alpha_C_0);                 /* 9 */

	      
	      else {
		printf(" Error in include.Parameter_Model.argumentControl_fprintf.c\n");
		printf(" The program will exit\n");
		exit(0);
	      }
	      skip++;
	      break;

              #include <include.Type_of_Model.argumentControl_fprintf.c>
