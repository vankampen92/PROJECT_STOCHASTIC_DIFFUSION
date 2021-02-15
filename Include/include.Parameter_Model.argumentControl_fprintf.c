/* Human Parameters */

	    case 'H': /* Maximum and Minimum Transmission Rate */
              if(argv[argcount][2]=='M')         fprintf(fp,"-HM  %d  ",
						  Table->No_of_CELLS);       /* M */

	      else if(argv[argcount][2]=='u')    fprintf(fp,"-Hu  %g  ",
							 Table->Mu);       /* u */
		
	      else if(argv[argcount][2]=='N')    fprintf(fp,"-HN  %d  ",
							 Table->No_of_INDIVIDUALS);    /* N */

              else if(argv[argcount][2]=='X')    fprintf(fp,"-HX  %d  ",
							 Table->No_of_CELLS_X);       /* X */
		
	      else if(argv[argcount][2]=='Y')    fprintf(fp,"-HY  %d  ",
							 Table->No_of_CELLS_Y);       /* Y */

	      else if(argv[argcount][2]=='S')    fprintf(fp,"-HS  %d  ",
							 Table->No_of_RESOURCES);     /* S */

              else if(argv[argcount][2]=='0')   fprintf(fp, "-H0 %g  ",
							&Lambda_R_0);                 /* 6 */

              else if(argv[argcount][2]=='1') {

		if(argv[argcount][3]=='\0')     fprintf(fp, "-H1 %g  ",
							&Delta_R_0);                 /* 1 */

		else if(argv[argcount][3]=='0')   fprintf(fp, "-H10  %g",
						        &Nu_C_0);                   /* 10 */
		
		else if(argv[argcount][3]=='1')   fprintf(fp, "-H11  %g",
							&Chi_C_0);                  /* 11 */

		else if(argv[argcount][3]=='2')   fprintf(fp, "-H12  %g",
							 &Eta_C_0);                 /* 12 */

		else if(argv[argcount][3]=='2')   fprintf(fp, "-H13  %g",
							 &Mu_C);                 /* 12 */
		
		else {
		  printf(" Error in include.Parameter_Model.argumentControl_fprintf.c\n");
		  printf(" The program will exit\n");
		  exit(0); 
		}
	      }
	      else if(argv[argcount][2]=='2')   fprintf(fp, "-H2 %g  ",
						        &Lambda_R_1);                /* 2 */
		
	      else if(argv[argcount][2]=='3')   fprintf(fp, "-H3 %g  ",
							&Delta_R_1);                 /* 3 */

              else if(argv[argcount][2]=='K')   fprintf(fp, "-HK %d  ",
							&K_R);                      /*  K */

              else if(argv[argcount][2]=='4')   fprintf(fp, "-H4  %g",
							&Beta_R);                   /* 4 */

              else if(argv[argcount][2]=='5')   fprintf(fp, "-H5  %g",
							&Lambda_C_0);               /* 5 */

	      else if(argv[argcount][2]=='6')   fprintf(fp, "-H6  %g",
						        &Delta_C_0);                /* 6 */
		
	      else if(argv[argcount][2]=='7')   fprintf(fp, "-H7  %g",
							&Lambda_C_1);                /* 7 */

              else if(argv[argcount][2]=='8')   fprintf(fp, "-H8  %g",
							&Delta_C_1);                 /* 8 */

              else if(argv[argcount][2]=='9')   fprintf(fp, "-H9  %g",
							&Alpha_C_0);                 /* 9 */

	      
	      else {
		printf(" Error in include.Parameter_Model.argumentControl_fprintf.c\n");
		printf(" The program will exit\n");
		exit(0);
	      }
	      skip++;
	      break;

              #include <include.Type_of_Model.argumentControl_fprintf.c>
