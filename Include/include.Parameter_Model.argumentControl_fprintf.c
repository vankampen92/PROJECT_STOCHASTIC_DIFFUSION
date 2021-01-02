/* Human Parameters */

	    case 'H': /* Maximum and Minimum Transmission Rate */
              if(argv[argcount][2]=='M')         fprintf(fp,"-HM  %d  ",
						  Table->No_of_CELLS);       /* N */

	      else if(argv[argcount][2]=='u')    fprintf(fp,"-Hu  %g  ",
							 Table->Mu);       /* 0 */
		
	      else if(argv[argcount][2]=='N')    fprintf(fp,"-HN  %d  ",
							 Table->No_of_INDIVIDUALS);    /* 1 */

              else if(argv[argcount][2]=='X')    fprintf(fp,"-HX  %d  ",
							 Table->No_of_CELLS_X);       /* 0 */
		
	      else if(argv[argcount][2]=='Y')    fprintf(fp,"-HY  %d  ",
							 Table->No_of_CELLS_Y);    /* 1 */

	      else if(argv[argcount][2]=='S')    fprintf(fp,"-HS  %d  ",
							 Table->No_of_RESOURCES);    /* 1 */

              else if(argv[argcount][2]=='0')   fprintf(fp, "-H0 %g  ",
							&Lambda_R_0);                 /* 6 */

              else if(argv[argcount][2]=='1')   fprintf(fp, "-H1 %g  ",
							&Delta_R_0);                 /* 7 */

	      else if(argv[argcount][2]=='2')   fprintf(fp, "-H2 %g  ",
						        &Lambda_R_1);                /* 8 */
		
	      else if(argv[argcount][2]=='3')   fprintf(fp, "-H3 %g  ",
							&Delta_R_1);                /* 9 */

              else if(argv[argcount][2]=='K')   fprintf(fp, "-HK %d  ",
							&K_R);                      /* 10 */
	      
	      else {
		printf(" Error at reading input arguments: -  %s  \n", argv[argcount]);
		exit(0);
	      }
	      skip++;
	      break;

              #include <include.Type_of_Model.argumentControl_fprintf.c>
