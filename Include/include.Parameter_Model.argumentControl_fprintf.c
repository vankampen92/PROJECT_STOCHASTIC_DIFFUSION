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
							 Table->No_of_SPECIES);    /* 1 */
	      
	      else {
		printf(" Error at reading input arguments: -  %s  \n", argv[argcount]);
		exit(0);
	      }
	      skip++;
	      break;

              #include <include.Type_of_Model.argumentControl_fprintf.c>
