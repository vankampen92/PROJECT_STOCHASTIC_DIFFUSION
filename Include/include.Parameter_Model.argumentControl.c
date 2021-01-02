/* Human Parameters */

	    case 'H': /* Maximum and Minimum Transmission Rate */
              if(argv[argcount][2]=='M')        sscanf(argv[argcount+1],"%d",
							&No_of_CELLS);                   /* 0 */

              else if(argv[argcount][2]=='X')   sscanf(argv[argcount+1],"%d",
							&No_of_CELLS_X);                 /* 1 */

              else if(argv[argcount][2]=='Y')   sscanf(argv[argcount+1],"%d",
							&No_of_CELLS_Y);                 /* 2 */

	      else if(argv[argcount][2]=='u')   sscanf(argv[argcount+1],"%lf",
						       &Mu);                             /* 3 */
		
	      else if(argv[argcount][2]=='N')   sscanf(argv[argcount+1],"%d",
							&No_of_INDIVIDUALS);             /* 4 */

              else if(argv[argcount][2]=='S')   sscanf(argv[argcount+1],"%d",
							&No_of_RESOURCES);                 /* 5 */

              else if(argv[argcount][2]=='0')   sscanf(argv[argcount+1],"%lf",
							&Lambda_R_0);                 /* 6 */

              else if(argv[argcount][2]=='1')   sscanf(argv[argcount+1],"%lf",
							&Delta_R_0);                 /* 7 */

	      else if(argv[argcount][2]=='2')   sscanf(argv[argcount+1],"%lf",
						        &Lambda_R_1);                /* 8 */
		
	      else if(argv[argcount][2]=='3')   sscanf(argv[argcount+1],"%lf",
							&Delta_R_1);                /* 9 */

              else if(argv[argcount][2]=='K')   sscanf(argv[argcount+1],"%d",
							&K_R);                      /* 10 */

	      else {
		printf(" Error at reading input arguments: -  %s  \n", argv[argcount]);
		exit(0);
	      }
	      skip++;
	      break;

              #include <include.Type_of_Model.argumentControl.c>
