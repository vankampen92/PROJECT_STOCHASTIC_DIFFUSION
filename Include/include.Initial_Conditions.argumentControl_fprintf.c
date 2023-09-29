/* Scan Parameters */
 
/* Initial Conditions */
	case 'x':
              
		if(argv[argcount][2]=='n') 
		fprintf(fp, "-xn  %d  ",Table->TYPE_of_INITIAL_CONDITION);
		
		else if(argv[argcount][2]=='R') 
		  fprintf(fp, "-xR  %d  ",Table->RESCALING_INITIAL_TOTAL_POPULATION);
		
		else if(argv[argcount][2]=='N') 
		  fprintf(fp, "-xN  %g  ",Table->INITIAL_TOTAL_POPULATION);
		
		else{
		  printf("xn -xN -xR -x0 to -x9 are the only allowable keys.\n");
		  printf("Something goes very wrong (in an argumentfprintf function)\n");
		  exit(0);
		}
	                     
        skip++;
        break;

	case 'i':
		
		if(argv[argcount][2]=='P')        fprintf(fp, "-iP  %d  ", Table->No_of_IC);
		
		else if(argv[argcount][2]=='0')    fprintf(fp, "-i0  %d  ",Table->IC_0);
		
		else{
		  printf("-iP -i0 are the only allowable keys.\n");
		  printf("Something goes very wrong (in an argumentfprintf function)\n");
		  exit(0);
		}
		
        skip++;
        break;

        case 'u': //Minimum values 
              if (TYPE_of_MODEL == 0 ) {
		
		if(argv[argcount][2]=='0')         fprintf(fp, "-u0  %g  ", Table->IC_min_0);
		else{
		  printf("-u0 to -u9 are the only allowable keys.\n");
		  printf("Something goes very wrong (in an argumentfprintf function)\n");
		  exit(0);
		}
	      }
        skip++;
        break;

        case 'U':  //Maximum values

		  if(argv[argcount][2]=='0')  fprintf(fp, "-U0  %g  ", Table->IC_MAX_0);
		  
		  else{
		    printf("-U0 to -U9 are the only allowable keys.\n");
		    printf("Something goes very wrong (in an argumentfprintf function)\n");
		    exit(0);
		  }
		              
        skip++;
        break;


