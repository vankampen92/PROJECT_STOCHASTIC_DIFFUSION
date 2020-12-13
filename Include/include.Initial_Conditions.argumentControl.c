/* Scan Parameters */
 
/* Initial Conditions */
	case 'x':
             
		if(argv[argcount][2]=='n') 
		sscanf(argv[argcount+1], "%d",&TYPE_of_INITIAL_CONDITION);
		
		else if(argv[argcount][2]=='R') 
		  sscanf(argv[argcount+1], "%d",&RESCALING_INITIAL_TOTAL_POPULATION);
		
		else if(argv[argcount][2]=='N') 
		  sscanf(argv[argcount+1], "%lf",&INITIAL_TOTAL_POPULATION);
		
		else{
		  printf("xn -xN -xR are the only allowable keys.\n");
		  printf("Something goes very wrong (in an argumentsscanf function)\n");
		  exit(0);
		}
   
        skip++;
        break;

	case 'i':
              	
		if(argv[argcount][2]=='P')        sscanf(argv[argcount+1], "%d", &No_of_IC);
		
		else if(argv[argcount][2]=='0')    sscanf(argv[argcount+1], "%d", &IC_0);
		
		else{
		  printf("-iP -i0 are the only allowable keys.\n");
		  printf("Something goes very wrong (in an argumentsscanf function)\n");
		  exit(0);
		}
		
        skip++;
        break;

        case 'u': //Minimum values 
              
		
		if(argv[argcount][2]=='0')  sscanf(argv[argcount+1], "%lf", &IC_min_0);
                else{
		  printf("-u0 to -u9 are the only allowable keys.\n");
		  printf("Something goes very wrong (in an argumentfprintf function)\n");
		  exit(0);
		}

        skip++;
        break;

        case 'U':  //Maximum values

         
		
		  if(argv[argcount][2]=='0')         sscanf(argv[argcount+1], "%lf", &IC_MAX_0);
	     
		  else {
		    printf(" This TYPE_of_MODEL (%d) code is not defined.\n", TYPE_of_MODEL);
		    printf("Check input argument list\n");
		    exit(0);
		  }

        skip++;
        break;


