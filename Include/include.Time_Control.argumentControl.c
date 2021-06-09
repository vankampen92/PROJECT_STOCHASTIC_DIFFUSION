/* Time Control Definition */

	    case 't':
	      if(argv[argcount][2]=='n')
		sscanf(argv[argcount+1],"%d",&I_Time);
              else if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%lf",&Time_0);
              else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%lf",&Time_1);
              else if(argv[argcount][2]=='U')
		sscanf(argv[argcount+1],"%lf",&Time_Scale_Unit);
              else if(argv[argcount][2]=='R')
		sscanf(argv[argcount+1],"%d",&Realizations);
              else if(argv[argcount][2]=='3')
		sscanf(argv[argcount+1],"%lf",&Delta_T);
	      else if(argv[argcount][2]=='4')
                sscanf(argv[argcount+1],"%d",&TYPE_of_TIME_DEPENDENCE);
	      else if(argv[argcount][2]=='E')
		sscanf(argv[argcount+1],"%lf",&EPSILON);
	      else if(argv[argcount][2]=='D')
		sscanf(argv[argcount+1],"%d", &DISCARTING_EXTINCTIONS);
              else exit(0);
	      skip++;
	      break;

      
