
case 'y': /* Type Control Definition */ 
	  if(argv[argcount][2]=='0')
		sscanf(argv[argcount+1],"%d",&TYPE_of_MODEL);
          else if(argv[argcount][2]=='1')
		sscanf(argv[argcount+1],"%d",&Growth_Function_Type);
          else if(argv[argcount][2]=='2')
		sscanf(argv[argcount+1],"%d",&TYPE_of_NETWORK);
          else {
		printf(" Error at reading input arguments: -  %s  \n", argv[argcount]);
		printf(" but you can only have:\n");
		printf(" -y0 [TYPE_of_MODEL] -y1 [Growth_Function_Type] -y2 [TYPE_of_NETWORK]\n");
		exit(0);
	  }
	  skip++;
	  break;

