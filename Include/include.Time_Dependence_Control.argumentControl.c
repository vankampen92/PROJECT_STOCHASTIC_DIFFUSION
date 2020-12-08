case 'D':
        if(argv[argcount][2]=='P')
        sscanf(argv[argcount+1],"%d", &TIME_DEPENDENT_PARAMETERS);
        else if(argv[argcount][2]=='C')
        sscanf(argv[argcount+1],"%d", &No_of_COVARIATES);
        else if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%d", &TYPE_0_PARAMETERS);
        else if(argv[argcount][2]=='1')
        sscanf(argv[argcount+1],"%d", &TYPE_1_PARAMETERS);
        else if(argv[argcount][2]=='2')
        sscanf(argv[argcount+1],"%d", &TYPE_2_PARAMETERS);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
      "-D...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'P':  //Dependent Parameters 
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%d", &dep_0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &dep_10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &dep_11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%d", &dep_12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"%d", &dep_13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%d", &dep_14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%d", &dep_15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%d", &dep_16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%d", &dep_17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%d", &dep_18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%d", &dep_19);
          else
          sscanf(argv[argcount+1],"%d", &dep_1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &dep_20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &dep_21);
          else
          sscanf(argv[argcount+1],"%d", &dep_2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%d", &dep_3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%d", &dep_4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%d", &dep_5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%d", &dep_6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%d", &dep_7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%d", &dep_8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%d", &dep_9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
      "-P...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'a':  //Pattern of time dependence 
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%d", &pat_0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &pat_10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &pat_11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%d", &pat_12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"-%d", &pat_13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%d", &pat_14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%d", &pat_15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%d", &pat_16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%d", &pat_17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%d", &pat_18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%d", &pat_19);
          else
          sscanf(argv[argcount+1],"%d", &pat_1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &pat_20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &pat_21);
          else
          sscanf(argv[argcount+1],"%d", &pat_2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%d", &pat_3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%d", &pat_4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%d", &pat_5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%d", &pat_6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%d", &pat_7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%d", &pat_8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%d", &pat_9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
      "-a...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

