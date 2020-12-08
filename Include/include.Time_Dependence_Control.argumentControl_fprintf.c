case 'D':
        if(argv[argcount][2]=='P')
        fprintf(fp,"-DP  %d  ", Table->TDC->TIME_DEPENDENT_PARAMETERS);
        else if(argv[argcount][2]=='C')
        fprintf(fp,"-DC  %d  ", Table->TDC->No_of_COVARIATES);
        else if(argv[argcount][2]=='0')
        fprintf(fp,"-D0  %d  ", Table->TDC->TYPE_0_PARAMETERS);
        else if(argv[argcount][2]=='1')
        fprintf(fp,"-D1  %d  ", Table->TDC->TYPE_1_PARAMETERS);
        else if(argv[argcount][2]=='2')
        fprintf(fp,"-D2  %d  ", Table->TDC->TYPE_2_PARAMETERS);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-D...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'P':  //Dependent Parameters 
        if(argv[argcount][2]=='0')
        fprintf(fp,"-P0  %d  ", dep_0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-P10  %d  ", dep_10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-P11  %d  ", dep_11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-P12  %d  ", dep_12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-P13  %d  ", dep_13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-P14  %d  ", dep_14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-P15  %d  ", dep_15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-P16  %d  ", dep_16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-P17  %d  ", dep_17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-P18  %d  ", dep_18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-P19  %d  ", dep_19);
          else
          fprintf(fp,"-P1  %d  ", dep_1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-P20  %d  ", dep_20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-P21  %d  ", dep_21);
          else
          fprintf(fp,"-P2  %d  ", dep_2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-P3  %d  ", dep_3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-P4  %d  ", dep_4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-P5  %d  ", dep_5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-P6  %d  ", dep_6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-P7  %d  ", dep_7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-P8  %d  ", dep_8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-P9  %d  ", dep_9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-P...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'a':  //Pattern of time dependence 
        if(argv[argcount][2]=='0')
        fprintf(fp,"-a0  %d  ", pat_0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-a10  %d  ", pat_10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-a11  %d  ", pat_11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-a12  %d  ", pat_12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-13  %d  ", pat_13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-a14  %d  ", pat_14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-a15  %d  ", pat_15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-a16  %d  ", pat_16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-a17  %d  ", pat_17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-a18  %d  ", pat_18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-a19  %d  ", pat_19);
          else
          fprintf(fp,"-a1  %d  ", pat_1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-a20  %d  ", pat_20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-a21  %d  ", pat_21);
          else
          fprintf(fp,"-a2  %d  ", pat_2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-a3  %d  ", pat_3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-a4  %d  ", pat_4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-a5  %d  ", pat_5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-a6  %d  ", pat_6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-a7  %d  ", pat_7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-a8  %d  ", pat_8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-a9  %d  ", pat_9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-a...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

