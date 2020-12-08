case 's':
        if(argv[argcount][2]=='P')
        sscanf(argv[argcount+1],"%d", &No_of_PARAMETERS);
        else if(argv[argcount][2]=='T')
        sscanf(argv[argcount+1],"%lf", &TOLERANCE);
        else if(argv[argcount][2]=='N')
        sscanf(argv[argcount+1],"%d", &MAX_No_of_ITERATIONS);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
"s...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'd':  //Discretization
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%d", &d0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &d10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &d11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%d", &d12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"%d", &d13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%d", &d14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%d", &d15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%d", &d16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%d", &d17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%d", &d18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%d", &d19);
          else
          sscanf(argv[argcount+1],"%d", &d1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &d20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &d21);
          else
          sscanf(argv[argcount+1],"%d", &d2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%d", &d3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%d", &d4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%d", &d5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%d", &d6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%d", &d7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%d", &d8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%d", &d9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
"d...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'I':  //Parameter Labels
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%d", &I0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &I10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &I11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%d", &I12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"%d", &I13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%d", &I14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%d", &I15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%d", &I16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%d", &I17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%d", &I18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%d", &I19);
          else
          sscanf(argv[argcount+1],"%d", &I1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%d", &I20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%d", &I21);
          else
          sscanf(argv[argcount+1],"%d", &I2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%d", &I3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%d", &I4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%d", &I5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%d", &I6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%d", &I7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%d", &I8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%d", &I9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
"I...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'm': //Minimum values
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%lf", &min_P0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%lf", &min_P10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%lf", &min_P11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%lf", &min_P12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"%lf", &min_P13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%lf", &min_P14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%lf", &min_P15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%lf", &min_P16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%lf", &min_P17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%lf", &min_P18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%lf", &min_P19);
          else
          sscanf(argv[argcount+1],"%lf", &min_P1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%lf", &min_P20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%lf", &min_P21);
          else
          sscanf(argv[argcount+1],"%lf", &min_P2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%lf", &min_P3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%lf", &min_P4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%lf", &min_P5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%lf", &min_P6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%lf", &min_P7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%lf", &min_P8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%lf", &min_P9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
"m...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;


case 'M':  //Maximum values
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%lf", &MAX_P0);
        else if(argv[argcount][2]=='1'){
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%lf", &MAX_P10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"M11 %lf", &MAX_P11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%lf", &MAX_P12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"%lf", &MAX_P13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%lf", &MAX_P14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%lf", &MAX_P15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%lf", &MAX_P16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%lf", &MAX_P17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%lf", &MAX_P18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%lf", &MAX_P19);
          else
          sscanf(argv[argcount+1],"%lf", &MAX_P1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%lf", &MAX_P20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%lf", &MAX_P21);
          else
          sscanf(argv[argcount+1],"%lf", &MAX_P2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%lf", &MAX_P3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%lf", &MAX_P4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%lf", &MAX_P5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%lf", &MAX_P6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%lf", &MAX_P7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%lf", &MAX_P8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%lf", &MAX_P9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
"M...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;


case 'A': // Accuracies
        if(argv[argcount][2]=='0')
        sscanf(argv[argcount+1],"%lf", &Acc_P0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%lf", &Acc_P10);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%lf", &Acc_P11);
          else if(argv[argcount][3]=='2')
          sscanf(argv[argcount+1],"%lf", &Acc_P12);
          else if(argv[argcount][3]=='3')
          sscanf(argv[argcount+1],"%lf", &Acc_P13);
          else if(argv[argcount][3]=='4')
          sscanf(argv[argcount+1],"%lf", &Acc_P14);
          else if(argv[argcount][3]=='5')
          sscanf(argv[argcount+1],"%lf", &Acc_P15);
          else if(argv[argcount][3]=='6')
          sscanf(argv[argcount+1],"%lf", &Acc_P16);
          else if(argv[argcount][3]=='7')
          sscanf(argv[argcount+1],"%lf", &Acc_P17);
          else if(argv[argcount][3]=='8')
          sscanf(argv[argcount+1],"%lf", &Acc_P18);
          else if(argv[argcount][3]=='9')
          sscanf(argv[argcount+1],"%lf", &Acc_P19);
          else
          sscanf(argv[argcount+1],"%lf", &Acc_P1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          sscanf(argv[argcount+1],"%lf", &Acc_P20);
          else if(argv[argcount][3]=='1')
          sscanf(argv[argcount+1],"%lf", &Acc_P21);
          else
          sscanf(argv[argcount+1],"%lf", &Acc_P2);
        }
        else if(argv[argcount][2]=='3')
        sscanf(argv[argcount+1],"%lf", &Acc_P3);
        else if(argv[argcount][2]=='4')
        sscanf(argv[argcount+1],"%lf", &Acc_P4);
        else if(argv[argcount][2]=='5')
        sscanf(argv[argcount+1],"%lf", &Acc_P5);
        else if(argv[argcount][2]=='6')
        sscanf(argv[argcount+1],"%lf", &Acc_P6);
        else if(argv[argcount][2]=='7')
        sscanf(argv[argcount+1],"%lf", &Acc_P7);
        else if(argv[argcount][2]=='8')
        sscanf(argv[argcount+1],"%lf", &Acc_P8);
        else if(argv[argcount][2]=='9')
        sscanf(argv[argcount+1],"%lf", &Acc_P9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
"A...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;
