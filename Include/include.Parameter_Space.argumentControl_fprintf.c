case 's':
        if(argv[argcount][2]=='P')
        fprintf(fp,"-sP  %d  ", Table->S->No_of_PARAMETERS);
        else if(argv[argcount][2]=='T')
        fprintf(fp,"-sT  %g  ", Table->S->TOLERANCE);
        else if(argv[argcount][2]=='N')
        fprintf(fp,"-sN  %d  ", Table->S->MAX_No_of_ITERATIONS);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-s...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'd':  //Discretization
        if(argv[argcount][2]=='0')
        fprintf(fp,"-d0  %d  ", d0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-d10  %d  ", d10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-d11  %d  ", d11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-d12  %d  ", d12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-d13  %d  ", d13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-d14  %d  ", d14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-d15  %d  ", d15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-d16  %d  ", d16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-d17  %d  ", d17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-d18  %d  ", d18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-d19  %d  ", d19);
          else
          fprintf(fp,"-d1  %d  ", d1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-d20  %d  ", d20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-d21  %d  ", d21);
          else
          fprintf(fp,"-d2  %d  ", d2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-d3  %d  ", d3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-d4  %d  ", d4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-d5  %d  ", d5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-d6  %d  ", d6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-d7  %d  ", d7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-d8  %d  ", d8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-d9  %d  ", d9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-d...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'I':  //Parameter Labels
        if(argv[argcount][2]=='0')
        fprintf(fp,"-I0  %d  ", I0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-I10  %d  ", I10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-I11  %d  ", I11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-I12  %d  ", I12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-I13  %d  ", I13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-I14  %d  ", I14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-I15  %d  ", I15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-I16  %d  ", I16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-I17  %d  ", I17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-I18  %d  ", I18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-I19  %d  ", I19);
          else
          fprintf(fp,"-I1  %d  ", I1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-I20  %d  ", I20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-I21  %d  ", I21);
          else
          fprintf(fp,"-I2  %d  ", I2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-I3  %d  ", I3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-I4  %d  ", I4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-I5  %d  ", I5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-I6  %d  ", I6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-I7  %d  ", I7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-I8  %d  ", I8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-I9  %d  ", I9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-I...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;

case 'm': //Minimum values
        if(argv[argcount][2]=='0')
        fprintf(fp,"-m0  %g  ", min_P0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-m10  %g  ", min_P10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-m11  %g  ", min_P11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-m12  %g  ", min_P12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-m13  %g  ", min_P13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-m14  %g  ", min_P14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-m15  %g  ", min_P15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-m16  %g  ", min_P16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-m17  %g  ", min_P17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-m18  %g  ", min_P18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-m19  %g  ", min_P19);
          else
          fprintf(fp,"-m1  %g  ", min_P1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-m20  %g  ", min_P20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-m21  %g  ", min_P21);
          else
          fprintf(fp,"-m2  %g  ", min_P2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-m3  %g  ", min_P3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-m4  %g  ", min_P4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-m5  %g  ", min_P5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-m6  %g  ", min_P6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-m7  %g  ", min_P7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-m8  %g  ", min_P8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-m9  %g  ", min_P9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-m...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;


case 'M':  //Maximum values
        if(argv[argcount][2]=='0')
        fprintf(fp,"-M0  %g  ", MAX_P0);
        else if(argv[argcount][2]=='1'){
          if(argv[argcount][3]=='0')
          fprintf(fp,"-M10  %g  ", MAX_P10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-M11 %g  ", MAX_P11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-M12  %g  ", MAX_P12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-M13  %g  ", MAX_P13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-M14  %g  ", MAX_P14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-M15  %g  ", MAX_P15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-M16  %g  ", MAX_P16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-M17  %g  ", MAX_P17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-M18  %g  ", MAX_P18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-M19  %g  ", MAX_P19);
          else
          fprintf(fp,"-M1  %g  ", MAX_P1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-M20  %g  ", MAX_P20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-M21  %g  ", MAX_P21);
          else
          fprintf(fp,"-M2  %g  ", MAX_P2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-M3  %g  ", MAX_P3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-M4  %g  ", MAX_P4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-M5  %g  ", MAX_P5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-M6  %g  ", MAX_P6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-M7  %g  ", MAX_P7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-M8  %g  ", MAX_P8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-M9  %g  ", MAX_P9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-M...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;


case 'A': // Accuracies
        if(argv[argcount][2]=='0')
        fprintf(fp,"-A0  %g  ", Acc_P0);
        else if(argv[argcount][2]=='1') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-A10  %g  ", Acc_P10);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-A11  %g  ", Acc_P11);
          else if(argv[argcount][3]=='2')
          fprintf(fp,"-A12  %g  ", Acc_P12);
          else if(argv[argcount][3]=='3')
          fprintf(fp,"-A13  %g  ", Acc_P13);
          else if(argv[argcount][3]=='4')
          fprintf(fp,"-A14  %g  ", Acc_P14);
          else if(argv[argcount][3]=='5')
          fprintf(fp,"-A15  %g  ", Acc_P15);
          else if(argv[argcount][3]=='6')
          fprintf(fp,"-A16  %g  ", Acc_P16);
          else if(argv[argcount][3]=='7')
          fprintf(fp,"-A17  %g  ", Acc_P17);
          else if(argv[argcount][3]=='8')
          fprintf(fp,"-A18  %g  ", Acc_P18);
          else if(argv[argcount][3]=='9')
          fprintf(fp,"-A19  %g  ", Acc_P19);
          else
          fprintf(fp,"-A1  %g  ", Acc_P1);
        }
        else if(argv[argcount][2]=='2') {
          if(argv[argcount][3]=='0')
          fprintf(fp,"-A20  %g  ", Acc_P20);
          else if(argv[argcount][3]=='1')
          fprintf(fp,"-A21  %g  ", Acc_P21);
          else
          fprintf(fp,"-A2  %g  ", Acc_P2);
        }
        else if(argv[argcount][2]=='3')
        fprintf(fp,"-A3  %g  ", Acc_P3);
        else if(argv[argcount][2]=='4')
        fprintf(fp,"-A4  %g  ", Acc_P4);
        else if(argv[argcount][2]=='5')
        fprintf(fp,"-A5  %g  ", Acc_P5);
        else if(argv[argcount][2]=='6')
        fprintf(fp,"-A6  %g  ", Acc_P6);
        else if(argv[argcount][2]=='7')
        fprintf(fp,"-A7  %g  ", Acc_P7);
        else if(argv[argcount][2]=='8')
        fprintf(fp,"-A8  %g  ", Acc_P8);
        else if(argv[argcount][2]=='9')
        fprintf(fp,"-A9  %g  ", Acc_P9);
        else {
          printf("**invalid command line argument: -  %s  %c\n",
          "-A...", argv[argcount][2]);
          exit(0);
        }
        skip++;
        break;
