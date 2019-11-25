Save results
    filename[0] = '\0';  // Empty filename string
    strcat(filename, dir);
    strcat(filename, "magnetization.dat");
    f = fopen(filename, "w");
    for (int i = 0; i < Ts; i++)
    {
        fprintf(f, "%.16f \t %.16f \n", T[i], m[i]);
    }
    fclose(f);