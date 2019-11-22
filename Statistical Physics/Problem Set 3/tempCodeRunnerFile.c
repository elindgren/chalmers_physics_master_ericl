    prev_m = new_m;  // save old value
            printf("Prev: %.2f \n", prev_m);
            new_m = magnetization(N, lat);
            printf("New: %.2f \n", new_m);
            m[steps] = new_m;
            printf("magnetization: %.2f \n", new_m);