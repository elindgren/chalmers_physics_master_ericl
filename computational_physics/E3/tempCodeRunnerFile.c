double d_kd_i = 0;
    for(int k=0; k<N; k++){
        if(k%100000 == 0){
            printf("k=%d \n", k);
        }
        d_kd_i = 0;
        for(int i=0; i<N-k; i++){
            d_kd_i += d[i+k]*d[i];
        }
        d_kd_i /= N;
        phi[k] = d_kd_i / avg_d_sq;
    }