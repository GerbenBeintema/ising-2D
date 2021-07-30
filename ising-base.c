#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

//Grid:
#define NXP2 (1002)
#define NYP2 (1002)

//upping or lowering temp:
#define lowT (2.0)
#define highT (2.4)
#define up

#define Tdefault (2.265)
#define Nbreak_up (2000)
#define kstepstotal (7200000000000)


#define NX (NXP2-2)
#define NY (NYP2-2)
#define ID(i,j) ((j) + NYP2*(i))
#define IDV(i,j) (*(lattice+ID(i,j)))

float J = 1.;
float beta = 1/Tdefault; //1/2.269; //crit = 2.2691853142130219681144
float mu = 0.;
float h = 0;
int nthreads = 6; //1 = 47750.509319, 239591.429527
long int transitiontable[2][2][2][2][2];
const char digits[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";


struct threadargstruc
{
    int *lattice;
    long int kmax;
    pthread_t tid;
    unsigned int seed;
};

void printlattice(int *lattice){
    for (int i = 1; i < NXP2-1; ++i)
    {
        for (int j = 1; j < NYP2-1; ++j)
        {
            if (IDV(i,j)==1){
                printf("#");
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }
}

void init_lattice(int *lattice){
    for (int i = 0; i < NXP2; ++i)
    {
        for (int j = 0; j < NYP2; ++j)
        {
#ifdef up
            IDV(i,j) = 1;
#else
            IDV(i,j) = (rand()%2)*2-1;
#endif
        }
    }
}

float cal_mean(int *lattice){
    long int S = 0;
    for (int i = 1; i < NXP2-1; ++i)
    {
        for (int j = 1; j < NYP2-1; ++j)
        {
            S = S + IDV(i,j);
        }
    }
    return ((float) S) / (NX*NY);
}

void *step_lattice_ktimes_threaded(void *argvoid)
{
    struct threadargstruc *arg = (struct threadargstruc *) argvoid;
    int *lattice = arg->lattice;
    long int kmax = arg->kmax;
    pthread_t tid = arg->tid;
    srand(arg->seed);

    int i, j;
    for (long int k = 0; k < kmax; ++k)
    {
        i = (rand()%NX)+1;
        j = (rand()%NY)+1;
        if (transitiontable[IDV(i,j)>0][IDV(i+1,j)>0][IDV(i-1,j)>0][IDV(i,j+1)>0][IDV(i,j-1)>0] >= rand())
        {
            IDV(i,j) = -IDV(i,j);

            //update boundaries for periodic boundary conditions
            if (i==1){
                IDV(NX+1,j) = IDV(i,j);
            }
            else if (i==NX){
                IDV(0,j) = IDV(i,j);
            }
            if (j==1)
            {
                IDV(i,NY+1) = IDV(i,j);
            }
            else if (j==NY)
            {
                IDV(i,0) = IDV(i,j);
            }
        }
    }
}

void create_transitiontable() {
    float p;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                for (int l = 0; l < 2; ++l)
                {
                    for (int m = 0; m < 2; ++m)
                    {
                        p = exp(-beta*2*(i*2-1)*(-mu*h  + J*((j*2-1) + (k*2-1) + (l*2-1) + (m*2-1))));
                        if (p>=1)
                        {
                            transitiontable[i][j][k][l][m] = RAND_MAX;
                        }
                        else {
                            transitiontable[i][j][k][l][m] = p*RAND_MAX;
                        }
                    }
                }
            }
        }
    }
}


void step_lattice_ktimes(int *lattice, long int ksteps, int nthreadsnow, long int know){
    pthread_t tid;

    struct threadargstruc args[nthreadsnow];
    pthread_t tids[nthreadsnow];

    for (int i = 0; i < nthreadsnow; i++){
        args[i].lattice = lattice;
        args[i].kmax = ksteps/nthreadsnow;
        args[i].tid = tid;
        args[i].seed = i + (know % 12345678);
        (void) pthread_create(&tid, NULL, step_lattice_ktimes_threaded, (void *)&args[i]);
        tids[i] = tid;
    }

    for (int i = 0; i < nthreadsnow; ++i)
    {
        pthread_join(tids[i], NULL);
    }
}

int inttochar(int val, int nplace)
{
    return digits[val];
}

void to_file(int *lattice, signed char *outarray, int id){
    
    for (int i = 1; i < NXP2-1; ++i)
    {
        for (int j = 1; j < NYP2-1; ++j)
        {
            *(outarray + (i-1) + (j-1)*NX) = IDV(i,j);
        }
    }

    int z;
#ifndef down
#ifndef up
    char filename[] = "images/arrayXYZW.data";
    z = 12;
#endif
#endif  
#ifdef down
    char filename[] = "imagesdown/arrayXYZW.data";
    z = 12+4;
#endif
#ifdef up
    char filename[] = "imagesup/arrayXYZW.data";
    z = 12+2;
#endif
    filename[z+3] = digits[(id/1)%10];
    filename[z+2] = digits[(id/10)%10];
    filename[z+1] = digits[(id/100)%10];
    filename[z+0] = digits[(id/1000)%10];

    FILE *fp;
    fp = fopen(filename, "wb");
    fwrite((void *) outarray, sizeof(signed char), NX*NY, fp);
    fclose(fp);
}

int main() {
    //INIT
    struct timespec start2, finish;
    double elapsed;

    signed char *outarray = malloc(sizeof(signed char)*NX*NY);
    int *lattice = malloc(sizeof(int)*NXP2*NYP2);
    init_lattice(lattice);
    create_transitiontable();
    clock_gettime(CLOCK_MONOTONIC, &start2);

#ifdef up
    system("exec mkdir ./imagesup/");
    system("exec rm -r ./imagesup/*");
#endif
#ifdef down
    system("exec mkdir ./imagesdown/");
    system("exec rm -r ./imagesdown/*");
#endif
#ifndef down
#ifndef up
    system("exec mkdir ./images/");
    system("exec rm -r ./images/*");
#endif
#endif
    // long long int kstepstotalnow = kstepstotal;


    //Main Loop
    float T;
    for (int i = 0; i < Nbreak_up; ++i)
    {
#ifdef down
        T = highT - i/(Nbreak_up - 1.)*(highT-lowT);
        beta = 1./T;
        printf("down T=%f\n", T);
        create_transitiontable();
#endif
#ifdef  up 
        T = lowT + i/(Nbreak_up - 1.)*(highT-lowT);
        beta = 1./T;
        printf("up T=%f\n", T);
        create_transitiontable();
#endif

        printf("step %d/%d, T=%f NX*NY=%d\n", i, Nbreak_up, 1./beta, NX*NY);
        step_lattice_ktimes(lattice, kstepstotal/Nbreak_up, nthreads, i*kstepstotal/Nbreak_up);
        to_file(lattice, outarray, i);

    }

    //Finish Up
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start2.tv_sec);
    elapsed += (finish.tv_nsec - start2.tv_nsec) / 1000000000.0;
    printf("elapsed=%f\n", elapsed);


    double freq = ((double) kstepstotal/1000000.)/elapsed;
    printf("million steps per second = %f\n", freq);

    return 0;
}
