#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

//Grid:
#define NXP2 (362)
#define NYP2 (482)

//upping or lowering temp:
// #define up
#define movie

#define Tdefault (2.0)
#define Nbreak_up (6572)
#define stepsbetweenframes (25000000) //2.5 mill
#define mu (0.05)


#define NX (NXP2-2)
#define NY (NYP2-2)
#define ID(i,j) ((j) + NYP2*(i))
#define IDV(i,j) (*(lattice+ID(i,j)))


float J = 1.;
float beta = 1/Tdefault; //1/2.269; //crit = 2.2691853142130219681144
signed char *h; //array
int nthreads = 6; 

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
            IDV(i,j) = (rand()%2)*2-1;
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

        if (exp(-beta*2*IDV(i,j)*(-mu*((float) *(h + (i-1)*NY + (j-1)))/127.  + J*(IDV(i+1,j) + IDV(i-1,j) + IDV(i,j+1) + IDV(i,j-1))))*RAND_MAX >= rand())
        {
            IDV(i,j) = -IDV(i,j);
            if (i==1){
                IDV(NX+1,j) = IDV(i,j);
            }
            if (j==1)
            {
                IDV(i,NY+1) = IDV(i,j);
            }
            if (i==NX){
                IDV(0,j) = IDV(i,j);
            }
            if (j==NY)
            {
                IDV(i,0) = IDV(i,j);
            }
        }
    }
}

void step_lattice_ktimes(int *lattice, long int ksteps, int nthreadsnow, long int know){
    pthread_t tid;

    struct threadargstruc args[nthreadsnow];
    pthread_t tids[nthreadsnow];

    for (int i = 0; i < nthreadsnow; i++){
        // struct threadargstruc test;
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

void load_image_to_h(int id) {
    char filename[] = "sourceframes/frame0000.data";
    filename[21] = digits[(id/1)%10]; // 21
    filename[20] = digits[(id/10)%10];
    filename[19] = digits[(id/100)%10];
    filename[18] = digits[(id/1000)%10];

    FILE *fp;
    fp = fopen(filename, "rb");

    fread((void *) h, sizeof(signed char), NX*NY, fp);
    fclose(fp);
}


int inttochar(int val, int nplace)
{
    return digits[val];
}

void to_file(int *lattice, signed char *outarray, int i){
    
    for (int i = 1; i < NXP2-1; ++i)
    {
        for (int j = 1; j < NYP2-1; ++j)
        {
            *(outarray + (i-1) + (j-1)*NX) = IDV(i,j);
        }
    }

    int z;
    char filename[] = "images/arrayXYZW.data";
    z = 12;
    filename[z+0] = digits[(i/1000)%10];
    filename[z+1] = digits[(i/100)%10];
    filename[z+2] = digits[(i/10)%10];
    filename[z+3] = digits[(i/1)%10];

    FILE *fp;
    fp = fopen(filename, "wb");
    fwrite((void *) outarray, sizeof(signed char), NX*NY, fp);
    fclose(fp);
}

int main() {
    //INIT
    struct timespec start2, finish;
    double elapsed;

    h = malloc(sizeof(signed char)*NX*NY);

    signed char *outarray = malloc(sizeof(signed char)*NX*NY);
    int *lattice = malloc(sizeof(int)*NXP2*NYP2);
    init_lattice(lattice);

    clock_gettime(CLOCK_MONOTONIC, &start2);


    system("exec mkdir ./images/");
    system("exec rm -r ./images/*");

    // long long int kstepstotalnow = kstepstotal;

    float T;
    for (int i = 0; i < Nbreak_up; ++i)
    {
        load_image_to_h(i);
        printf("step %d/%d %f %d\n",i, Nbreak_up, 1./beta, NX*NY);
        step_lattice_ktimes(lattice, stepsbetweenframes, nthreads, i*stepsbetweenframes);
        to_file(lattice, outarray, i);
        // printlattice(lattice);
    }
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start2.tv_sec);
    elapsed += (finish.tv_nsec - start2.tv_nsec) / 1000000000.0;
    printf("elapsed=%f\n", elapsed);


    double freq = ((double) stepsbetweenframes*Nbreak_up/1000000.)/elapsed;
    printf("million steps per second = %f\n", freq);

    return 0;
}
