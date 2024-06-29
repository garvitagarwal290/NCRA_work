//#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdio.h>
#include <errno.h>
  
int main()
{
    // ftok to generate unique key
    //key_t key = ftok("garvit",65);
    key_t key= 2000;	
    
    // shmget returns an identifier in shmid
    int shmid = shmget(key,1000,0666|IPC_CREAT);
  
    // shmat to attach to shared memory
    float *array=(float*) shmat(shmid,NULL,0);
    
    if(array==(float*) -1) printf("%d \n",errno);
  
    float* pointer=array;
    printf("Write Data :\n");
    for(int i=0;i<5;i++,pointer++) scanf("%f",pointer);
  
    //detach from shared memory 
    shmdt(array);
    //shmctl(shmid,IPC_RMID,NULL);
  
    return 0;
}
