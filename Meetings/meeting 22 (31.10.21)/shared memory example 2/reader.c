#include <sys/ipc.h>
#include <sys/shm.h>
#include <stdio.h>
  
int main()
{
    // ftok to generate unique key
    //key_t key = ftok("shmfile",65);
    key_t key=2000;
  
    // shmget returns an identifier in shmid
    int shmid = shmget(key,1000,0666);
  
    // shmat to attach to shared memory
    float *array=(float*) shmat(shmid,NULL,0);	
  
    float* pointer=array;
    printf("Data read from memory:\n");
    for(int i=0;i<5;i++,pointer++) printf("%f\n",*pointer);
      
    //detach from shared memory 
    shmdt(array);
    
    // destroy the shared memory
    shmctl(shmid,IPC_RMID,NULL);
     
    return 0;
}
