#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/vlimit.h>

int g09run_()
{
        char* comm = "g09 geom.com";
        vlimit(LIM_CORE,0);
        return system(comm);
}

int g16run_()
{
        char* comm = "g16 geom.com";
        vlimit(LIM_CORE,0);
        return system(comm);
}


int copytmps_()
{
        char* comm = "cp tmp.chk tmps.chk";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int copytmpsr_()
{
        char* comm = "cp tmps.chk tmp.chk";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int copytmps2_()
{
        char* comm = "cp tmps.chk tmps2.chk";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int copytmps2r_()
{
        char* comm = "cp tmps2.chk tmp.chk";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int copytmpb_()
{
        char* comm = "cp tmp.chk tmpb.chk";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int copytmpbr_()
{
        char* comm = "cp tmpb.chk tmp.chk";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int copygeomirc_()
{
        char* comm = "cp geom.log irc_g09.log";
        vlimit(LIM_CORE,0);
        return system(comm);
}
int trimtemp_()
{
  char filename[20];
  char buffer[200];
  FILE *fp;

  sprintf(filename,"temp.log");
  fp=fopen(filename,"r");
  fscanf(fp,"%s ",buffer);
  fclose(fp);
  //  printf(" I have read %s",buffer);
  sprintf(filename,"temp1.log");
  fp=fopen(filename,"w");
  fprintf(fp,"%s=",buffer);
  fclose(fp);
  return 0;
}



void molprorun_(int* nproc)
{
  char comm[70];
  sprintf(comm," molprop -n %d molpro.inp",*nproc);
  printf("command is %s \n",comm);
  vlimit(LIM_CORE,0);
  system(comm);
  //  sleep(120);
  return;
  //return system(comm);
}



int commrun_(char *command, int ll)
{
  command[ll--]='\0';
  printf("command is %s \n",command);
  vlimit(LIM_CORE,0);
  return system(command);
}


