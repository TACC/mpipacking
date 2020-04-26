/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%
   %%%% This program file is part of a set of performance tests
   %%%% by Victor Eijkhout, copyright 2018-2020
   %%%%
   %%%% packperf.c : performance of packing schemes
   %%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <mpi.h>

#define nschemes 9
enum schemes { CONTIGUOUS=0,COPYING=1,VECTOR=2,BUFFERED=3,PERSISTENT=4,
	       SUBARRAY=5,ONESIDED=6,PACKVECTOR=7,PACKELEMENT=8 };
  /*
   * 0: contiguous
   * 1: manual copying
   * 2: vector type
   * 3: buffered
   * 4: persistent
   * 5: subarray type
   * 6: one-sided
   * 7: packing by vector type
   * 8: packing by element WE SKIP THIS ONE
   */

#define NREPEATS 20
#ifndef L3SIZE
#define L3SIZE 55000000
#endif

float pure_rate;

void report_scheme(int scheme,int nrepeats,double *runtimes,double data_size) {
  double runtime = runtimes[scheme];

  // total data size divided by total runtime
  float drate = 1.*data_size/runtime;
  int // for printing only:
    mrate = 1.e-6*drate,
    grate = 1.e-9*drate;

  //  runtime /= nrepeats; // but now we want the time for a single transfer
  //runtimes[scheme] = runtime;

  if (scheme==CONTIGUOUS) pure_rate = drate;
  printf("Runtime=%6.1e (for %6.1e bytes), ",runtime,data_size);
  printf("effective data rate=%9.6fGb/s = %4.2f pure\n",
	 1.e-9*drate,drate/pure_rate);
  if (scheme>CONTIGUOUS) {
    int overhead_percentage = 100*(runtime-runtimes[0])/runtimes[0];
    printf("  overhead percentage over contiguous: %d\n",overhead_percentage);
  }
  if (scheme>COPYING) {
    int overhead_percentage = 100*(runtime-runtimes[0])/runtimes[0];
    overhead_percentage = 100*(runtime-runtimes[1])/runtimes[1];
    printf("  overhead percentage over copying:    %d\n",overhead_percentage);
  }

}

int test_one_scheme(int scheme, int nrepeats,
		     double *runtimes, int buffersize,double *rdata_size)
{

  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs,procno;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&procno);
  int process1 = 0, process2 = nprocs-1;
#ifdef MULTI
  if (nprocs%2==1) {printf("Even number of ranks for multi, please!!!\n"); return -1; }
  if (procno<nprocs/2) {
    process1 = procno; process2 = nprocs-procno-1;
  } else {
    process2 = procno; process1 = nprocs-procno-1;
  }
#endif  

  /*
   * Allocation
   */
  size_t buffer_bytes = buffersize*sizeof(double);
  // we need memory to copy from
  double *memory=NULL;
  posix_memalign((void**)&memory,64,2*buffer_bytes);
  if (memory==NULL) {
    printf("Could not allocate memory\n"); return 1; }
  memset(memory,37,2*buffer_bytes);
#ifdef DEBUG
  for (int i=0; i<2*buffersize; i++)
    memory[i] = (double)(i+1);
#endif
  double *buffer=NULL;
#ifndef MALLOC
  // we need a buffer to send. most of the time.
  posix_memalign((void**)&buffer,64,buffer_bytes);
  if (buffer==NULL) {
    printf("Could not allocate buffer\n"); return 1; }
  memset(buffer,37,buffer_bytes);
  //printf("buffer\n");
#endif

#ifdef FLUSH
  // we also allocate a large buffer to flush the L3 with
  char *l3buffer=NULL;
  posix_memalign((void**)&l3buffer,64,L3SIZE);
  if (l3buffer==NULL) {
    printf("Could not allocate l3buffer\n"); return 1; }
  memset(l3buffer,37,L3SIZE);
#endif

  // some cases have a special treatment of the buffer
  MPI_Win the_window;
  MPI_Request send_request;
  double *nicbuffer=NULL; int nicbuffer_size = buffer_bytes+MPI_BSEND_OVERHEAD;
  if (scheme==BUFFERED) {
    posix_memalign((void**)&nicbuffer,64,nicbuffer_size);
    if (nicbuffer==NULL) {
      printf("Could not allocate nicbuffer\n"); return 1; }
    MPI_Buffer_attach(nicbuffer,nicbuffer_size);
  } else if (scheme==ONESIDED) {
#ifdef MALLOC
    // the send buffer is inside the timing loop,
    // so we need a separate buffer for the window
    posix_memalign((void**)&buffer,64,buffer_bytes);
    if (buffer==NULL) {
      printf("Could not allocate buffer for window\n"); return 1; }
#endif
    MPI_Aint winsize=1;
    if (procno==process2)
      winsize = buffer_bytes;
    MPI_Win_create( buffer,winsize,sizeof(double),
		    MPI_INFO_NULL,comm,&the_window );
  }

  // receive buffer for pong message
  double zerobuf[1];

  if (procno==process1)
    printf("==== Scheme %d ====\n",scheme);

  double data_size=0.;
  float individual_times[nrepeats], individual_rates[nrepeats];

  /****************
   ****************
   **************** Repeated experiments to get reliable performance
   ****************
   ****************/

  double runtime,start_run = MPI_Wtime();
  int maxsize;
#ifdef MULTI
  MPI_Barrier(comm);
#endif

  /*
   * Define send type if needed
   * (this used to be inside the iterative loop.
   *  not sure why)
   */
  MPI_Datatype sendtype = MPI_DATATYPE_NULL;
  if (scheme==VECTOR || scheme==BUFFERED || scheme==PERSISTENT
      || scheme==ONESIDED || scheme==PACKVECTOR) {
    // vector type
    MPI_Type_vector(buffersize,1,2,MPI_DOUBLE,&sendtype);
    MPI_Type_commit(&sendtype);
  } else if (scheme==SUBARRAY) {
    // subarray type
    int sizes[2]={buffersize,2}, subsizes[2]={buffersize,1}, starts[2]={0,0};
    MPI_Type_create_subarray
      (2,sizes,subsizes,starts, MPI_ORDER_C,MPI_DOUBLE,&sendtype);
    MPI_Type_commit(&sendtype);
  }

  /*
   * Set up persistent sends
   */
  MPI_Request persistent_requests[2] = {MPI_REQUEST_NULL,MPI_REQUEST_NULL};
  if (scheme==PERSISTENT) {
    if (procno==process1) {
      // big buffer send
      MPI_Send_init( memory,1,sendtype,
		     process2,0,comm,persistent_requests+0);
      // ack receive
      MPI_Recv_init( zerobuf,0,MPI_DOUBLE,
		     process2,0,comm,persistent_requests+1);
    } else if (procno==process2) {
      // big buffer receive
      MPI_Recv_init( buffer,buffersize,MPI_DOUBLE,
		     process1,0,comm,persistent_requests+0);
      // ack send
      MPI_Send_init( zerobuf,0,MPI_DOUBLE,
		     process1,0,comm,persistent_requests+1);
    }
  }

  for (int repeat=0; repeat<nrepeats; repeat++) {

#ifdef MALLOC
    // send buffer inside the timing loop
    posix_memalign((void**)&buffer,64,buffer_bytes);
    if (buffer==NULL) {
      printf("Could not allocate buffer\n"); return 1; }
    //memset(buffer,37,buffer_bytes);
#endif

    if (scheme==PACKELEMENT)
      continue;

    // stuff for MPI_PACKED schemes
    int pack_position = 0, pack_buflen = buffer_bytes;

    // Send part
    if (procno==process1) {

      /*
       * Flush the L3
       */
#ifdef FLUSH
      for (int il3=0; il3<L3SIZE; il3++)
	l3buffer[il3] = l3buffer[il3] || (char)repeat;
#endif
      double this_runtime,this_starttime = MPI_Wtime();

      /*
       * First pack the data
       */
      if (scheme==CONTIGUOUS) {
	// contiguous buffer, no action required
	;
      } else if (scheme==COPYING) {
	// manual copying
	for (int i=0; i<buffersize; i++)
	  buffer[i] = memory[2*i];
      } else if ( scheme==VECTOR || scheme==BUFFERED || scheme==PERSISTENT
		  || scheme==SUBARRAY ) {
	// no buffering required for derived types
      } else if (scheme==ONESIDED) {
	// no buffering required for one-sided
      } else if (scheme==PACKELEMENT) {
#ifdef DEBUG
	{ int checksize;
	  MPI_Pack_size(buffersize,MPI_DOUBLE,comm,&checksize);
	  if (checksize>buffer_bytes)
	    printf("Packing by element overflows buffer: %d>%d\n",checksize,buffer_bytes);
	}
#endif
	for (int i=0; i<buffersize; i++) {
	  MPI_Pack(memory+2*i,1,MPI_DOUBLE,buffer,pack_buflen,&pack_position,comm);
	}
	maxsize = pack_position;
      } else if (scheme==PACKVECTOR) {
#ifdef DEBUG
	{ int checksize;
	  MPI_Pack_size(1,sendtype,comm,&checksize);
	  if (checksize>buffer_bytes)
	    printf("Packing by vector overflows buffer: %d>%d\n",checksize,buffer_bytes);
	}
#endif
	MPI_Pack(memory,1,sendtype,buffer,pack_buflen,&pack_position,comm);
	maxsize = pack_position;
      } else {
	printf("Scheme %d not yet implemented for packing\n",scheme);
      }

      /*
       * Now actually send the data
       */
      if (scheme<VECTOR)
	// send the buffer as such
	MPI_Send( buffer,buffersize,MPI_DOUBLE,process2,0,comm);
      else if (scheme==VECTOR || scheme==SUBARRAY)
	// send one derived type
	MPI_Send( memory,1,sendtype,process2,0,comm);
      else if (scheme==BUFFERED)
	// buffered send derived type
	MPI_Bsend( memory,1,sendtype,process2,0,comm);
      else if (scheme==PERSISTENT)
	MPI_Startall(2,persistent_requests);
      else if (scheme==ONESIDED) {
	// write derived datatype into window
	MPI_Aint zero_offset=0;
	MPI_Win_fence(0,the_window);
	MPI_Put( /* send: */ memory,1,sendtype,
		 /* where */ process2,zero_offset,
		 /* what */ buffersize,MPI_DOUBLE,
		 the_window);
	MPI_Win_fence(0,the_window);
      } else if (scheme==PACKELEMENT || scheme==PACKVECTOR) {
	// send as PACKED
	MPI_Send( buffer,pack_position,MPI_PACKED,process2,0,comm);
      } else {
	printf("Scheme %d not yet implemented for send\n",scheme); 
      }
      data_size = buffer_bytes;

      /*
       * Pong: receive the acknowledgement if you're not one-sided
       */
      if (scheme==PERSISTENT) {
	MPI_Waitall(2,persistent_requests,MPI_STATUSES_IGNORE);
      } else if (scheme!=ONESIDED) {
	// receive the pong
	MPI_Recv( zerobuf,0,MPI_DOUBLE,process2,0,comm,MPI_STATUS_IGNORE);
      }

      this_runtime = MPI_Wtime()-this_starttime;
      individual_times[repeat] = this_runtime;
    } else if (procno==process2) { // Receive part

      /*
       * Receiving process always gets a contiguous buffer
       * (except when packing)
       */
      MPI_Status status;
      if (scheme<ONESIDED && scheme!=PERSISTENT) {
	// usually receive contiguous buffer
	MPI_Recv( buffer,buffersize,MPI_DOUBLE,process1,0,comm,&status);
      } else if (scheme==PERSISTENT) {
	MPI_Startall(2,persistent_requests);
      } else if (scheme==ONESIDED) {
	// one-sided: no action required, just synchronization
	MPI_Win_fence(0,the_window);
	MPI_Win_fence(0,the_window);
      } else if (scheme==PACKELEMENT || scheme==PACKVECTOR) {
	// packing: receive as MPI_PACKED
	MPI_Recv( buffer,pack_buflen,MPI_PACKED, process1,0,comm,MPI_STATUS_IGNORE);
      } else {
	printf("Scheme %d not yet implemented for recv\n",scheme);
      }
#ifdef DEBUG
      if (repeat==0) {
	int recvcount;
	MPI_Get_count(&status,MPI_DOUBLE,&recvcount);
	if (recvcount!=buffersize)
	  printf("Recvcount %d /= buffersize %d\n",recvcount,buffersize);
	double checkval = 2 * ( recvcount * (recvcount+1.) /2 ) -1, checksum=0;
	for (int i=0; i<recvcount; i++)
	  checksum += buffer[i];
	if (fabs(checkval-checksum)/checkval>1.e-2)
	  printf("Received sum %e <> %e\n",checksum,checkval);
      }
#endif

      /*
       * Pong!
       */
      if (scheme==PERSISTENT) {
	MPI_Waitall(2,persistent_requests,MPI_STATUSES_IGNORE);
      } else if (scheme!=ONESIDED) {
	// in all but the one-sided scheme, acknowledge back
	MPI_Send( zerobuf,0,MPI_DOUBLE,process1,0,comm);
      }

#ifdef DEBUG
      if (scheme>0) {
	int errors=0;
	for (int i=0; i<buffersize; i++) {
	  double should = (double)(2*(i+1));
	  errors += fabs(buffer[i]-should)/should;
	}
	printf("Transmission errors in scheme %d: %d\n",scheme,errors);
      }
#endif

    }

#ifdef MALLOC
      free(buffer);
#endif

  } // end of repeat loop

  runtime = MPI_Wtime()-start_run;

  // free stuff
  if (scheme==VECTOR || scheme==SUBARRAY || scheme==ONESIDED || scheme==PACKVECTOR) {
    MPI_Type_free(&sendtype);
  }
      
  // report back
  if (procno==process1) {
    double avgtime = 0.;
    printf("Measured times:");
    for (int repeat=0; repeat<nrepeats; repeat++) {
      printf(" %7.3e",individual_times[repeat]);
      avgtime += individual_times[repeat];
    }
    avgtime /= nrepeats;
    printf(" : average=%7.3e\n",avgtime);
    if (nrepeats>1) {
      double deviate = 0.;
      for (int repeat=0; repeat<nrepeats; repeat++) {
	double t = avgtime - individual_times[repeat];
	deviate += t*t;
      }
      deviate = sqrt(deviate/(nrepeats-1));
      int nreject=0, accepted=-1;
      for (int repeat=0; repeat<nrepeats; repeat++) {
	double t = avgtime - individual_times[repeat];
	if (fabs(t)>deviate) {
	  /* printf("%2d: avg=%e tim=%e dev=%e, sig=%e\n", */
	  /* 	 repeat,avgtime,individual_times[repeat],fabs(t),deviate); */
	  individual_times[repeat] = -1.;
	  nreject++;
	} else
	  accepted = repeat;
      }
      avgtime = 0.; int nvalid=0; double mintime = individual_times[accepted];
      for (int repeat=0; repeat<nrepeats; repeat++)
	if (individual_times[repeat]>0) {
	  double thistime = individual_times[repeat];
	  avgtime += thistime;
	  if (thistime<mintime) mintime = thistime;
	  nvalid++;
	}
      avgtime /= nvalid;
      printf(".. rejecting %d runs;\n",nreject);
      printf(".. new average=%7.3e\n",avgtime);
      printf(".. minimum measurement: %7.3e\n",mintime);
    }
    runtimes[scheme] = avgtime;

    *rdata_size = data_size;
  }

  // free specific stuff
  if (scheme==BUFFERED) {
    MPI_Buffer_detach(nicbuffer,&nicbuffer_size);
    free(nicbuffer);
  } else if (scheme==ONESIDED)
    MPI_Win_free(&the_window);
  // free general stuff
  free(memory);
#ifndef MALLOC
  free(buffer);
#endif
#ifdef FLUSH
  free(l3buffer);
#endif
  return 0;
}

int main(int argc,char **argv) {

  MPI_Init(0,0);
  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs,procno;
  MPI_Comm_size(comm,&nprocs);
  MPI_Comm_rank(comm,&procno);
  int process1 = 0, process2 = nprocs-1;

#define Ndefault 100
  int N = Ndefault;
  if (argc==2) {
    if (!strcmp(argv[1],"-h")) {
      if (procno==0)
	printf("\nUsage: %s [number]\n\n",argv[0]);
      return 0;
    } else
      N = atoi(argv[1]);
  }

  int buffersize = N;
#ifdef DEBUG
  int nrepeats = 1;
#else
  int nrepeats = NREPEATS;
#endif
  if (procno==0) {
#ifdef FLUSH
    printf("\n\nUsing artificial L3 cache flush\n\n\n");
#endif
    printf("MPI timer resolution: %e\n\n\n",MPI_Wtick());
  }

  /*
   * Loop over all schemes
   */
  double runtimes[nschemes];
  for (int scheme=0; scheme<=PACKVECTOR; scheme++) {
    double data_size;
    int err =
      test_one_scheme(scheme,nrepeats,runtimes,buffersize,&data_size);
    if (err) return err;
    if (procno==0)
      report_scheme(scheme,nrepeats,runtimes,data_size);
  } // end scheme loop

  MPI_Finalize();
  return 0;
}
