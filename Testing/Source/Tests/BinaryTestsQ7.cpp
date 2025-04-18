#include "BinaryTestsQ7.h"
#include <stdio.h>
#include "Error.h"
#include "dsp/basic_math_functions.h"

#define SNR_THRESHOLD 20

/* 

Reference patterns are generated with
a double precision computation.

*/
#define ABS_ERROR_Q7 ((q7_t)5)

#define ABS_ERROR_Q63 ((q63_t)(1<<16))



/* Upper bound of maximum matrix dimension used by Python */
#define MAXMATRIXDIM 47

static void checkInnerTail(q7_t *b)
{
    ASSERT_TRUE(b[0] == 0);
    ASSERT_TRUE(b[1] == 0);
    ASSERT_TRUE(b[2] == 0);
    ASSERT_TRUE(b[3] == 0);
    ASSERT_TRUE(b[4] == 0);
    ASSERT_TRUE(b[5] == 0);
    ASSERT_TRUE(b[6] == 0);
    ASSERT_TRUE(b[7] == 0);
    ASSERT_TRUE(b[8] == 0);
    ASSERT_TRUE(b[9] == 0);
    ASSERT_TRUE(b[10] == 0);
    ASSERT_TRUE(b[11] == 0);
    ASSERT_TRUE(b[12] == 0);
    ASSERT_TRUE(b[13] == 0);
    ASSERT_TRUE(b[14] == 0);
    ASSERT_TRUE(b[15] == 0);

}

#define LOADDATA2()                         \
      const q7_t *inp1=input1.ptr();       \
      const q7_t *inp2=input2.ptr();       \
                                            \
      q7_t *ap=a.ptr();                    \
      q7_t *bp=b.ptr();                    \
                                            \
      q7_t *outp=output.ptr();             \
      q7_t *refp=ref.ptr();             \
      q7_t *tmpPtr=tmp.ptr();              \
      int16_t *dimsp = dims.ptr();          \
      int nbMatrixes = dims.nbSamples() / 4;\
      int rows,internal,columns,shift;            \
      int i;


#define PREPAREDATA2()                                                   \
      in1.numRows=rows;                                                  \
      in1.numCols=internal;                                               \
      memcpy((void*)ap,(const void*)inp1,2*sizeof(q7_t)*rows*internal);\
      in1.pData = ap;                                                    \
                                                                         \
      in2.numRows=internal;                                                  \
      in2.numCols=columns;                                               \
      memcpy((void*)bp,(const void*)inp2,2*sizeof(q7_t)*internal*columns);\
      in2.pData = bp;                                                    \
                                                                         \
      out.numRows=rows;                                                  \
      out.numCols=columns;                                               \
      out.pData = outp;


      

    void BinaryTestsQ7::test_mat_mult_q7()
    {     
      LOADDATA2();
      (void)shift;
      arm_status status;

      for(i=0;i < nbMatrixes ; i ++)
      {
          rows = *dimsp++;
          internal = *dimsp++;
          columns = *dimsp++;
          shift = *dimsp++;

          PREPAREDATA2();
          try
          {
              memset(tmpPtr,0,sizeof(q7_t)*internal*columns + 16);
              checkInnerTail(tmpPtr + internal*columns);
              
              shift = 0;
              if (shift!=0)
                arm_shift_q7(inp1,-shift,ap,rows*internal);
              status=arm_mat_mult_q7(&this->in1,&this->in2,&this->out,tmpPtr);
              if (shift!=0)
                arm_shift_q7(outp,shift,outp,rows*columns);
              ASSERT_TRUE(status==ARM_MATH_SUCCESS);
              ASSERT_NEAR_EQ_NB(outp,refp,ABS_ERROR_Q7,rows*columns);
    
              outp += (rows * columns);
              refp += (rows * columns);
              checkInnerTail(outp);
              checkInnerTail(tmpPtr + internal*columns);
          }
          catch(Client::Error &err)
          {
            char tmp[256];
            snprintf(tmp,256," (%d x %d x %d)\n",rows,internal,columns);
            strcat(err.details,tmp);
            throw(err);
          }
      }

      ASSERT_EMPTY_TAIL(output);

      ASSERT_SNR(output,ref,(q7_t)SNR_THRESHOLD);

      ASSERT_NEAR_EQ(output,ref,ABS_ERROR_Q7);

    } 





    void BinaryTestsQ7::setUp(Testing::testID_t id,std::vector<Testing::param_t>& params,Client::PatternMgr *mgr)
    {


      (void)params;
      switch(id)
      {
         case TEST_MAT_MULT_Q7_1:
            input1.reload(BinaryTestsQ7::INPUTS1_Q7_ID,mgr);
            input2.reload(BinaryTestsQ7::INPUTS2_Q7_ID,mgr);
            dims.reload(BinaryTestsQ7::DIMSBINARY1_S16_ID,mgr);

            ref.reload(BinaryTestsQ7::REFMUL1_Q7_ID,mgr);

            output.create(ref.nbSamples(),BinaryTestsQ7::OUT_Q7_ID,mgr);
            a.create(MAXMATRIXDIM*MAXMATRIXDIM,BinaryTestsQ7::TMPA_Q7_ID,mgr);
            b.create(MAXMATRIXDIM*MAXMATRIXDIM,BinaryTestsQ7::TMPB_Q7_ID,mgr);
            tmp.create(MAXMATRIXDIM*MAXMATRIXDIM,BinaryTestsQ7::TMP_Q7_ID,mgr);
         break;


    
      }
       

    
    }

    void BinaryTestsQ7::tearDown(Testing::testID_t id,Client::PatternMgr *mgr)
    {
       (void)id;
       output.dump(mgr);
    }
