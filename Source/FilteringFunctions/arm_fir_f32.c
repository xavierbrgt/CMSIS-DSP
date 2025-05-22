/* ----------------------------------------------------------------------
 * Project:      CMSIS DSP Library
 * Title:        arm_fir_f32.c
 * Description:  Floating-point FIR filter processing function
 *
 * $Date:        23 April 2021
 * $Revision:    V1.9.0
 *
 * Target Processor: Cortex-M and Cortex-A cores
 * -------------------------------------------------------------------- */
/*
 * Copyright (C) 2010-2021 ARM Limited or its affiliates. All rights reserved.
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the License); you may
 * not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "dsp/filtering_functions.h"
#include <stdio.h>
/**
  @ingroup groupFilters
 */

/**
  @defgroup FIR Finite Impulse Response (FIR) Filters

  This set of functions implements Finite Impulse Response (FIR) filters
  for Q7, Q15, Q31, and floating-point data types.  Fast versions of Q15 and Q31 are also provided.
  The functions operate on blocks of input and output data and each call to the function processes
  <code>blockSize</code> samples through the filter.  <code>pSrc</code> and
  <code>pDst</code> points to input and output arrays containing <code>blockSize</code> values.

  @par           Algorithm
                   The FIR filter algorithm is based upon a sequence of multiply-accumulate (MAC) operations.
                   Each filter coefficient <code>b[n]</code> is multiplied by a state variable which equals a previous input sample <code>x[n]</code>.
  <pre>
      y[n] = b[0] * x[n] + b[1] * x[n-1] + b[2] * x[n-2] + ...+ b[numTaps-1] * x[n-numTaps+1]
  </pre>
  @par
                   \image html FIR.GIF "Finite Impulse Response filter"
  @par
                   <code>pCoeffs</code> points to a coefficient array of size <code>numTaps</code>.
                   Coefficients are stored in time reversed order.
  @par
  <pre>
      {b[numTaps-1], b[numTaps-2], b[N-2], ..., b[1], b[0]}
  </pre>
  @par
                   <code>pState</code> points to a state array of size <code>numTaps + blockSize - 1</code>.
                   Samples in the state buffer are stored in the following order.
  @par
  <pre>
      {x[n-numTaps+1], x[n-numTaps], x[n-numTaps-1], x[n-numTaps-2]....x[n](==pSrc[0]), x[n+1](==pSrc[1]), ..., x[n+blockSize-1](==pSrc[blockSize-1])}
  </pre>

  @par
                   Note that the length of the state buffer exceeds the length of the coefficient array by <code>blockSize-1</code>.
                   The increased state buffer length allows circular addressing, which is traditionally used in the FIR filters,
                   to be avoided and yields a significant speed improvement.
                   The state variables are updated after each block of data is processed; the coefficients are untouched.

  @par           Instance Structure
                   The coefficients and state variables for a filter are stored together in an instance data structure.
                   A separate instance structure must be defined for each filter.
                   Coefficient arrays may be shared among several instances while state variable arrays cannot be shared.
                   There are separate instance structure declarations for each of the 4 supported data types.

  @par           Initialization Functions
                   There is also an associated initialization function for each data type.
                   The initialization function performs the following operations:
                   - Sets the values of the internal structure fields.
                   - Zeros out the values in the state buffer.
                   To do this manually without calling the init function, assign the follow subfields of the instance structure:
                   numTaps, pCoeffs, pState. Also set all of the values in pState to zero.

  @par
                   Use of the initialization function is optional.
                   However, if the initialization function is used, then the instance structure cannot be placed into a const data section.
                   To place an instance structure into a const data section, the instance structure must be manually initialized.
                   Set the values in the state buffer to zeros before static initialization.
                   The code below statically initializes each of the 4 different data type filter instance structures
  <pre>
      arm_fir_instance_f32 S = {numTaps, pState, pCoeffs};
      arm_fir_instance_q31 S = {numTaps, pState, pCoeffs};
      arm_fir_instance_q15 S = {numTaps, pState, pCoeffs};
      arm_fir_instance_q7 S =  {numTaps, pState, pCoeffs};
  </pre>
                   where <code>numTaps</code> is the number of filter coefficients in the filter; <code>pState</code> is the address of the state buffer;
                   <code>pCoeffs</code> is the address of the coefficient buffer.

  @par          Initialization of Helium version
                 For Helium version the array of coefficients must be padded with zero to contain
                 a full number of lanes.

                 The array length L must be a multiple of x. L = x * a :
                 - x is 4  for f32
                 - x is 4  for q31
                 - x is 4  for f16 (so managed like the f32 version and not like the q15 one)
                 - x is 8  for q15
                 - x is 16 for q7

                 The additional coefficients 
                 (x * a - numTaps) must be set to 0.
                 numTaps is still set to its right value in the init function. It means that
                 the implementation may require to read more coefficients due to the vectorization and
                 to avoid having to manage too many different cases in the code.

  @par          Helium state buffer
                 The state buffer must contain some additional temporary data
                 used during the computation but which is not the state of the FIR.
                 The first A samples are temporary data.
                 The remaining samples are the state of the FIR filter.

  @par
                 So the state buffer has size <code> numTaps + A + blockSize - 1 </code> :
                 - A is blockSize for f32
                 - A is 8*ceil(blockSize/8) for f16
                 - A is 8*ceil(blockSize/4) for q31
                 - A is 0 for other datatypes (q15 and q7)


  @par           Fixed-Point Behavior
                   Care must be taken when using the fixed-point versions of the FIR filter functions.
                   In particular, the overflow and saturation behavior of the accumulator used in each function must be considered.
                   Refer to the function specific documentation below for usage guidelines.

 */

/**
  @addtogroup FIR
  @{
 */

/**
  @brief         Processing function for floating-point FIR filter.
  @param[in]     S          points to an instance of the floating-point FIR filter structure
  @param[in]     pSrc       points to the block of input data
  @param[out]    pDst       points to the block of output data
  @param[in]     blockSize  number of samples to process
 */

#if defined(ARM_MATH_MVEF) && !defined(ARM_MATH_AUTOVECTORIZE)

#define FIR_F32_MAX_COEF_BLK        8

#define FIR_F32_CORE(pSamples, c, NB_TAPS)                                 \
        vecAcc0 = vdupq_n_f32(0.0f);                                       \
        for (int i = 0; i < NB_TAPS; i++) {                                \
            vecIn0 = vld1q(&pSamples[i]);                                  \
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c[i]);                        \
        }


#define NB_TAPS 4
__STATIC_INLINE void arm_fir_f32_1_4_mve(const arm_fir_instance_f32 * S, 
  const float32_t * __restrict pSrc, 
  float32_t * __restrict pDst, uint32_t blockSize)
{
    float32_t *pRefStatePtr = S->pState + blockSize;
    float32_t      *pState = pRefStatePtr; /* State pointer */
    const float32_t *pCoeffs = S->pCoeffs;      /* Coefficient pointer */
    float32_t      *pStateCur;  /* Points to the current sample of the state */
    const float32_t *pSamples;  /* Temporary pointer to the sample buffer */
    float32_t      *pOutput;    /* Temporary pointer to the output buffer */
    const float32_t *pTempSrc;  /* Temporary pointer to the source data */
    float32_t      *pTempDest;  /* Temporary pointer to the destination buffer */
    uint32_t        numTaps = S->numTaps;       /* Number of filter coefficients in the filter */
    int32_t         blkCnt;
    float32x4_t         vecIn0;
    float32x4_t         vecAcc0;
    float32_t       c[NB_TAPS];
    const float32_t *pCoeffsCur = pCoeffs;

    /*
     * pState points to state array which contains previous frame (numTaps - 1) samples
     * pStateCur points to the location where the new input data should be written
     */
    pStateCur = &(pState[(numTaps - 1u)]);
    pTempSrc = pSrc;

    pSamples = pState;
    pOutput = pDst;

    for (int i = 0; i < NB_TAPS; i++)
        c[i] = *pCoeffsCur++;

    blkCnt = blockSize >> 2;
    while (blkCnt > 0) {
        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        FIR_F32_CORE(pSamples, c, NB_TAPS);

        vst1q(pOutput, vecAcc0);

        pOutput += 4;
        pSamples += 4;

        blkCnt--;
    }

    blkCnt = blockSize & 3;
    if (blkCnt)
    {
        mve_pred16_t    p0 = vctp32q(blkCnt);

        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        FIR_F32_CORE(pSamples, c, NB_TAPS);

        vstrwq_p_f32(pOutput, vecAcc0, p0);
    }

    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pState[blockSize];
    pTempDest = pState;

    blkCnt = numTaps - 1;
    do {
        mve_pred16_t    p = vctp32q(blkCnt);

        vstrwq_p_f32(pTempDest, vldrwq_z_f32(pTempSrc, p), p);
        pTempSrc += 4;
        pTempDest += 4;
        blkCnt -= 4;
    }
    while (blkCnt > 0);
}
#undef NB_TAPS

__STATIC_INLINE void arm_fir_f32_5_8_mve(const arm_fir_instance_f32 * S, 
  const float32_t * __restrict pSrc, 
  float32_t * __restrict pDst, uint32_t blockSize)
{
    float32_t *pRefStatePtr = S->pState + blockSize;
    float32_t *pState = pRefStatePtr;      /* State pointer */
    const float32_t *pCoeffs = S->pCoeffs;    /* Coefficient pointer */
    const float32_t *pSamples;          /* Temporary pointer to the sample buffer */
    const float32_t *pTempSrc;          /* Temporary pointer to the source data */
    float32_t *pTempDest;               /* Temporary pointer to the destination buffer */
    uint32_t  numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    int32_t  blkCnt;
    float32_t c0, c1, c2, c3;
    float32_t c4, c5, c6, c7;


    pTempSrc = pSrc;
    pTempDest = &(pState[(numTaps - 1u)]);
    int cnt = blockSize;
    do {
        mve_pred16_t p0 = vctp32q(cnt);
        vstrwq_p_f32(pTempDest, vld1q(pTempSrc), p0);
        pTempDest += 4;
        pTempSrc += 4;
        cnt -= 4;
    } while(cnt > 0);



    pSamples = pState;
    c0 = *pCoeffs++;
    c1 = *pCoeffs++;
    c2 = *pCoeffs++;
    c3 = *pCoeffs++;
    c4 = *pCoeffs++;
    c5 = *pCoeffs++;
    c6 = *pCoeffs++;
    c7 = *pCoeffs++;

    cnt = blockSize >> 2;
    while(cnt > 0) 
    {
        float32x4_t vecAcc0;
        float32x4_t vecIn0;

        vecIn0 = vld1q(pSamples);
        vecAcc0 = vmulq(vecIn0, c0);
        vecIn0 = vld1q(&pSamples[1]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
        vecIn0 = vld1q(&pSamples[2]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
        vecIn0 = vld1q(&pSamples[3]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
        vecIn0 = vld1q(&pSamples[4]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
        vecIn0 = vld1q(&pSamples[5]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
        vecIn0 = vld1q(&pSamples[6]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
        vecIn0 = vld1q(&pSamples[7]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
        pSamples += 4;
        vst1q(pDst, vecAcc0);
        cnt--;
        pDst += 4;
    }

    cnt = blockSize & 3;
    if (cnt > 0) 
    {
        float32x4_t vecAcc0;
        float32x4_t vecIn0;

        mve_pred16_t p0 = vctp32q(cnt);

        vecIn0 = vld1q(pSamples);
        vecAcc0 = vmulq(vecIn0, c0);
        vecIn0 = vld1q(&pSamples[1]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
        vecIn0 = vld1q(&pSamples[2]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
        vecIn0 = vld1q(&pSamples[3]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
        vecIn0 = vld1q(&pSamples[4]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
        vecIn0 = vld1q(&pSamples[5]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
        vecIn0 = vld1q(&pSamples[6]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
        vecIn0 = vld1q(&pSamples[7]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
        vstrwq_p_f32(pDst, vecAcc0,p0);
    }


    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pState[blockSize];
    pTempDest = pState;
    blkCnt = numTaps;
    while (blkCnt > 0)
    {
        *pTempDest++ = *pTempSrc++;
        blkCnt--;
    }
}



ARM_DSP_ATTRIBUTE void arm_fir_f32(
const arm_fir_instance_f32 * S,
const float32_t * pSrc,
float32_t * pDst,
uint32_t blockSize)
{
    /* 
       S->pState is the arm_fir_partial_accu
       S->pState + blockSize is the FIR state
    */
    float32_t *pRefStatePtr = S->pState + blockSize;
    float32_t *pState = pRefStatePtr ;      /* State pointer */
    const float32_t *pCoeffs = S->pCoeffs;    /* Coefficient pointer */
    const float32_t *pSamples;          /* Temporary pointer to the sample buffer */
    float32_t *pOutput;                 /* Temporary pointer to the output buffer */
    const float32_t *pTempSrc;          /* Temporary pointer to the source data */
    float32_t *pTempDest;               /* Temporary pointer to the destination buffer */
    uint32_t  numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t  blkCnt;
    float32_t c0, c1, c2, c3;
    float32_t c4, c5, c6, c7;

    /*
     * [1 to 8 taps] specialized routines
     */
    if (numTaps <= 4)
    {
        arm_fir_f32_1_4_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 8)
    {
        arm_fir_f32_5_8_mve(S, pSrc, pDst, blockSize);
        return;
    }

    pTempSrc = pSrc;
    pTempDest = &(pState[(numTaps - 1u)]);
    int cnt = blockSize;
    do {
        mve_pred16_t p0 = vctp32q(cnt);
        vstrwq_p_f32(pTempDest, vld1q(pTempSrc), p0);
        pTempDest += 4;
        pTempSrc += 4;
        cnt -= 4;
    } while(cnt > 0);

    float32_t *partial_accu_ptr = S->pState;

    pSamples = pState;
    c0 = *pCoeffs++;
    c1 = *pCoeffs++;
    c2 = *pCoeffs++;
    c3 = *pCoeffs++;
    c4 = *pCoeffs++;
    c5 = *pCoeffs++;
    c6 = *pCoeffs++;
    c7 = *pCoeffs++;

    cnt = blockSize >> 2;
    while(cnt > 0) {
        float32x4_t vecAcc0;
        float32x4_t vecIn0;

        vecIn0 = vld1q(pSamples);
        vecAcc0 = vmulq(vecIn0, c0);
        vecIn0 = vld1q(&pSamples[1]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
        vecIn0 = vld1q(&pSamples[2]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
        vecIn0 = vld1q(&pSamples[3]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
        vecIn0 = vld1q(&pSamples[4]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
        vecIn0 = vld1q(&pSamples[5]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
        vecIn0 = vld1q(&pSamples[6]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
        vecIn0 = vld1q(&pSamples[7]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
        pSamples += 4;
        vst1q(partial_accu_ptr, vecAcc0);
        cnt--;
        partial_accu_ptr += 4;
    }

    cnt = blockSize & 3;
    if (cnt > 0) 
    {
        float32x4_t vecAcc0;
        float32x4_t vecIn0;

        mve_pred16_t p0 = vctp32q(cnt);

        vecIn0 = vld1q(pSamples);
        vecAcc0 = vmulq(vecIn0, c0);
        vecIn0 = vld1q(&pSamples[1]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
        vecIn0 = vld1q(&pSamples[2]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
        vecIn0 = vld1q(&pSamples[3]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
        vecIn0 = vld1q(&pSamples[4]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
        vecIn0 = vld1q(&pSamples[5]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
        vecIn0 = vld1q(&pSamples[6]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
        vecIn0 = vld1q(&pSamples[7]);
        vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
        vstrwq_p_f32(partial_accu_ptr, vecAcc0,p0);
    }

    int localTaps = numTaps - FIR_F32_MAX_COEF_BLK;
    int sample_offset = FIR_F32_MAX_COEF_BLK;
    while (localTaps > FIR_F32_MAX_COEF_BLK) {
        c0 = *pCoeffs++;
        c1 = *pCoeffs++;
        c2 = *pCoeffs++;
        c3 = *pCoeffs++;
        c4 = *pCoeffs++;
        c5 = *pCoeffs++;
        c6 = *pCoeffs++;
        c7 = *pCoeffs++;

        partial_accu_ptr = S->pState;
        pSamples = pState + sample_offset;
        int cnt = blockSize >> 2;
        while(cnt > 0) {
            float32x4_t vecAcc0;
            float32x4_t vecIn0;

            vecIn0 = vld1q(pSamples);
            vecAcc0 = vmulq(vecIn0, c0);
            vecIn0 = vld1q(&pSamples[1]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
            vecIn0 = vld1q(&pSamples[2]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
            vecIn0 = vld1q(&pSamples[3]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
            vecIn0 = vld1q(&pSamples[4]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
            vecIn0 = vld1q(&pSamples[5]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
            vecIn0 = vld1q(&pSamples[6]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
            vecIn0 = vld1q(&pSamples[7]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
            pSamples += 4;
            vecAcc0 += vld1q_f32(partial_accu_ptr);
            vst1q(partial_accu_ptr, vecAcc0);
            cnt--;
            partial_accu_ptr += 4;
        }

        cnt = blockSize & 3;
        if (cnt > 0) {
            float32x4_t vecAcc0;
            float32x4_t vecIn0;

            mve_pred16_t p0 = vctp32q(cnt);

            vecIn0 = vld1q(pSamples);
            vecAcc0 = vmulq(vecIn0, c0);
            vecIn0 = vld1q(&pSamples[1]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
            vecIn0 = vld1q(&pSamples[2]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
            vecIn0 = vld1q(&pSamples[3]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
            vecIn0 = vld1q(&pSamples[4]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
            vecIn0 = vld1q(&pSamples[5]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
            vecIn0 = vld1q(&pSamples[6]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
            vecIn0 = vld1q(&pSamples[7]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
            vecAcc0 += vld1q_f32(partial_accu_ptr);
            vstrwq_p_f32(partial_accu_ptr, vecAcc0,p0);
        }

        localTaps -= FIR_F32_MAX_COEF_BLK;
        sample_offset += FIR_F32_MAX_COEF_BLK;
    }

    pSamples = pState + sample_offset;

    if (localTaps > 4) {
        c0 = *pCoeffs++;
        c1 = *pCoeffs++;
        c2 = *pCoeffs++;
        c3 = *pCoeffs++;
        c4 = *pCoeffs++;
        c5 = *pCoeffs++;
        c6 = *pCoeffs++;
        c7 = *pCoeffs++;
        pOutput = pDst;

        partial_accu_ptr = S->pState;
        cnt = blockSize  >> 2;
        while(cnt > 0) {
            float32x4_t vecAcc0;
            float32x4_t vecIn0;

            vecIn0 = vld1q(pSamples);
            vecAcc0 = vmulq(vecIn0, c0);
            vecIn0 = vld1q(&pSamples[1]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
            vecIn0 = vld1q(&pSamples[2]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
            vecIn0 = vld1q(&pSamples[3]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
            vecIn0 = vld1q(&pSamples[4]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
            vecIn0 = vld1q(&pSamples[5]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
            vecIn0 = vld1q(&pSamples[6]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
            vecIn0 = vld1q(&pSamples[7]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
            pSamples += 4;
            float32x4_t pap = vld1q_f32(partial_accu_ptr);
            vst1q(pOutput, vecAcc0+pap);
            cnt--;
            partial_accu_ptr += 4;
            pOutput += 4;
        }

        cnt = blockSize  & 3;
        if (cnt > 0) {
            float32x4_t vecAcc0;
            float32x4_t vecIn0;

            mve_pred16_t p0 = vctp32q(cnt);

            vecIn0 = vld1q(pSamples);
            vecAcc0 = vmulq(vecIn0, c0);
            vecIn0 = vld1q(&pSamples[1]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
            vecIn0 = vld1q(&pSamples[2]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
            vecIn0 = vld1q(&pSamples[3]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
            vecIn0 = vld1q(&pSamples[4]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c4);
            vecIn0 = vld1q(&pSamples[5]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c5);
            vecIn0 = vld1q(&pSamples[6]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c6);
            vecIn0 = vld1q(&pSamples[7]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c7);
            float32x4_t pap = vld1q_f32(partial_accu_ptr);
            vstrwq_p_f32(pOutput, vecAcc0+pap,p0);
            pOutput += cnt;
        }
    }
    else {
        c0 = *pCoeffs++;
        c1 = *pCoeffs++;
        c2 = *pCoeffs++;
        c3 = *pCoeffs++;
        pOutput = pDst;

        partial_accu_ptr = S->pState;
        cnt = blockSize >> 2;
        while(cnt > 0) {
            float32x4_t vecAcc0;
            float32x4_t vecIn0;

            vecIn0 = vld1q(pSamples);
            vecAcc0 = vmulq(vecIn0, c0);
            vecIn0 = vld1q(&pSamples[1]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
            vecIn0 = vld1q(&pSamples[2]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
            vecIn0 = vld1q(&pSamples[3]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
            pSamples += 4;
            float32x4_t pap = vld1q_f32(partial_accu_ptr);
            vst1q(pOutput, vecAcc0+pap);
            cnt--;
            partial_accu_ptr += 4;
            pOutput += 4;
        }

        cnt = blockSize & 3;
        if (cnt > 0) {
            float32x4_t vecAcc0;
            float32x4_t vecIn0;

            mve_pred16_t p0 = vctp32q(cnt);

            vecIn0 = vld1q(pSamples);
            vecAcc0 = vmulq(vecIn0, c0);
            vecIn0 = vld1q(&pSamples[1]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c1);
            vecIn0 = vld1q(&pSamples[2]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c2);
            vecIn0 = vld1q(&pSamples[3]);
            vecAcc0 = vfmaq(vecAcc0, vecIn0, c3);
            float32x4_t pap = vld1q_f32(partial_accu_ptr);
            vstrwq_p_f32(pOutput, vecAcc0+pap,p0);
            pOutput += cnt;
        }
    }

    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pRefStatePtr[blockSize];
    pTempDest = pRefStatePtr;

    blkCnt = numTaps >> 2;
    while (blkCnt > 0)
    {
        vst1q(pTempDest, vld1q(pTempSrc));
        pTempSrc += 4;
        pTempDest += 4;
        blkCnt--;
    }
    blkCnt = numTaps & 3;
    if (blkCnt > 0)
    {
        mve_pred16_t p0 = vctp32q(blkCnt);
        vstrwq_p_f32(pTempDest, vld1q(pTempSrc), p0);
    }
}

#else
#if defined(ARM_MATH_NEON)

static void arm_update_fir_f32_1(float32x4_t *accv0, float32x4_t *accv1, float32x4_t x0, float32x4_t x1, float32x4_t b)
{
    float32x4_t xa = x0;
    float32x4_t xb = x1;
    *accv0 = vfmaq_n_f32(*accv0, xa, vgetq_lane_f32(b, 0));
    *accv1 = vfmaq_n_f32(*accv1, xb, vgetq_lane_f32(b, 0));
}
static void arm_update_fir_f32_2(float32x4_t *accv0, float32x4_t *accv1, float32x4_t *accv2, float32x4_t *accv3,
                                 float32x4_t x0, float32x4_t x1, float32x4_t x2, float32x4_t b)
{
    float32x4_t xa = x0;
    float32x4_t xb = x1;
    *accv0 = vfmaq_n_f32(*accv0, xa, vgetq_lane_f32(b, 0));
    *accv1 = vfmaq_n_f32(*accv1, xb, vgetq_lane_f32(b, 0));

    xa = vextq_f32(x0, x1, 1);
    xb = vextq_f32(x1, x2, 1);

    *accv2 = vfmaq_n_f32(*accv2, xa, vgetq_lane_f32(b, 1));
    *accv3 = vfmaq_n_f32(*accv3, xb, vgetq_lane_f32(b, 1));
}
static void arm_update_fir_f32_3(float32x4_t *accv0, float32x4_t *accv1, float32x4_t *accv2, float32x4_t *accv3,
                                 float32x4_t x0, float32x4_t x1, float32x4_t x2, float32x4_t b)
{
    float32x4_t xa = x0;
    float32x4_t xb = x1;
    *accv0 = vfmaq_n_f32(*accv0, xa, vgetq_lane_f32(b, 0));
    *accv1 = vfmaq_n_f32(*accv1, xb, vgetq_lane_f32(b, 0));

    xa = vextq_f32(x0, x1, 1);
    xb = vextq_f32(x1, x2, 1);

    *accv2 = vfmaq_n_f32(*accv2, xa, vgetq_lane_f32(b, 1));
    *accv3 = vfmaq_n_f32(*accv3, xb, vgetq_lane_f32(b, 1));

    xa = vextq_f32(x0, x1, 2);
    xb = vextq_f32(x1, x2, 2);

    *accv0 = vfmaq_n_f32(*accv0, xa, vgetq_lane_f32(b, 2));
    *accv1 = vfmaq_n_f32(*accv1, xb, vgetq_lane_f32(b, 2));
}

#define ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState2)                                                                  \
    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */                          \
    /* pStateCurnt points to the location where the new input data should be written */                                \
    pStateCurnt = &(S->pState[(numTaps - 1U)]);                                                                        \
                                                                                                                       \
    /* Loop unrolling */                                                                                               \
    blkCnt = blockSize >> 3;                                                                                           \
    /* Initialize coefficient pointer */                                                                               \
    int nTaps = S->numTaps >> 2;                                                                                       \
    px = pState2;                                                                                                      \
    while ((blkCnt > 0U))                                                                                              \
    {                                                                                                                  \
        /* Copy 8 samples at a time into state buffers */                                                              \
        samples0 = vld1q_f32(pSrc);                                                                                    \
        vst1q_f32(pStateCurnt, samples0);                                                                              \
                                                                                                                       \
        pStateCurnt += 4;                                                                                              \
        pSrc += 4;                                                                                                     \
                                                                                                                       \
        samples1 = vld1q_f32(pSrc);                                                                                    \
        vst1q_f32(pStateCurnt, samples1);                                                                              \
                                                                                                                       \
        pStateCurnt += 4;                                                                                              \
        pSrc += 4;                                                                                                     \
        pb = pCoeffs;                                                                                                  \
        /* Set the accumulators to zero */                                                                             \
        accv0 = vdupq_n_f32(0);                                                                                        \
        accv1 = vdupq_n_f32(0);                                                                                        \
        accv2 = vdupq_n_f32(0);                                                                                        \
        accv3 = vdupq_n_f32(0);                                                                                        \
                                                                                                                       \
    /* Initialize state pointer */

#define ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pState2)                                                             \
    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */                          \
    /* pStateCurnt points to the location where the new input data should be written */                                \
    /*pStateCurnt = &(S->pState[(numTaps - 1U)]);*/                                                                    \
                                                                                                                       \
    /* Loop unrolling */                                                                                               \
    blkCnt = blockSize >> 3;                                                                                           \
    /* Initialize coefficient pointer */                                                                               \
    px = pState2;                                                                                                      \
    while ((blkCnt > 0U))                                                                                              \
    {                                                                                                                  \
        pb = pCoeffs;                                                                                                  \
        /* Set the accumulators to zero */                                                                             \
        accv0 = vdupq_n_f32(0);                                                                                        \
        accv1 = vdupq_n_f32(0);                                                                                        \
        accv2 = vdupq_n_f32(0);                                                                                        \
        accv3 = vdupq_n_f32(0);                                                                                        \
                                                                                                                       \
    /* Initialize state pointer */

#define ARM_END_OF_FUNCTION_COMMON_FIR_F32()                                                                           \
    accv0 = vaddq_f32(accv0, accv2);                                                                                   \
    accv1 = vaddq_f32(accv1, accv3);                                                                                   \
    /* The result is stored in the destination buffer. */                                                              \
    vst1q_f32(pDst, accv0);                                                                                            \
    pDst += 4;                                                                                                         \
    vst1q_f32(pDst, accv1);                                                                                            \
    pDst += 4;                                                                                                         \
    /* Advance state pointer by 8 for the next 8 samples */                                                            \
    pState = pState + 8;                                                                                               \
    px += 8;                                                                                                           \
    blkCnt--;                                                                                                          \
    }                                                                                                                  \
                                                                                                                       \
    /* Tail */                                                                                                         \
    blkCnt = blockSize & 0x7;                                                                                          \
    while (blkCnt > 0U)                                                                                                \
    {                                                                                                                  \
        /* Copy one sample at a time into state buffer */                                                              \
        *pStateCurnt++ = *pSrc++;                                                                                      \
                                                                                                                       \
        /* Set the accumulator to zero */                                                                              \
        acc = 0.0f;                                                                                                    \
                                                                                                                       \
        /* Initialize state pointer */                                                                                 \
        px = pState;                                                                                                   \
                                                                                                                       \
        /* Initialize Coefficient pointer */                                                                           \
        pb = pCoeffs;                                                                                                  \
                                                                                                                       \
        i = numTaps;                                                                                                   \
                                                                                                                       \
        /* Perform the multiply-accumulates */                                                                         \
        do                                                                                                             \
        {                                                                                                              \
            /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3]    \
             * +...+ b[0] * x[0] */                                                                                    \
            acc += *px++ * *pb++;                                                                                      \
            i--;                                                                                                       \
                                                                                                                       \
        } while (i > 0U);                                                                                              \
                                                                                                                       \
        /* The result is stored in the destination buffer. */                                                          \
        *pDst++ = acc;                                                                                                 \
        /* Advance state pointer by 1 for the next sample */                                                           \
        pState = pState + 1;                                                                                           \
                                                                                                                       \
        blkCnt--;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    /* Processing is complete.                                                                                         \
    ** Now copy the last numTaps - 1 samples to the starting of the state buffer.                                      \
    ** This prepares the state buffer for the next function call. */                                                   \
                                                                                                                       \
    /* Points to the start of the state buffer */                                                                      \
    pStateCurnt = S->pState;                                                                                           \
                                                                                                                       \
    /* Copy numTaps number of values */                                                                                \
    tapCnt = numTaps - 1U;                                                                                             \
                                                                                                                       \
    /* Copy data */                                                                                                    \
    while (tapCnt > 0U)                                                                                                \
    {                                                                                                                  \
        *pStateCurnt++ = *pState++;                                                                                    \
                                                                                                                       \
        /* Decrement the loop counter */                                                                               \
        tapCnt--;                                                                                                      \
    }

#define ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP()                                                                      \
    /*end of inner loop 32*/                                                                                           \
    accv0 = vaddq_f32(accv0, accv2);                                                                                   \
    accv1 = vaddq_f32(accv1, accv3);                                                                                   \
    x1 = vld1q_f32(pState);                                                                                            \
    vst1q_f32(pDst, vaddq_f32(accv0, x1));                                                                             \
    pState += 4;                                                                                                       \
    pDst += 4;                                                                                                         \
    x2 = vld1q_f32(pState);                                                                                            \
    vst1q_f32(pDst, vaddq_f32(accv1, x2));                                                                             \
    pDst += 4;                                                                                                         \
    pState += 4;                                                                                                       \
    px += 8;                                                                                                           \
    /* Advance state pointer by 8 for the next 8 samples */                                                            \
    /*pState = pState + 8;*/                                                                                           \
                                                                                                                       \
    blkCnt--;                                                                                                          \
    }                                                                                                                  \
    /* Tail */                                                                                                         \
    blkCnt = blockSize & 0x7;                                                                                          \
    pStateacc = S->pState + 2 * blockSize - blkCnt;                                                                    \
    while (blkCnt > 0U)                                                                                                \
    {                                                                                                                  \
        /* Copy one sample at a time into state buffer */                                                              \
        /*printf("in %f suppose 7.374\n", pSrc);*/                                                                     \
        *pStateCurnt++ = *pSrc++;                                                                                      \
                                                                                                                       \
        /* Set the accumulator to zero */                                                                              \
        acc = 0.0f;                                                                                                    \
                                                                                                                       \
        /* Initialize state pointer */                                                                                 \
        px = pStateacc;                                                                                                \
                                                                                                                       \
        /* Initialize Coefficient pointer */                                                                           \
        pb = pCoeffs;                                                                                                  \
                                                                                                                       \
        i = numTaps;                                                                                                   \
                                                                                                                       \
        /* Perform the multiply-accumulates */                                                                         \
        do                                                                                                             \
        {                                                                                                              \
            /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3]    \
             * +...+ b[0] * x[0] */                                                                                    \
            /*printf("px %f pb %f ", px[0], pb[0]);*/                                                                  \
            acc += *px++ * *pb++;                                                                                      \
            /*printf("acc %f\n", acc);*/                                                                               \
            i--;                                                                                                       \
                                                                                                                       \
        } while (i > 0U);                                                                                              \
                                                                                                                       \
        /* The result is stored in the destination buffer. */                                                          \
        *pDst++ = acc;                                                                                                 \
        /* Advance state pointer by 1 for the next sample */                                                           \
        pStateacc = pStateacc + 1;                                                                                     \
                                                                                                                       \
        blkCnt--;                                                                                                      \
    }                                                                                                                  \
                                                                                                                       \
    /* Processing is complete.                                                                                         \
    ** Now copy the last numTaps - 1 samples to the starting of the state buffer.                                      \
    ** This prepares the state buffer for the next function call. */                                                   \
                                                                                                                       \
    /* Points to the start of the state buffer */                                                                      \
    pStateCurnt = S->pState + blockSize;                                                                               \
    pState = S->pState + 2 * blockSize;                                                                                \
                                                                                                                       \
    /* Copy numTaps number of values */                                                                                \
    tapCnt = numTaps - 1U;                                                                                             \
                                                                                                                       \
    /* Copy data */                                                                                                    \
    while (tapCnt > 0U)                                                                                                \
    {                                                                                                                  \
        *pStateCurnt++ = *pState++;                                                                                    \
                                                                                                                       \
        /* Decrement the loop counter */                                                                               \
        tapCnt--;                                                                                                      \
    }

static void arm_update_fir_f32(float32x4_t *accv0, float32x4_t *accv1, float32x4_t *accv2, float32x4_t *accv3,
                               float32x4_t x0, float32x4_t x1, float32x4_t x2, float32x4_t b)
{
    float32x4_t xa = x0;
    float32x4_t xb = x1;
    *accv0 = vfmaq_n_f32(*accv0, xa, vgetq_lane_f32(b, 0));
    *accv1 = vfmaq_n_f32(*accv1, xb, vgetq_lane_f32(b, 0));

    xa = vextq_f32(x0, x1, 1);
    xb = vextq_f32(x1, x2, 1);

    *accv2 = vfmaq_n_f32(*accv2, xa, vgetq_lane_f32(b, 1));
    *accv3 = vfmaq_n_f32(*accv3, xb, vgetq_lane_f32(b, 1));

    xa = vextq_f32(x0, x1, 2);
    xb = vextq_f32(x1, x2, 2);

    *accv0 = vfmaq_n_f32(*accv0, xa, vgetq_lane_f32(b, 2));
    *accv1 = vfmaq_n_f32(*accv1, xb, vgetq_lane_f32(b, 2));

    xa = vextq_f32(x0, x1, 3);
    xb = vextq_f32(x1, x2, 3);

    *accv2 = vfmaq_n_f32(*accv2, xa, vgetq_lane_f32(b, 3));
    *accv3 = vfmaq_n_f32(*accv3, xb, vgetq_lane_f32(b, 3));
}

#define ARM_INNER_LOOP_COEF_32_35_FIR_F32()                                                                            \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 16);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 20);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 24);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 28);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 32);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[6]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 36);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[7]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x2, x0, b[8]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x1 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[8]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[8]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_28_31_FIR_F32()                                                                            \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 16);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 20);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 24);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 28);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 32);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[6]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x0 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x1, x2, b[7]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x0 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[7]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x0 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[7]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_24_27_FIR_F32()                                                                            \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 16);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 20);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 24);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 28);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x2 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x0, x1, b[6]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x2 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[6]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x2 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[6]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_20_23_FIR_F32()                                                                            \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 16);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 20);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 24);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x2, x0, b[5]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x1 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_16_19_FIR_F32()                                                                            \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 16);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                              \
                                                                                                                       \
    x2 = vld1q_f32(px + 20);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x0 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x1, x2, b[4]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x0 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x0 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_12_15_FIR_F32()                                                                            \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    x1 = vld1q_f32(px + 16);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x2 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x0, x1, b[3]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x2 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x2 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_8_11_FIR_F32()                                                                             \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
                                                                                                                       \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
                                                                                                                       \
    x0 = vld1q_f32(px + 12);                                                                                           \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                              \
                                                                                                                       \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x2, x0, b[2]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x1 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_COEF_4_7_FIR_F32()                                                                              \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
    arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                              \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        x0 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x1, x2, b[1]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        x0 = vld1q_f32(px + (nTaps) * 4 + 8);                                                                          \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        x0 = vld1q_f32(px + nTaps * 4 + 8);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);                                        \
        break;                                                                                                         \
    default:                                                                                                           \
        break;                                                                                                         \
    }

#define ARM_INNER_LOOP_GENERIC_FIR_F32()                                                                               \
    x0 = vld1q_f32(px);                                                                                                \
    x1 = vld1q_f32(px + 4);                                                                                            \
    x2 = vld1q_f32(px + 8);                                                                                            \
    for (int i = 0; i < nTaps; i++)                                                                                    \
    {                                                                                                                  \
        b[0] = vld1q_f32(pb + i * 4);                                                                                  \
        arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                          \
        x0 = vld1q_f32(px + (i + 1) * 4);                                                                              \
        x1 = vld1q_f32(px + (i + 1) * 4 + 4);                                                                          \
        x2 = vld1q_f32(px + (i + 1) * 4 + 8);                                                                          \
    }                                                                                                                  \
    switch (nTail)                                                                                                     \
    {                                                                                                                  \
    case 0:                                                                                                            \
        break;                                                                                                         \
    case 1:                                                                                                            \
        b[0] = vld1q_f32(pb + (nTaps) * 4);                                                                            \
        arm_update_fir_f32_1(&accv0, &accv1, x0, x1, b[0]);                                                            \
        break;                                                                                                         \
    case 2:                                                                                                            \
        b[0] = vld1q_f32(pb + (nTaps) * 4);                                                                            \
        arm_update_fir_f32_2(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                        \
        break;                                                                                                         \
    case 3:                                                                                                            \
        b[0] = vld1q_f32(pb + (nTaps) * 4);                                                                            \
        arm_update_fir_f32_3(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);                                        \
        break;                                                                                                         \
                                                                                                                       \
    default:                                                                                                           \
        break;                                                                                                         \
    }

ARM_DSP_ATTRIBUTE void arm_fir_f32(const arm_fir_instance_f32 *S, const float32_t *pSrc, float32_t *pDst,
                                   uint32_t blockSize)
{
    float32_t *pState = S->pState; /* State pointer */
    float32_t *pStateacc = S->pState + blockSize;
    const float32_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    float32_t *pStateCurnt;                /* Points to the current sample of the state */
    float32_t *px;                         /* Temporary pointers for state buffer */
    const float32_t *pb;                   /* Temporary pointers for coefficient buffer */
    uint32_t numTaps = S->numTaps;         /* Number of filter coefficients in the filter */
    uint32_t i, tapCnt, blkCnt;            /* Loop counters */

    float32x4_t accv0, accv1, samples0, samples1;
    float32x4_t accv2, accv3;
    float32_t acc;
    float32x4_t x0, x1, x2;
    int32_t nbTaps = (numTaps) >> 2;
    int nTail = S->numTaps & 3;
    /* Specialized routine for small number of taps
        trying to prevent reload and to add enough delay between the reuse of the vectors */
    if (numTaps > 35)
    {
        nbTaps = (numTaps & 31) >> 2;
        int REMAINING_BLOCK_TAPS_TO_PROCESS = numTaps;
        int blocksdone = 0;

        pb = pCoeffs;
        float32x4_t b[9];
        for (int i = 0; i < 8; i++)
            b[i] = vld1q_f32(pb + i * 4);
        blocksdone += 8;

        /* S->pState points to state array which contains previous frame (numTaps - 1) samples */
        /* pStateCurnt points to the location where the new input data should be written */
        pStateCurnt = &(pStateacc[(numTaps - 1U)]);
        /* Initialize state pointer */
        px = pStateacc;
        /* Loop unrolling */
        blkCnt = blockSize >> 3;
        /* Initialize coefficient pointer */
        while ((blkCnt > 0U))
        {
            /* Copy 8 samples at a time into state buffers */
            samples0 = vld1q_f32(pSrc);
            vst1q_f32(pStateCurnt, samples0);

            pStateCurnt += 4;
            pSrc += 4;

            samples1 = vld1q_f32(pSrc);
            vst1q_f32(pStateCurnt, samples1);

            pStateCurnt += 4;
            pSrc += 4;
            pb = pCoeffs;
            /* Set the accumulators to zero */
            accv0 = vdupq_n_f32(0);
            accv1 = vdupq_n_f32(0);
            accv2 = vdupq_n_f32(0);
            accv3 = vdupq_n_f32(0);

            x0 = vld1q_f32(px);
            x1 = vld1q_f32(px + 4);
            x2 = vld1q_f32(px + 8);

            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);

            x0 = vld1q_f32(px + 12);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);

            x1 = vld1q_f32(px + 16);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);

            x2 = vld1q_f32(px + 20);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);

            x0 = vld1q_f32(px + 24);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);

            x1 = vld1q_f32(px + 28);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);

            x2 = vld1q_f32(px + 32);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[6]);

            x0 = vld1q_f32(px + 36);
            arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[7]);

            accv0 = vaddq_f32(accv0, accv2);
            accv1 = vaddq_f32(accv1, accv3);

            /* The result is stored in the destination buffer. */
            vst1q_f32(pState, accv0);
            pState += 4;
            vst1q_f32(pState, accv1);
            pState += 4;
            /* Advance state pointer by 8 for the next 8 samples */
            px += 8;
            blkCnt--;
        }
        REMAINING_BLOCK_TAPS_TO_PROCESS -= 32;
        while (REMAINING_BLOCK_TAPS_TO_PROCESS > 35)
        {
            for (int i = 0; i < 8; i++)
                b[i] = vld1q_f32(pb + i * 4 + blocksdone * 4);

            /* S->pState points to state array which contains previous frame (numTaps - 1) samples */
            /* pStateCurnt points to the location where the new input data should be written */
            pState = S->pState;
            /* Initialize state pointer */
            px = pStateacc + (32 * (blocksdone >> 3));
            /* Loop unrolling */
            blocksdone += 8;
            blkCnt = blockSize >> 3;
            while ((blkCnt > 0U))
            {
                /* Set the accumulators to zero */
                accv0 = vdupq_n_f32(0);
                accv1 = vdupq_n_f32(0);
                accv2 = vdupq_n_f32(0);
                accv3 = vdupq_n_f32(0);

                x0 = vld1q_f32(px);
                x1 = vld1q_f32(px + 4);
                x2 = vld1q_f32(px + 8);

                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[0]);

                x0 = vld1q_f32(px + 12);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[1]);

                x1 = vld1q_f32(px + 16);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[2]);

                x2 = vld1q_f32(px + 20);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[3]);

                x0 = vld1q_f32(px + 24);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[4]);

                x1 = vld1q_f32(px + 28);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x2, x0, x1, b[5]);

                x2 = vld1q_f32(px + 32);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x0, x1, x2, b[6]);

                x0 = vld1q_f32(px + 36);
                arm_update_fir_f32(&accv0, &accv1, &accv2, &accv3, x1, x2, x0, b[7]);

                accv0 = vaddq_f32(accv0, accv2);
                accv1 = vaddq_f32(accv1, accv3);
                /* The result is stored in the destination buffer. */
                x1 = vld1q_f32(pState);
                vst1q_f32(pState, vaddq_f32(accv0, x1));
                pState += 4;
                x2 = vld1q_f32(pState);
                vst1q_f32(pState, vaddq_f32(accv1, x2));
                pState += 4;

                /* Advance state pointer by 8 for the next 8 samples */
                px += 8;
                blkCnt--;
            }
            REMAINING_BLOCK_TAPS_TO_PROCESS -= 32;
        }
        if (REMAINING_BLOCK_TAPS_TO_PROCESS <= 35)
        {
            switch (nbTaps)
            {
            case 1: {
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                b[0] = vld1q_f32(pb + blocksdone * 4);
                if (nTail > 0)
                    b[1] = vld1q_f32(pb + blocksdone * 4 + 4);
                pState = S->pState;
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_4_7_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            case 2: {
                pState = S->pState;
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 2; i++)
                    b[i] = vld1q_f32(pb + blocksdone * 4 + i * 4);
                if (nTail > 0)
                    b[2] = vld1q_f32(pb + blocksdone * 4 + 8);
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_8_11_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            case 3: {
                pState = S->pState;
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 3; i++)
                    b[i] = vld1q_f32(pb + blocksdone * 4 + i * 4);
                if (nTail > 0)
                    b[3] = vld1q_f32(pb + blocksdone * 4 + 4 * 3);
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_12_15_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            case 4: {
                pState = S->pState;
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 4; i++)
                    b[i] = vld1q_f32(pb + blocksdone * 4 + i * 4);
                if (nTail > 0)
                    b[4] = vld1q_f32(pb + blocksdone * 4 + 4 * 4);
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_16_19_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            case 5: {
                pState = S->pState;
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 5; i++)
                    b[i] = vld1q_f32(pb + blocksdone * 4 + i * 4);
                if (nTail > 0)
                    b[5] = vld1q_f32(pb + blocksdone * 4 + 5 * 4);
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_20_23_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            case 6: {
                pState = S->pState;
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 6; i++)
                    b[i] = vld1q_f32(pb + blocksdone * 4 + i * 4);
                if (nTail > 0)
                    b[6] = vld1q_f32(pb + blocksdone * 4 + 6 * 4);
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_24_27_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            case 7: {
                pState = S->pState;
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 7; i++)
                    b[i] = vld1q_f32(pb + blocksdone * 4 + i * 4);
                if (nTail > 0)
                    b[7] = vld1q_f32(pb + blocksdone * 4 + 7 * 4);
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_28_31_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            default: {
                int nTaps = REMAINING_BLOCK_TAPS_TO_PROCESS >> 2;

                for (int i = 0; i < 8; i++)
                    b[i] = vld1q_f32(pb + i * 4 + blocksdone * 4);
                if (nTail > 0)
                    b[8] = vld1q_f32(pb + blocksdone * 4 + 8 * 4);
                pState = S->pState;
                ARM_START_OF_FUNCTION_COMMON_FIR_F32_LOOP(pStateacc + (32 * (blocksdone >> 3)));
                ARM_INNER_LOOP_COEF_32_35_FIR_F32();
                ARM_END_OF_FUNCTION_COMMON_FIR_F32_LOOP();
                return;
            }
            }
        }
    }
    else
    {
        switch (nbTaps)
        {
        case 1: {
            pb = pCoeffs;
            float32x4_t b[2];

            b[0] = vld1q_f32(pb);
            if (nTail > 0)
                b[1] = vld1q_f32(pb + 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_4_7_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 2: {
            pb = pCoeffs;
            float32x4_t b[3];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[2] = vld1q_f32(pb + nbTaps * 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_8_11_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 3: {
            pb = pCoeffs;
            float32x4_t b[4];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[3] = vld1q_f32(pb + nbTaps * 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_12_15_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 4: {
            pb = pCoeffs;
            float32x4_t b[5];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[4] = vld1q_f32(pb + nbTaps * 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_16_19_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 5: {
            pb = pCoeffs;
            float32x4_t b[6];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[5] = vld1q_f32(pb + nbTaps * 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_20_23_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 6: {
            pb = pCoeffs;
            float32x4_t b[7];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[6] = vld1q_f32(pb + nbTaps * 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_24_27_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 7: {
            pb = pCoeffs;
            float32x4_t b[8];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[7] = vld1q_f32(pb + nbTaps * 4);
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_28_31_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        case 8: {
            pb = pCoeffs;
            float32x4_t b[9];

            for (int i = 0; i < nbTaps; i++)
                b[i] = vld1q_f32(pb + i * 4);
            if (nTail > 0)
                b[8] = vld1q_f32(pb + nbTaps * 4);

            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_COEF_32_35_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        default: {
            float32x4_t b[1];
            ARM_START_OF_FUNCTION_COMMON_FIR_F32(pState);
            ARM_INNER_LOOP_GENERIC_FIR_F32();
            ARM_END_OF_FUNCTION_COMMON_FIR_F32();
            return;
        }
        }
    }
}

#else
ARM_DSP_ATTRIBUTE void arm_fir_f32(
  const arm_fir_instance_f32 * S,
  const float32_t * pSrc,
        float32_t * pDst,
        uint32_t blockSize)
{
        float32_t *pState = S->pState;                 /* State pointer */
  const float32_t *pCoeffs = S->pCoeffs;               /* Coefficient pointer */
        float32_t *pStateCurnt;                        /* Points to the current sample of the state */
        float32_t *px;                                 /* Temporary pointer for state buffer */
  const float32_t *pb;                                 /* Temporary pointer for coefficient buffer */
        float32_t acc0;                                /* Accumulator */
        uint32_t numTaps = S->numTaps;                 /* Number of filter coefficients in the filter */
        uint32_t i, tapCnt, blkCnt;                    /* Loop counters */

#if defined (ARM_MATH_LOOPUNROLL)
        float32_t acc1, acc2, acc3, acc4, acc5, acc6, acc7;     /* Accumulators */
        float32_t x0, x1, x2, x3, x4, x5, x6, x7;               /* Temporary variables to hold state values */
        float32_t c0;                                           /* Temporary variable to hold coefficient value */
#endif

  /* S->pState points to state array which contains previous frame (numTaps - 1) samples */
  /* pStateCurnt points to the location where the new input data should be written */
  pStateCurnt = &(S->pState[(numTaps - 1U)]);

#if defined (ARM_MATH_LOOPUNROLL)

  /* Loop unrolling: Compute 8 output values simultaneously.
   * The variables acc0 ... acc7 hold output values that are being computed:
   *
   *    acc0 =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0]
   *    acc1 =  b[numTaps-1] * x[n-numTaps]   + b[numTaps-2] * x[n-numTaps-1] + b[numTaps-3] * x[n-numTaps-2] +...+ b[0] * x[1]
   *    acc2 =  b[numTaps-1] * x[n-numTaps+1] + b[numTaps-2] * x[n-numTaps]   + b[numTaps-3] * x[n-numTaps-1] +...+ b[0] * x[2]
   *    acc3 =  b[numTaps-1] * x[n-numTaps+2] + b[numTaps-2] * x[n-numTaps+1] + b[numTaps-3] * x[n-numTaps]   +...+ b[0] * x[3]
   */

  blkCnt = blockSize >> 3U;

  while (blkCnt > 0U)
  {
    /* Copy 4 new input samples into the state buffer. */
    *pStateCurnt++ = *pSrc++;
    *pStateCurnt++ = *pSrc++;
    *pStateCurnt++ = *pSrc++;
    *pStateCurnt++ = *pSrc++;

    /* Set all accumulators to zero */
    acc0 = 0.0f;
    acc1 = 0.0f;
    acc2 = 0.0f;
    acc3 = 0.0f;
    acc4 = 0.0f;
    acc5 = 0.0f;
    acc6 = 0.0f;
    acc7 = 0.0f;

    /* Initialize state pointer */
    px = pState;

    /* Initialize coefficient pointer */
    pb = pCoeffs;

    /* This is separated from the others to avoid
     * a call to __aeabi_memmove which would be slower
     */
    *pStateCurnt++ = *pSrc++;
    *pStateCurnt++ = *pSrc++;
    *pStateCurnt++ = *pSrc++;
    *pStateCurnt++ = *pSrc++;

    /* Read the first 7 samples from the state buffer:  x[n-numTaps], x[n-numTaps-1], x[n-numTaps-2] */
    x0 = *px++;
    x1 = *px++;
    x2 = *px++;
    x3 = *px++;
    x4 = *px++;
    x5 = *px++;
    x6 = *px++;

    /* Loop unrolling: process 8 taps at a time. */
    tapCnt = numTaps >> 3U;

    while (tapCnt > 0U)
    {
      /* Read the b[numTaps-1] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-3] sample */
      x7 = *(px++);

      /* acc0 +=  b[numTaps-1] * x[n-numTaps] */
      acc0 += x0 * c0;

      /* acc1 +=  b[numTaps-1] * x[n-numTaps-1] */
      acc1 += x1 * c0;

      /* acc2 +=  b[numTaps-1] * x[n-numTaps-2] */
      acc2 += x2 * c0;

      /* acc3 +=  b[numTaps-1] * x[n-numTaps-3] */
      acc3 += x3 * c0;

      /* acc4 +=  b[numTaps-1] * x[n-numTaps-4] */
      acc4 += x4 * c0;

      /* acc1 +=  b[numTaps-1] * x[n-numTaps-5] */
      acc5 += x5 * c0;

      /* acc2 +=  b[numTaps-1] * x[n-numTaps-6] */
      acc6 += x6 * c0;

      /* acc3 +=  b[numTaps-1] * x[n-numTaps-7] */
      acc7 += x7 * c0;

      /* Read the b[numTaps-2] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-4] sample */
      x0 = *(px++);

      /* Perform the multiply-accumulate */
      acc0 += x1 * c0;
      acc1 += x2 * c0;
      acc2 += x3 * c0;
      acc3 += x4 * c0;
      acc4 += x5 * c0;
      acc5 += x6 * c0;
      acc6 += x7 * c0;
      acc7 += x0 * c0;

      /* Read the b[numTaps-3] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-5] sample */
      x1 = *(px++);

      /* Perform the multiply-accumulates */
      acc0 += x2 * c0;
      acc1 += x3 * c0;
      acc2 += x4 * c0;
      acc3 += x5 * c0;
      acc4 += x6 * c0;
      acc5 += x7 * c0;
      acc6 += x0 * c0;
      acc7 += x1 * c0;

      /* Read the b[numTaps-4] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-6] sample */
      x2 = *(px++);

      /* Perform the multiply-accumulates */
      acc0 += x3 * c0;
      acc1 += x4 * c0;
      acc2 += x5 * c0;
      acc3 += x6 * c0;
      acc4 += x7 * c0;
      acc5 += x0 * c0;
      acc6 += x1 * c0;
      acc7 += x2 * c0;

      /* Read the b[numTaps-4] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-6] sample */
      x3 = *(px++);
      /* Perform the multiply-accumulates */
      acc0 += x4 * c0;
      acc1 += x5 * c0;
      acc2 += x6 * c0;
      acc3 += x7 * c0;
      acc4 += x0 * c0;
      acc5 += x1 * c0;
      acc6 += x2 * c0;
      acc7 += x3 * c0;

      /* Read the b[numTaps-4] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-6] sample */
      x4 = *(px++);

      /* Perform the multiply-accumulates */
      acc0 += x5 * c0;
      acc1 += x6 * c0;
      acc2 += x7 * c0;
      acc3 += x0 * c0;
      acc4 += x1 * c0;
      acc5 += x2 * c0;
      acc6 += x3 * c0;
      acc7 += x4 * c0;

      /* Read the b[numTaps-4] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-6] sample */
      x5 = *(px++);

      /* Perform the multiply-accumulates */
      acc0 += x6 * c0;
      acc1 += x7 * c0;
      acc2 += x0 * c0;
      acc3 += x1 * c0;
      acc4 += x2 * c0;
      acc5 += x3 * c0;
      acc6 += x4 * c0;
      acc7 += x5 * c0;

      /* Read the b[numTaps-4] coefficient */
      c0 = *(pb++);

      /* Read x[n-numTaps-6] sample */
      x6 = *(px++);

      /* Perform the multiply-accumulates */
      acc0 += x7 * c0;
      acc1 += x0 * c0;
      acc2 += x1 * c0;
      acc3 += x2 * c0;
      acc4 += x3 * c0;
      acc5 += x4 * c0;
      acc6 += x5 * c0;
      acc7 += x6 * c0;

      /* Decrement loop counter */
      tapCnt--;
    }

    /* Loop unrolling: Compute remaining outputs */
    tapCnt = numTaps % 0x8U;

    while (tapCnt > 0U)
    {
      /* Read coefficients */
      c0 = *(pb++);

      /* Fetch 1 state variable */
      x7 = *(px++);

      /* Perform the multiply-accumulates */
      acc0 += x0 * c0;
      acc1 += x1 * c0;
      acc2 += x2 * c0;
      acc3 += x3 * c0;
      acc4 += x4 * c0;
      acc5 += x5 * c0;
      acc6 += x6 * c0;
      acc7 += x7 * c0;

      /* Reuse the present sample states for next sample */
      x0 = x1;
      x1 = x2;
      x2 = x3;
      x3 = x4;
      x4 = x5;
      x5 = x6;
      x6 = x7;

      /* Decrement loop counter */
      tapCnt--;
    }

    /* Advance the state pointer by 8 to process the next group of 8 samples */
    pState = pState + 8;

    /* The results in the 8 accumulators, store in the destination buffer. */
    *pDst++ = acc0;
    *pDst++ = acc1;
    *pDst++ = acc2;
    *pDst++ = acc3;
    *pDst++ = acc4;
    *pDst++ = acc5;
    *pDst++ = acc6;
    *pDst++ = acc7;


    /* Decrement loop counter */
    blkCnt--;
  }

  /* Loop unrolling: Compute remaining output samples */
  blkCnt = blockSize % 0x8U;

#else

  /* Initialize blkCnt with number of taps */
  blkCnt = blockSize;

#endif /* #if defined (ARM_MATH_LOOPUNROLL) */

  while (blkCnt > 0U)
  {
    /* Copy one sample at a time into state buffer */
    *pStateCurnt++ = *pSrc++;

    /* Set the accumulator to zero */
    acc0 = 0.0f;

    /* Initialize state pointer */
    px = pState;

    /* Initialize Coefficient pointer */
    pb = pCoeffs;

    i = numTaps;

    /* Perform the multiply-accumulates */
    while (i > 0U)
    {
      /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */
      acc0 += *px++ * *pb++;

      i--;
    }

    /* Store result in destination buffer. */
    *pDst++ = acc0;

    /* Advance state pointer by 1 for the next sample */
    pState = pState + 1U;

    /* Decrement loop counter */
    blkCnt--;
  }

  /* Processing is complete.
     Now copy the last numTaps - 1 samples to the start of the state buffer.
     This prepares the state buffer for the next function call. */

  /* Points to the start of the state buffer */
  pStateCurnt = S->pState;

#if defined (ARM_MATH_LOOPUNROLL)

  /* Loop unrolling: Compute 4 taps at a time */
  tapCnt = (numTaps - 1U) >> 2U;

  /* Copy data */
  while (tapCnt > 0U)
  {
    *pStateCurnt++ = *pState++;
    *pStateCurnt++ = *pState++;
    *pStateCurnt++ = *pState++;
    *pStateCurnt++ = *pState++;

    /* Decrement loop counter */
    tapCnt--;
  }

  /* Calculate remaining number of copies */
  tapCnt = (numTaps - 1U) % 0x4U;

#else

  /* Initialize tapCnt with number of taps */
  tapCnt = (numTaps - 1U);

#endif /* #if defined (ARM_MATH_LOOPUNROLL) */

  /* Copy remaining data */
  while (tapCnt > 0U)
  {
    *pStateCurnt++ = *pState++;

    /* Decrement loop counter */
    tapCnt--;
  }

}

#endif /* #if defined(ARM_MATH_NEON) */
#endif /* defined(ARM_MATH_MVEF) && !defined(ARM_MATH_AUTOVECTORIZE) */

/**
* @} end of FIR group
*/
