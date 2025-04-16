/* ----------------------------------------------------------------------
 * Project:      CMSIS DSP Library
 * Title:        arm_fir_q31.c
 * Description:  Q31 FIR filter processing function
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

/**
  @ingroup groupFilters
 */

/**
  @addtogroup FIR
  @{
 */

/**
  @brief         Processing function for Q31 FIR filter.
  @param[in]     S          points to an instance of the Q31 FIR filter structure
  @param[in]     pSrc       points to the block of input data
  @param[out]    pDst       points to the block of output data
  @param[in]     blockSize  number of samples to process

  @par           Scaling and Overflow Behavior
                   The function is implemented using an internal 64-bit accumulator.
                   The accumulator has a 2.62 format and maintains full precision of the intermediate multiplication results but provides only a single guard bit.
                   Thus, if the accumulator result overflows it wraps around rather than clip.
                   In order to avoid overflows completely the input signal must be scaled down by log2(numTaps) bits.
                   After all multiply-accumulates are performed, the 2.62 accumulator is right shifted by 31 bits and saturated to 1.31 format to yield the final result.

 @remark
                   Refer to \ref arm_fir_fast_q31() for a faster but less precise implementation of this filter.
 */
#if defined(ARM_MATH_MVEI) && !defined(ARM_MATH_AUTOVECTORIZE)

#include "arm_helium_utils.h"

#define FIR_Q31_CORE(nbAcc, nbVecTaps, pSample, vecCoeffs)       \
    for (int j = 0; j < nbAcc; j++)                              \
    {                                                            \
        const q31_t *pSmp = &pSamples[j];                        \
        q31x4_t vecIn0;                                          \
        q63_t acc[4];                                            \
                                                                 \
        acc[j] = 0;                                              \
        for (int i = 0; i < nbVecTaps; i++)                      \
        {                                                        \
            vecIn0 = vld1q(pSmp + 4 * i);                        \
            acc[j] = vrmlaldavhaq(acc[j], vecIn0, vecCoeffs[i]); \
        }                                                        \
        *pOutput++ = (q31_t)asrl(acc[j], 23);                    \
    }

#define FIR_Q31_CORE_STR_PARTIAL(nbAcc, nbVecTaps, pSample, vecCoeffs) \
    for (int j = 0; j < nbAcc; j++)                                    \
    {                                                                  \
        const q31_t *pSmp = &pSamples[j];                              \
        q31x4_t vecIn0;                                                \
                                                                       \
        acc[j] = 0;                                                    \
        for (int i = 0; i < nbVecTaps; i++)                            \
        {                                                              \
            vecIn0 = vld1q(pSmp + 4 * i);                              \
            acc[j] = vrmlaldavhaq(acc[j], vecIn0, vecCoeffs[i]);       \
        }                                                              \
        *arm_fir_partial_accu_ptr++ = acc[j];                          \
    }

#define FIR_Q31_CORE_LD_PARTIAL(nbAcc, nbVecTaps, pSample, vecCoeffs) \
    for (int j = 0; j < nbAcc; j++)                                   \
    {                                                                 \
        const q31_t *pSmp = &pSamples[j];                             \
        q31x4_t vecIn0;                                               \
                                                                      \
        acc[j] = *arm_fir_partial_accu_ptr++;                         \
                                                                      \
        for (int i = 0; i < nbVecTaps; i++)                           \
        {                                                             \
            vecIn0 = vld1q(pSmp + 4 * i);                             \
            acc[j] = vrmlaldavhaq(acc[j], vecIn0, vecCoeffs[i]);      \
        }                                                             \
        *pOutput++ = (q31_t)asrl(acc[j], 23);                         \
    }

#define FIR_Q31_MAIN_CORE()                                                                  \
    {                                                                                        \
        q31_t *pRefStatePtr = S->pState + 2 * ARM_ROUND_UP(blockSize, 4);                    \
        q31_t *pState = pRefStatePtr;      /* State pointer */                               \
        const q31_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */                         \
        q31_t *pStateCur;                  /* Points to the current sample of the state */   \
        const q31_t *pSamples;             /* Temporary pointer to the sample buffer */      \
        q31_t *pOutput;                    /* Temporary pointer to the output buffer */      \
        const q31_t *pTempSrc;             /* Temporary pointer to the source data */        \
        q31_t *pTempDest;                  /* Temporary pointer to the destination buffer */ \
        uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */ \
        int32_t blkCnt;                                                                      \
                                                                                             \
        /*                                                                                   \
         * load coefs                                                                        \
         */                                                                                  \
        q31x4_t vecCoeffs[NBVECTAPS];                                                        \
                                                                                             \
        for (int i = 0; i < NBVECTAPS; i++)                                                  \
            vecCoeffs[i] = vld1q(pCoeffs + 4 * i);                                           \
                                                                                             \
        /*                                                                                   \
         * pState points to state array which contains previous frame (numTaps - 1) samples  \
         * pStateCur points to the location where the new input data should be written       \
         */                                                                                  \
        pStateCur = &(pState[(numTaps - 1u)]);                                               \
        pTempSrc = pSrc;                                                                     \
        pSamples = pState;                                                                   \
        pOutput = pDst;                                                                      \
                                                                                             \
        blkCnt = blockSize >> 2;                                                             \
        while (blkCnt > 0)                                                                   \
        {                                                                                    \
            /*                                                                               \
             * Save 4 input samples in the history buffer                                    \
             */                                                                              \
            vstrwq_s32(pStateCur, vldrwq_s32(pTempSrc));                                     \
            pStateCur += 4;                                                                  \
            pTempSrc += 4;                                                                   \
                                                                                             \
            FIR_Q31_CORE(4, NBVECTAPS, pSamples, vecCoeffs);                                 \
                                                                                             \
            pSamples += 4;                                                                   \
            /*                                                                               \
             * Decrement the sample block loop counter                                       \
             */                                                                              \
            blkCnt--;                                                                        \
        }                                                                                    \
                                                                                             \
        /* tail */                                                                           \
        int32_t residual = blockSize & 3;                                                    \
        switch (residual)                                                                    \
        {                                                                                    \
        case 3:                                                                              \
        {                                                                                    \
            for (int i = 0; i < residual; i++)                                               \
                *pStateCur++ = *pTempSrc++;                                                  \
                                                                                             \
            FIR_Q31_CORE(3, NBVECTAPS, pSamples, vecCoeffs);                                 \
        }                                                                                    \
        break;                                                                               \
                                                                                             \
        case 2:                                                                              \
        {                                                                                    \
            for (int i = 0; i < residual; i++)                                               \
                *pStateCur++ = *pTempSrc++;                                                  \
                                                                                             \
            FIR_Q31_CORE(2, NBVECTAPS, pSamples, vecCoeffs);                                 \
        }                                                                                    \
        break;                                                                               \
                                                                                             \
        case 1:                                                                              \
        {                                                                                    \
            for (int i = 0; i < residual; i++)                                               \
                *pStateCur++ = *pTempSrc++;                                                  \
                                                                                             \
            FIR_Q31_CORE(1, NBVECTAPS, pSamples, vecCoeffs);                                 \
        }                                                                                    \
        break;                                                                               \
        }                                                                                    \
                                                                                             \
        /*                                                                                   \
         * Copy the samples back into the history buffer start                               \
         */                                                                                  \
        pTempSrc = &pState[blockSize];                                                       \
        pTempDest = pState;                                                                  \
                                                                                             \
        blkCnt = (numTaps - 1) >> 2;                                                         \
        while (blkCnt > 0)                                                                   \
        {                                                                                    \
            vstrwq_s32(pTempDest, vldrwq_s32(pTempSrc));                                     \
            pTempSrc += 4;                                                                   \
            pTempDest += 4;                                                                  \
            blkCnt--;                                                                        \
        }                                                                                    \
        blkCnt = (numTaps - 1) & 3;                                                          \
        if (blkCnt > 0)                                                                      \
        {                                                                                    \
            mve_pred16_t p0 = vctp32q(blkCnt);                                               \
            vstrwq_p_s32(pTempDest, vldrwq_z_s32(pTempSrc, p0), p0);                         \
        }                                                                                    \
    }

static void arm_fir_q31_1_4_mve(const arm_fir_instance_q31 *S,
                                const q31_t *__restrict pSrc,
                                q31_t *__restrict pDst, uint32_t blockSize)
{
    q31_t *pRefStatePtr = S->pState + 2 * ARM_ROUND_UP(blockSize, 4);
    q31_t *pState = pRefStatePtr;      /* State pointer */
    const q31_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q31_t *pStateCur;                  /* Points to the current sample of the state */
    const q31_t *pSamples;             /* Temporary pointer to the sample buffer */
    q31_t *pOutput;                    /* Temporary pointer to the output buffer */
    const q31_t *pTempSrc;             /* Temporary pointer to the source data */
    q31_t *pTempDest;                  /* Temporary pointer to the destination buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t blkCnt;
    q31x4_t vecIn0;

    /*
     * pState points to state array which contains previous frame (numTaps - 1) samples
     * pStateCur points to the location where the new input data should be written
     */
    pStateCur = &(pState[(numTaps - 1u)]);
    pTempSrc = pSrc;
    pSamples = pState;
    pOutput = pDst;

    q63_t acc0 = 0, acc1 = 0, acc2 = 0, acc3 = 0;
    /*
     * load 4 coefs
     */
    q31x4_t vecCoeffs = *(q31x4_t *)pCoeffs;

    blkCnt = blockSize >> 2;
    while (blkCnt > 0U)
    {
        const q31_t *pSamplesTmp = pSamples;

        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        vecIn0 = vld1q(pSamplesTmp);
        acc0 = vrmlaldavhq(vecIn0, vecCoeffs);

        vecIn0 = vld1q(&pSamplesTmp[1]);
        acc1 = vrmlaldavhq(vecIn0, vecCoeffs);

        vecIn0 = vld1q(&pSamplesTmp[2]);
        acc2 = vrmlaldavhq(vecIn0, vecCoeffs);

        vecIn0 = vld1q(&pSamplesTmp[3]);
        acc3 = vrmlaldavhq(vecIn0, vecCoeffs);

        acc0 = asrl(acc0, 23);
        acc1 = asrl(acc1, 23);
        acc2 = asrl(acc2, 23);
        acc3 = asrl(acc3, 23);

        *pOutput++ = (q31_t)acc0;
        *pOutput++ = (q31_t)acc1;
        *pOutput++ = (q31_t)acc2;
        *pOutput++ = (q31_t)acc3;

        pSamples += 4;
        /*
         * Decrement the sample block loop counter
         */
        blkCnt--;
    }

    uint32_t residual = blockSize & 3;
    switch (residual)
    {
    case 3:
    {
        /*
         * Save 4 input samples in the history buffer
         */
        *(q31x4_t *)pStateCur = *(q31x4_t *)pTempSrc;
        pStateCur += 4;
        pTempSrc += 4;

        vecIn0 = vld1q(pSamples);
        acc0 = vrmlaldavhq(vecIn0, vecCoeffs);

        vecIn0 = vld1q(&pSamples[1]);
        acc1 = vrmlaldavhq(vecIn0, vecCoeffs);

        vecIn0 = vld1q(&pSamples[2]);
        acc2 = vrmlaldavhq(vecIn0, vecCoeffs);

        acc0 = asrl(acc0, 23);
        acc1 = asrl(acc1, 23);
        acc2 = asrl(acc2, 23);

        *pOutput++ = (q31_t)acc0;
        *pOutput++ = (q31_t)acc1;
        *pOutput++ = (q31_t)acc2;
    }
    break;

    case 2:
    {
        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        vecIn0 = vld1q(pSamples);
        acc0 = vrmlaldavhq(vecIn0, vecCoeffs);

        vecIn0 = vld1q(&pSamples[1]);
        acc1 = vrmlaldavhq(vecIn0, vecCoeffs);

        acc0 = asrl(acc0, 23);
        acc1 = asrl(acc1, 23);

        *pOutput++ = (q31_t)acc0;
        *pOutput++ = (q31_t)acc1;
    }
    break;

    case 1:
    {
        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        vecIn0 = vld1q(pSamples);
        acc0 = vrmlaldavhq(vecIn0, vecCoeffs);

        acc0 = asrl(acc0, 23);

        *pOutput++ = (q31_t)acc0;
    }
    break;
    }

    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pState[blockSize];
    pTempDest = pState;

    blkCnt = (numTaps - 1) >> 2;
    while (blkCnt > 0U)
    {
        vst1q(pTempDest, vld1q(pTempSrc));
        pTempSrc += 4;
        pTempDest += 4;
        blkCnt--;
    }
    blkCnt = (numTaps - 1) & 3;
    if (blkCnt > 0U)
    {
        mve_pred16_t p0 = vctp32q(blkCnt);
        vstrwq_p_s32(pTempDest, vld1q(pTempSrc), p0);
    }
}

static void arm_fir_q31_5_8_mve(const arm_fir_instance_q31 *S,
                                const q31_t *__restrict pSrc,
                                q31_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 8
#define NBVECTAPS (NBTAPS / 4)
    FIR_Q31_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q31_9_12_mve(const arm_fir_instance_q31 *S,
                                 const q31_t *__restrict pSrc,
                                 q31_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 12
#define NBVECTAPS (NBTAPS / 4)
    FIR_Q31_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q31_13_16_mve(const arm_fir_instance_q31 *S,
                                  const q31_t *__restrict pSrc,
                                  q31_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 16
#define NBVECTAPS (NBTAPS / 4)
    FIR_Q31_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q31_17_20_mve(const arm_fir_instance_q31 *S,
                                  const q31_t *__restrict pSrc,
                                  q31_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 20
#define NBVECTAPS (NBTAPS / 4)
    FIR_Q31_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q31_21_24_mve(const arm_fir_instance_q31 *S,
                                  const q31_t *__restrict pSrc,
                                  q31_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 24
#define NBVECTAPS (NBTAPS / 4)
    FIR_Q31_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q31_25_28_mve(const arm_fir_instance_q31 *S,
                                  const q31_t *__restrict pSrc,
                                  q31_t *__restrict pDst, uint32_t blockSize)
{
#define NBTAPS 28
#define NBVECTAPS (NBTAPS / 4)
    FIR_Q31_MAIN_CORE();
#undef NBVECTAPS
#undef NBTAPS
}

static void arm_fir_q31_29_32_mve(const arm_fir_instance_q31 *S,
                                  const q31_t *__restrict pSrc,
                                  q31_t *__restrict pDst,
                                  uint32_t blockSize)
{
    q31_t *pRefStatePtr = S->pState + 2 * ARM_ROUND_UP(blockSize, 4);
    q31_t *pState = pRefStatePtr;      /* State pointer */
    const q31_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q31_t *pStateCur;                  /* Points to the current sample of the state */
    const q31_t *pSamples;             /* Temporary pointer to the sample buffer */
    q31_t *pOutput;                    /* Temporary pointer to the output buffer */
    const q31_t *pTempSrc;             /* Temporary pointer to the source data */
    q31_t *pTempDest;                  /* Temporary pointer to the destination buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    int32_t blkCnt;
    q63_t acc0, acc1, acc2, acc3;

#define MAX_VECT_BATCH 7

    /*
     * pre-load 28 1st coefs
     */
    q31x4_t vecCoeffs0 = vld1q(pCoeffs + 4 * 0);
    q31x4_t vecCoeffs1 = vld1q(pCoeffs + 4 * 1);
    q31x4_t vecCoeffs2 = vld1q(pCoeffs + 4 * 2);
    q31x4_t vecCoeffs3 = vld1q(pCoeffs + 4 * 3);
    q31x4_t vecCoeffs4 = vld1q(pCoeffs + 4 * 4);
    q31x4_t vecCoeffs5 = vld1q(pCoeffs + 4 * 5);
    q31x4_t vecCoeffs6 = vld1q(pCoeffs + 4 * 6);

    /*
     * pState points to state array which contains previous frame (numTaps - 1) samples
     * pStateCur points to the location where the new input data should be written
     */
    pStateCur = &(pState[(numTaps - 1u)]);
    pTempSrc = pSrc;
    pSamples = pState;

    q63_t *arm_fir_partial_accu_ptr = (q63_t *)S->pState;

    blkCnt = blockSize >> 2;
    while (blkCnt > 0)
    {
        /*
         * Save 4 input samples in the history buffer
         */
        vstrwq_s32(pStateCur, vldrwq_s32(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        const q31_t *pSmp;
        q31x4_t vecIn0;

        pSmp = &pSamples[0];

        vecIn0 = vld1q(pSmp);
        acc0 = vrmlaldavhq(vecIn0, vecCoeffs0);
        vecIn0 = vld1q(pSmp + 4 * 1);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs1);
        vecIn0 = vld1q(pSmp + 4 * 2);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs2);
        vecIn0 = vld1q(pSmp + 4 * 3);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs3);
        vecIn0 = vld1q(pSmp + 4 * 4);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs4);
        vecIn0 = vld1q(pSmp + 4 * 5);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs5);
        vecIn0 = vld1q(pSmp + 4 * 6);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs6);

        *arm_fir_partial_accu_ptr++ = acc0;

        pSmp = &pSamples[1];

        vecIn0 = vld1q(pSmp);
        acc1 = vrmlaldavhq(vecIn0, vecCoeffs0);
        vecIn0 = vld1q(pSmp + 4 * 1);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs1);
        vecIn0 = vld1q(pSmp + 4 * 2);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs2);
        vecIn0 = vld1q(pSmp + 4 * 3);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs3);
        vecIn0 = vld1q(pSmp + 4 * 4);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs4);
        vecIn0 = vld1q(pSmp + 4 * 5);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs5);
        vecIn0 = vld1q(pSmp + 4 * 6);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs6);

        *arm_fir_partial_accu_ptr++ = acc1;

        pSmp = &pSamples[2];

        vecIn0 = vld1q(pSmp);
        acc2 = vrmlaldavhq(vecIn0, vecCoeffs0);
        vecIn0 = vld1q(pSmp + 4 * 1);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs1);
        vecIn0 = vld1q(pSmp + 4 * 2);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs2);
        vecIn0 = vld1q(pSmp + 4 * 3);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs3);
        vecIn0 = vld1q(pSmp + 4 * 4);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs4);
        vecIn0 = vld1q(pSmp + 4 * 5);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs5);
        vecIn0 = vld1q(pSmp + 4 * 6);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs6);
        *arm_fir_partial_accu_ptr++ = acc2;

        pSmp = &pSamples[3];

        vecIn0 = vld1q(pSmp);
        acc3 = vrmlaldavhq(vecIn0, vecCoeffs0);
        vecIn0 = vld1q(pSmp + 4 * 1);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs1);
        vecIn0 = vld1q(pSmp + 4 * 2);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs2);
        vecIn0 = vld1q(pSmp + 4 * 3);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs3);
        vecIn0 = vld1q(pSmp + 4 * 4);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs4);
        vecIn0 = vld1q(pSmp + 4 * 5);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs5);
        vecIn0 = vld1q(pSmp + 4 * 6);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs6);

        *arm_fir_partial_accu_ptr++ = acc3;

        pSamples += 4;
        /*
         * Decrement the sample block loop counter
         */
        blkCnt--;
    }

    /* reminder */

    /* load last 4 coef */
    vecCoeffs0 = vld1q(pCoeffs + 4 * MAX_VECT_BATCH);
    arm_fir_partial_accu_ptr = (q63_t *)S->pState;
    pOutput = pDst;
    pSamples = pState + (MAX_VECT_BATCH * 4);

    blkCnt = blockSize >> 2;
    while (blkCnt > 0)
    {
        q31x4_t vecIn0;

        /* reload intermediate MAC */
        acc0 = *arm_fir_partial_accu_ptr++;
        acc1 = *arm_fir_partial_accu_ptr++;
        acc2 = *arm_fir_partial_accu_ptr++;
        acc3 = *arm_fir_partial_accu_ptr++;

        vecIn0 = vld1q(&pSamples[0]);
        acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs0);

        vecIn0 = vld1q(&pSamples[1]);
        acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs0);

        vecIn0 = vld1q(&pSamples[2]);
        acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs0);

        vecIn0 = vld1q(&pSamples[3]);
        acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs0);

        *pOutput++ = asrl(acc0, 23);
        *pOutput++ = asrl(acc1, 23);
        *pOutput++ = asrl(acc2, 23);
        *pOutput++ = asrl(acc3, 23);

        pSamples += 4;
        /*
         * Decrement the sample block loop counter
         */
        blkCnt--;
    }

    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pState[blockSize];
    pTempDest = pState;

    blkCnt = numTaps - 1;
    do
    {
        mve_pred16_t p = vctp32q(blkCnt);

        vstrwq_p_s32(pTempDest, vldrwq_z_s32(pTempSrc, p), p);
        pTempSrc += 4;
        pTempDest += 4;
        blkCnt -= 4;
    } while (blkCnt > 0);
}

ARM_DSP_ATTRIBUTE void arm_fir_q31(
    const arm_fir_instance_q31 *S,
    const q31_t *pSrc,
    q31_t *pDst,
    uint32_t blockSize)
{
    q31_t *pRefStatePtr = S->pState + 2 * ARM_ROUND_UP(blockSize, 4);
    q31_t *pState = pRefStatePtr;      /* State pointer */
    const q31_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q31_t *pStateCur;                  /* Points to the current sample of the state */
    const q31_t *pSamples;             /* Temporary pointer to the sample buffer */
    q31_t *pOutput;                    /* Temporary pointer to the output buffer */
    const q31_t *pTempSrc;             /* Temporary pointer to the source data */
    q31_t *pTempDest;                  /* Temporary pointer to the destination buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t blkCnt;
    q31x4_t vecIn0;
    uint32_t tapsBlkCnt = (numTaps + 3) / 4;
    q63_t acc0, acc1, acc2, acc3;
    q31x4_t vecCoeffs;

    /*
     * [1 to 32 taps] specialized routines
     */
    if (numTaps <= 4)
    {
        arm_fir_q31_1_4_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 8)
    {
        arm_fir_q31_5_8_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 12)
    {
        arm_fir_q31_9_12_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 16)
    {
        arm_fir_q31_13_16_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 20)
    {
        arm_fir_q31_17_20_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 24)
    {
        arm_fir_q31_21_24_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if (numTaps <= 28)
    {
        arm_fir_q31_25_28_mve(S, pSrc, pDst, blockSize);
        return;
    }
    else if ((numTaps <= 32) && (blockSize >= 32))
    {
        arm_fir_q31_29_32_mve(S, pSrc, pDst, blockSize);
        return;
    }

    /*
     * pState points to state array which contains previous frame (numTaps - 1) samples
     * pStateCur points to the location where the new input data should be written
     */
    pStateCur = &(pState[(numTaps - 1u)]);
    pSamples = pState;
    pTempSrc = pSrc;
    pOutput = pDst;
    blkCnt = blockSize >> 2;
    while (blkCnt > 0)
    {
        const q31_t *pCoeffsTmp = pCoeffs;
        const q31_t *pSamplesTmp = pSamples;

        acc0 = 0LL;
        acc1 = 0LL;
        acc2 = 0LL;
        acc3 = 0LL;

        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            /*
             * load 4 coefs
             */
            vecCoeffs = *(q31x4_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[1]);
            acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[2]);
            acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[3]);
            acc3 = vrmlaldavhaq(acc3, vecIn0, vecCoeffs);

            pSamplesTmp += 4;
            pCoeffsTmp += 4;
            /*
             * Decrement the taps block loop counter
             */
            i--;
        }

        /* .54-> .31 conversion and store accumulators */
        acc0 = asrl(acc0, 23);
        acc1 = asrl(acc1, 23);
        acc2 = asrl(acc2, 23);
        acc3 = asrl(acc3, 23);

        *pOutput++ = (q31_t)acc0;
        *pOutput++ = (q31_t)acc1;
        *pOutput++ = (q31_t)acc2;
        *pOutput++ = (q31_t)acc3;

        pSamples += 4;

        /*
         * Decrement the sample block loop counter
         */
        blkCnt--;
    }

    int32_t residual = blockSize & 3;
    switch (residual)
    {
    case 3:
    {
        const q31_t *pCoeffsTmp = pCoeffs;
        const q31_t *pSamplesTmp = pSamples;

        acc0 = 0LL;
        acc1 = 0LL;
        acc2 = 0LL;

        /*
         * Save 4 input samples in the history buffer
         */
        *(q31x4_t *)pStateCur = *(q31x4_t *)pTempSrc;
        pStateCur += 4;
        pTempSrc += 4;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            vecCoeffs = *(q31x4_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[1]);
            acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[2]);
            acc2 = vrmlaldavhaq(acc2, vecIn0, vecCoeffs);

            pSamplesTmp += 4;
            pCoeffsTmp += 4;
            i--;
        }

        acc0 = asrl(acc0, 23);
        acc1 = asrl(acc1, 23);
        acc2 = asrl(acc2, 23);

        *pOutput++ = (q31_t)acc0;
        *pOutput++ = (q31_t)acc1;
        *pOutput++ = (q31_t)acc2;
    }
    break;

    case 2:
    {
        const q31_t *pCoeffsTmp = pCoeffs;
        const q31_t *pSamplesTmp = pSamples;

        acc0 = 0LL;
        acc1 = 0LL;

        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            vecCoeffs = *(q31x4_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs);

            vecIn0 = vld1q(&pSamplesTmp[1]);
            acc1 = vrmlaldavhaq(acc1, vecIn0, vecCoeffs);

            pSamplesTmp += 4;
            pCoeffsTmp += 4;
            i--;
        }

        acc0 = asrl(acc0, 23);
        acc1 = asrl(acc1, 23);

        *pOutput++ = (q31_t)acc0;
        *pOutput++ = (q31_t)acc1;
    }
    break;

    case 1:
    {
        const q31_t *pCoeffsTmp = pCoeffs;
        const q31_t *pSamplesTmp = pSamples;

        acc0 = 0LL;

        /*
         * Save 4 input samples in the history buffer
         */
        vst1q(pStateCur, vld1q(pTempSrc));
        pStateCur += 4;
        pTempSrc += 4;

        int i = tapsBlkCnt;
        while (i > 0)
        {
            vecCoeffs = *(q31x4_t *)pCoeffsTmp;

            vecIn0 = vld1q(pSamplesTmp);
            acc0 = vrmlaldavhaq(acc0, vecIn0, vecCoeffs);

            pSamplesTmp += 4;
            pCoeffsTmp += 4;
            i--;
        }

        acc0 = asrl(acc0, 23);

        *pOutput++ = (q31_t)acc0;
    }
    break;
    }

    /*
     * Copy the samples back into the history buffer start
     */
    pTempSrc = &pState[blockSize];
    pTempDest = pState;

    blkCnt = (numTaps - 1U) >> 2;
    while (blkCnt > 0)
    {
        vst1q(pTempDest, vld1q(pTempSrc));
        pTempSrc += 4;
        pTempDest += 4;
        blkCnt--;
    }
    blkCnt = (numTaps - 1U) & 3;
    if (blkCnt > 0)
    {
        mve_pred16_t p0 = vctp32q(blkCnt);
        vstrwq_p_s32(pTempDest, vld1q(pTempSrc), p0);
    }
}

#else

#if defined(ARM_MATH_NEON)

#define UPDATE_SUM_FIR_Q31_1(b, x0, x1)                      \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(x0), b[0]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(x0), b[0]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(x1), b[0]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(x1), b[0]);

#define UPDATE_SUM_FIR_Q31_2(b, x0, x1, x2)                  \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(x0), b[0]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(x0), b[0]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(x1), b[0]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(x1), b[0]); \
                                                             \
    xa = vextq_u32(x0, x1, 1);                               \
    xb = vextq_u32(x1, x2, 1);                               \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(xa), b[1]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(xa), b[1]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(xb), b[1]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(xb), b[1]);

#define UPDATE_SUM_FIR_Q31_3(b, x0, x1, x2)                  \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(x0), b[0]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(x0), b[0]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(x1), b[0]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(x1), b[0]); \
                                                             \
    xa = vextq_u32(x0, x1, 1);                               \
    xb = vextq_u32(x1, x2, 1);                               \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(xa), b[1]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(xa), b[1]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(xb), b[1]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(xb), b[1]); \
                                                             \
    xa = vextq_u32(x0, x1, 2);                               \
    xb = vextq_u32(x1, x2, 2);                               \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(xa), b[2]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(xa), b[2]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(xb), b[2]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(xb), b[2]);

#define UPDATE_SUM_FIR_Q31(b, x0, x1, x2)                    \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(x0), b[0]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(x0), b[0]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(x1), b[0]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(x1), b[0]); \
                                                             \
    xa = vextq_u32(x0, x1, 1);                               \
    xb = vextq_u32(x1, x2, 1);                               \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(xa), b[1]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(xa), b[1]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(xb), b[1]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(xb), b[1]); \
                                                             \
    xa = vextq_u32(x0, x1, 2);                               \
    xb = vextq_u32(x1, x2, 2);                               \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(xa), b[2]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(xa), b[2]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(xb), b[2]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(xb), b[2]); \
                                                             \
    xa = vextq_u32(x0, x1, 3);                               \
    xb = vextq_u32(x1, x2, 3);                               \
    accv0LO = vmlal_n_s32(accv0LO, vget_low_s32(xa), b[3]);  \
    accv0HI = vmlal_n_s32(accv0HI, vget_high_s32(xa), b[3]); \
    accv1LO = vmlal_n_s32(accv1LO, vget_low_s32(xb), b[3]);  \
    accv1HI = vmlal_n_s32(accv1HI, vget_high_s32(xb), b[3]);

#define START_OF_FUNCTION_COMMON()                                                            \
    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */ \
    /* pStateCurnt points to the location where the new input data should be written */       \
    pStateCurnt = &(S->pState[(numTaps - 1U)]);                                               \
                                                                                              \
    /* Initialize blkCnt with blockSize */                                                    \
    blkCnt = blockSize >> 3;                                                                  \
    while (blkCnt > 0U)                                                                       \
    {                                                                                         \
        /* Copy one sample at a time into state buffer */                                     \
        samples0 = vld1q_s32(pSrc);                                                           \
        vst1q_s32(pStateCurnt, samples0);                                                     \
                                                                                              \
        pStateCurnt += 4;                                                                     \
        pSrc += 4;                                                                            \
                                                                                              \
        samples1 = vld1q_s32(pSrc);                                                           \
        vst1q_s32(pStateCurnt, samples1);                                                     \
                                                                                              \
        pStateCurnt += 4;                                                                     \
        pSrc += 4;                                                                            \
        /* Set the accumulator to zero */                                                     \
        accv0LO = vdupq_n_s64(0);                                                             \
        accv0HI = vdupq_n_s64(0);                                                             \
        accv1LO = vdupq_n_s64(0);                                                             \
        accv1HI = vdupq_n_s64(0);                                                             \
                                                                                              \
        /* Initialize state pointer */                                                        \
        px = pState;                                                                          \
                                                                                              \
        /* Initialize Coefficient pointer */                                                  \
        pb = pCoeffs;

#define END_OF_FUNCTION_COMMON()                                                                                                         \
    tempLO = vshrn_n_s64(accv0LO, 31);                                                                                                   \
    tempHI = vshrn_n_s64(accv0HI, 31);                                                                                                   \
    accv0 = vcombine_s32(tempLO, tempHI);                                                                                                \
                                                                                                                                         \
    /* The result is stored in the destination buffer. */                                                                                \
    vst1q_s32(pDst, accv0);                                                                                                              \
    pDst += 4;                                                                                                                           \
                                                                                                                                         \
    tempLO = vshrn_n_s64(accv1LO, 31);                                                                                                   \
    tempHI = vshrn_n_s64(accv1HI, 31);                                                                                                   \
    accv1 = vcombine_s32(tempLO, tempHI);                                                                                                \
                                                                                                                                         \
    vst1q_s32(pDst, accv1);                                                                                                              \
    pDst += 4;                                                                                                                           \
                                                                                                                                         \
    /* Advance state pointer by 1 for the next sample */                                                                                 \
    pState = pState + 8;                                                                                                                 \
                                                                                                                                         \
    blkCnt--;                                                                                                                            \
    }                                                                                                                                    \
                                                                                                                                         \
    blkCnt = blockSize & 0x7;                                                                                                            \
    while (blkCnt > 0U)                                                                                                                  \
    {                                                                                                                                    \
        /* Copy one sample at a time into state buffer */                                                                                \
        *pStateCurnt++ = *pSrc++;                                                                                                        \
                                                                                                                                         \
        /* Set the accumulator to zero */                                                                                                \
        acc = 0;                                                                                                                         \
                                                                                                                                         \
        /* Initialize state pointer */                                                                                                   \
        px = pState;                                                                                                                     \
                                                                                                                                         \
        /* Initialize Coefficient pointer */                                                                                             \
        pb = pCoeffs;                                                                                                                    \
                                                                                                                                         \
        uint32_t i = numTaps;                                                                                                            \
                                                                                                                                         \
        /* Perform the multiply-accumulates */                                                                                           \
        do                                                                                                                               \
        {                                                                                                                                \
            /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
            acc += (q63_t) * px++ * *pb++;                                                                                               \
            i--;                                                                                                                         \
                                                                                                                                         \
        } while (i > 0U);                                                                                                                \
                                                                                                                                         \
        /* The result is stored in the destination buffer. */                                                                            \
        *pDst++ = (q31_t)(acc >> 31U);                                                                                                   \
                                                                                                                                         \
        /* Advance state pointer by 1 for the next sample */                                                                             \
        pState = pState + 1;                                                                                                             \
                                                                                                                                         \
        blkCnt--;                                                                                                                        \
    }                                                                                                                                    \
                                                                                                                                         \
    /* Processing is complete.                                                                                                           \
    ** Now copy the last numTaps - 1 samples to the starting of the state buffer.                                                        \
    ** This prepares the state buffer for the next function call. */                                                                     \
                                                                                                                                         \
    /* Points to the start of the state buffer */                                                                                        \
    pStateCurnt = S->pState;                                                                                                             \
                                                                                                                                         \
    /* Copy numTaps number of values */                                                                                                  \
    tapCnt = numTaps - 1U;                                                                                                               \
                                                                                                                                         \
    /* Copy data */                                                                                                                      \
    while (tapCnt > 0U)                                                                                                                  \
    {                                                                                                                                    \
        *pStateCurnt++ = *pState++;                                                                                                      \
                                                                                                                                         \
        /* Decrement the loop counter */                                                                                                 \
        tapCnt--;                                                                                                                        \
    }

#define INNER_LOOP_COEF_32_35()                                                                                                  \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 16);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[2], x2, x0, x1);                                                                                        \
    x2 = vld1q_s32(px + 20);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[3], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 24);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[4], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 28);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[5], x2, x0, x1);                                                                                        \
    x2 = vld1q_s32(px + 32);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[6], x0, x1, x2);                                                                                        \
                                                                                                                                 \
    x0 = vld1q_s32(px + 36);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[7], x1, x2, x0);                                                                                        \
                                                                                                                                 \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x1 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[8] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[8], x2, x0);                                                                                      \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x1 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[8] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[8], x2, x0, x1);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x1 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[8] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[8], x2, x0, x1);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_28_31()                                                                                                  \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 16);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[2], x2, x0, x1);                                                                                        \
    x2 = vld1q_s32(px + 20);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[3], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 24);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[4], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 28);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[5], x2, x0, x1);                                                                                        \
    x2 = vld1q_s32(px + 32);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[6], x0, x1, x2);                                                                                        \
                                                                                                                                 \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x0 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[7] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[7], x1, x2)                                                                                       \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x0 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[7] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[7], x1, x2, x0);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x0 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[7] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[7], x1, x2, x0);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_24_27()                                                                                                  \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 16);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[2], x2, x0, x1);                                                                                        \
                                                                                                                                 \
    x2 = vld1q_s32(px + 20);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[3], x0, x1, x2);                                                                                        \
                                                                                                                                 \
    x0 = vld1q_s32(px + 24);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[4], x1, x2, x0);                                                                                        \
                                                                                                                                 \
    x1 = vld1q_s32(px + 28);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[5], x2, x0, x1);                                                                                        \
                                                                                                                                 \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x2 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[6] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[6], x0, x1)                                                                                       \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x2 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[6] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[6], x0, x1, x2);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x2 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[6] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[6], x0, x1, x2);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_20_23()                                                                                                  \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 16);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[2], x2, x0, x1);                                                                                        \
                                                                                                                                 \
    x2 = vld1q_s32(px + 20);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[3], x0, x1, x2);                                                                                        \
                                                                                                                                 \
    x0 = vld1q_s32(px + 24);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[4], x1, x2, x0);                                                                                        \
                                                                                                                                 \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x1 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[5] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[5], x2, x0)                                                                                       \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x1 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[5] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[5], x2, x0, x1);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x1 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[5] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[5], x2, x0, x1);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_16_19()                                                                                                  \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 16);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[2], x2, x0, x1);                                                                                        \
                                                                                                                                 \
    x2 = vld1q_s32(px + 20);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[3], x0, x1, x2);                                                                                        \
                                                                                                                                 \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x0 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[4] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[4], x1, x2)                                                                                       \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x0 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[4] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[4], x1, x2, x0);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x0 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[4] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[4], x1, x2, x0);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_12_15()                                                                                                  \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
    x1 = vld1q_s32(px + 16);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[2], x2, x0, x1);                                                                                        \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x2 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[3] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[3], x0, x1)                                                                                       \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x2 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[3] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[3], x0, x1, x2);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x2 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[3] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[3], x0, x1, x2);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_8_11()                                                                                                   \
    x0 = vld1q_s32(px);                                                                                                          \
    x1 = vld1q_s32(px + 4);                                                                                                      \
    x2 = vld1q_s32(px + 8);                                                                                                      \
                                                                                                                                 \
    /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */ \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);                                                                                        \
                                                                                                                                 \
    x0 = vld1q_s32(px + 12);                                                                                                     \
    UPDATE_SUM_FIR_Q31(b[1], x1, x2, x0);                                                                                        \
                                                                                                                                 \
    switch (nTail)                                                                                                               \
    {                                                                                                                            \
    case 0:                                                                                                                      \
        break;                                                                                                                   \
    case 1:                                                                                                                      \
        x1 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[2] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_1(b[2], x2, x0)                                                                                       \
        break;                                                                                                                   \
    case 2:                                                                                                                      \
        x1 = vld1q_s32(px + (nTaps) * 4 + 8);                                                                                    \
        b[2] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_2(b[2], x2, x0, x1);                                                                                  \
        break;                                                                                                                   \
    case 3:                                                                                                                      \
        x1 = vld1q_s32(px + nTaps * 4 + 8);                                                                                      \
        b[2] = vld1q_s32(pb + (nTaps) * 4);                                                                                      \
        UPDATE_SUM_FIR_Q31_3(b[2], x2, x0, x1);                                                                                  \
        break;                                                                                                                   \
    default:                                                                                                                     \
        break;                                                                                                                   \
    }

#define INNER_LOOP_COEF_4_7()                   \
    x0 = vld1q_s32(px);                         \
    x1 = vld1q_s32(px + 4);                     \
    x2 = vld1q_s32(px + 8);                     \
                                                \
    UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);       \
                                                \
    switch (nTail)                              \
    {                                           \
    case 0:                                     \
        break;                                  \
    case 1:                                     \
        x0 = vld1q_s32(px + nTaps * 4 + 8);     \
        b[1] = vld1q_s32(pb + (nTaps) * 4);     \
        UPDATE_SUM_FIR_Q31_1(b[1], x1, x2)      \
        break;                                  \
    case 2:                                     \
        x0 = vld1q_s32(px + (nTaps) * 4 + 8);   \
        b[1] = vld1q_s32(pb + (nTaps) * 4);     \
        UPDATE_SUM_FIR_Q31_2(b[1], x1, x2, x0); \
        break;                                  \
    case 3:                                     \
        x0 = vld1q_s32(px + nTaps * 4 + 8);     \
        b[1] = vld1q_s32(pb + (nTaps) * 4);     \
        UPDATE_SUM_FIR_Q31_3(b[1], x1, x2, x0); \
        break;                                  \
    default:                                    \
        break;                                  \
    }

#define INNER_LOOP_GENERIC()                    \
                                                \
    b[0] = vld1q_s32(pb);                       \
    x0 = vld1q_s32(px);                         \
    x1 = vld1q_s32(px + 4);                     \
    x2 = vld1q_s32(px + 8);                     \
    for (int i = 0; i < nTaps; i++)             \
    {                                           \
        UPDATE_SUM_FIR_Q31(b[0], x0, x1, x2);   \
        b[0] = vld1q_s32(pb + (i + 1) * 4);     \
        x0 = vld1q_s32(px + (i + 1) * 4);       \
        x1 = vld1q_s32(px + (i + 1) * 4 + 4);   \
        x2 = vld1q_s32(px + (i + 1) * 4 + 8);   \
    }                                           \
    switch (nTail)                              \
    {                                           \
    case 0:                                     \
        break;                                  \
    case 1:                                     \
        UPDATE_SUM_FIR_Q31_1(b[0], x0, x1);     \
        break;                                  \
    case 2:                                     \
        UPDATE_SUM_FIR_Q31_2(b[0], x0, x1, x2); \
        break;                                  \
    case 3:                                     \
        UPDATE_SUM_FIR_Q31_3(b[0], x0, x1, x2); \
        break;                                  \
    default:                                    \
        break;                                  \
    }

ARM_DSP_ATTRIBUTE void arm_fir_q31(const arm_fir_instance_q31 *S,
                                   const q31_t *pSrc,
                                   q31_t *pDst,
                                   uint32_t blockSize)
{
    q31_t *pState = S->pState;         /* State pointer */
    const q31_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q31_t *pStateCurnt;                /* Points to the current sample of the state */
    q31_t *px;                         /* Temporary pointers for state buffer */
    const q31_t *pb;                   /* Temporary pointers for coefficient buffer */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t tapCnt, blkCnt;           /* Loop counters */
    int64x2_t accv0LO, accv1LO, accv0HI, accv1HI;
    int32x4_t accv0, accv1;
    int32x4_t samples0, samples1, x0, x1, x2, xa, xb;
    int32x2_t tempLO, tempHI;
    q63_t acc;
    int32_t nTaps = (numTaps) >> 2;
    int32_t nTail = numTaps & 0x3;
    switch (nTaps)
    {
    case 1:
    {
        pb = pCoeffs;
        int32x4_t b[2];
        b[0] = vld1q_s32(pb);
        if (nTail > 0)
            b[1] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_4_7();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 2:
    {
        pb = pCoeffs;
        int32x4_t b[3];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[2] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_8_11();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 3:
    {
        pb = pCoeffs;
        int32x4_t b[4];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[3] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_12_15();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 4:
    {
        pb = pCoeffs;
        int32x4_t b[5];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[4] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_16_19();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 5:
    {
        pb = pCoeffs;
        int32x4_t b[6];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[5] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_20_23();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 6:
    {
        pb = pCoeffs;
        int32x4_t b[7];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[6] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_24_27();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 7:
    {
        pb = pCoeffs;
        int32x4_t b[8];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[7] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_28_31();
        END_OF_FUNCTION_COMMON();
        return;
    }
    case 8:
    {
        pb = pCoeffs;
        int32x4_t b[9];

        for (int i = 0; i < nTaps; i++)
            b[i] = vld1q_s32(pb + i * 4);
        if (nTail > 0)
            b[8] = vld1q_s32(pb + 4);
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_COEF_32_35();
        END_OF_FUNCTION_COMMON();
        return;
    }
    default:
    {
        int32x4_t b[1];
        START_OF_FUNCTION_COMMON();
        INNER_LOOP_GENERIC();
        END_OF_FUNCTION_COMMON();
        return;
    }
    }
}

#else

ARM_DSP_ATTRIBUTE void arm_fir_q31(
    const arm_fir_instance_q31 *S,
    const q31_t *pSrc,
    q31_t *pDst,
    uint32_t blockSize)
{
    q31_t *pState = S->pState;         /* State pointer */
    const q31_t *pCoeffs = S->pCoeffs; /* Coefficient pointer */
    q31_t *pStateCurnt;                /* Points to the current sample of the state */
    q31_t *px;                         /* Temporary pointer for state buffer */
    const q31_t *pb;                   /* Temporary pointer for coefficient buffer */
    q63_t acc0;                        /* Accumulator */
    uint32_t numTaps = S->numTaps;     /* Number of filter coefficients in the filter */
    uint32_t i, tapCnt, blkCnt;        /* Loop counters */

#if defined(ARM_MATH_LOOPUNROLL)
    q63_t acc1, acc2; /* Accumulators */
    q31_t x0, x1, x2; /* Temporary variables to hold state values */
    q31_t c0;         /* Temporary variable to hold coefficient value */
#endif

    /* S->pState points to state array which contains previous frame (numTaps - 1) samples */
    /* pStateCurnt points to the location where the new input data should be written */
    pStateCurnt = &(S->pState[(numTaps - 1U)]);

#if defined(ARM_MATH_LOOPUNROLL)

    /* Loop unrolling: Compute 4 output values simultaneously.
     * The variables acc0 ... acc3 hold output values that are being computed:
     *
     *    acc0 =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0]
     *    acc1 =  b[numTaps-1] * x[n-numTaps]   + b[numTaps-2] * x[n-numTaps-1] + b[numTaps-3] * x[n-numTaps-2] +...+ b[0] * x[1]
     *    acc2 =  b[numTaps-1] * x[n-numTaps+1] + b[numTaps-2] * x[n-numTaps]   + b[numTaps-3] * x[n-numTaps-1] +...+ b[0] * x[2]
     *    acc3 =  b[numTaps-1] * x[n-numTaps+2] + b[numTaps-2] * x[n-numTaps+1] + b[numTaps-3] * x[n-numTaps]   +...+ b[0] * x[3]
     */

    blkCnt = blockSize / 3;

    while (blkCnt > 0U)
    {
        /* Copy 3 new input samples into the state buffer. */
        *pStateCurnt++ = *pSrc++;
        *pStateCurnt++ = *pSrc++;
        *pStateCurnt++ = *pSrc++;

        /* Set all accumulators to zero */
        acc0 = 0;
        acc1 = 0;
        acc2 = 0;

        /* Initialize state pointer */
        px = pState;

        /* Initialize coefficient pointer */
        pb = pCoeffs;

        /* Read the first 2 samples from the state buffer: x[n-numTaps], x[n-numTaps-1] */
        x0 = *px++;
        x1 = *px++;

        /* Loop unrolling: process 3 taps at a time. */
        tapCnt = numTaps / 3;

        while (tapCnt > 0U)
        {
            /* Read the b[numTaps] coefficient */
            c0 = *pb;

            /* Read x[n-numTaps-2] sample */
            x2 = *(px++);

            /* Perform the multiply-accumulates */
            acc0 += ((q63_t)x0 * c0);
            acc1 += ((q63_t)x1 * c0);
            acc2 += ((q63_t)x2 * c0);

            /* Read the coefficient and state */
            c0 = *(pb + 1U);
            x0 = *(px++);

            /* Perform the multiply-accumulates */
            acc0 += ((q63_t)x1 * c0);
            acc1 += ((q63_t)x2 * c0);
            acc2 += ((q63_t)x0 * c0);

            /* Read the coefficient and state */
            c0 = *(pb + 2U);
            x1 = *(px++);

            /* update coefficient pointer */
            pb += 3U;

            /* Perform the multiply-accumulates */
            acc0 += ((q63_t)x2 * c0);
            acc1 += ((q63_t)x0 * c0);
            acc2 += ((q63_t)x1 * c0);

            /* Decrement loop counter */
            tapCnt--;
        }

        /* Loop unrolling: Compute remaining outputs */
        tapCnt = numTaps % 0x3U;

        while (tapCnt > 0U)
        {
            /* Read coefficients */
            c0 = *(pb++);

            /* Fetch 1 state variable */
            x2 = *(px++);

            /* Perform the multiply-accumulates */
            acc0 += ((q63_t)x0 * c0);
            acc1 += ((q63_t)x1 * c0);
            acc2 += ((q63_t)x2 * c0);

            /* Reuse the present sample states for next sample */
            x0 = x1;
            x1 = x2;

            /* Decrement loop counter */
            tapCnt--;
        }

        /* Advance the state pointer by 3 to process the next group of 3 samples */
        pState = pState + 3;

        /* The result is in 2.30 format. Convert to 1.31 and store in destination buffer. */
        *pDst++ = (q31_t)(acc0 >> 31U);
        *pDst++ = (q31_t)(acc1 >> 31U);
        *pDst++ = (q31_t)(acc2 >> 31U);

        /* Decrement loop counter */
        blkCnt--;
    }

    /* Loop unrolling: Compute remaining output samples */
    blkCnt = blockSize % 0x3U;

#else

    /* Initialize blkCnt with number of taps */
    blkCnt = blockSize;

#endif /* #if defined (ARM_MATH_LOOPUNROLL) */

    while (blkCnt > 0U)
    {
        /* Copy one sample at a time into state buffer */
        *pStateCurnt++ = *pSrc++;

        /* Set the accumulator to zero */
        acc0 = 0;

        /* Initialize state pointer */
        px = pState;

        /* Initialize Coefficient pointer */
        pb = pCoeffs;

        i = numTaps;

        /* Perform the multiply-accumulates */
        do
        {
            /* acc =  b[numTaps-1] * x[n-numTaps-1] + b[numTaps-2] * x[n-numTaps-2] + b[numTaps-3] * x[n-numTaps-3] +...+ b[0] * x[0] */
            acc0 += (q63_t)*px++ * *pb++;

            i--;
        } while (i > 0U);

        /* Result is in 2.62 format. Convert to 1.31 and store in destination buffer. */
        *pDst++ = (q31_t)(acc0 >> 31U);

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

#if defined(ARM_MATH_LOOPUNROLL)

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
#endif /* defined(ARM_MATH_MVEI) */
#endif /* defined(ARM_NEON) */
/**
  @} end of FIR group
 */
